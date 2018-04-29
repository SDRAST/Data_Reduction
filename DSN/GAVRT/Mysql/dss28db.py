# -*- coding: utf-8 -*-
"""
Support for GAVRT mysql database

Example::
 

The databases and their schemas are described in
http://gsc.lewiscenter.org/data_info/dss28_eac.php.

The server has these databases::
 'dss28_eac'
 'dss28_spec'
 'gavrt_sources'.

Database 'dss28_eac' has these tables::
  In [17]: dbplotter.get_public_tables()
  Out[17]: 
  (('angles',),
   ('chan_cfg',),
   ('conv_cfg',),
   ('fiber_cfg',),
   ('five_point',),
   ('pointing_cfg',),
   ('raster',),
   ('raster_cfg',),
   ('rf_cfg',),
   ('rss_cfg',),
   ('seti_cfg',),
   ('seti_frame',),
   ('tlog',),
   ('weather',),
   ('xpwr',),
   ('xpwr_cfg',),
   ('xscan',),
   ('zplot',),
   ('zplot_cfg',))

Database 'gavrt_sources' has these tables::
 'catalog',
 'class',
 'source'

'raster' columns::
  raster_id,
  raster_cfg_id, year, doy, utc, epoch, xdecoff, decoff, ha, dec, tsrc
'rss_cfg' columns::
  rss_cfg_id,
  year, doy, utc, chan, sky_freq, feed, pol, nd, if_mode, if_bw, bb_bw, fiber_chan
'source' columns::
  source_id, catalog_id, class_id,
  name, RA, Dec, size_dec, size_xdec, reference, aka
'weather' columns::
  weather_id,
  datetime, pressure, temp, humidity, wind_speed, wind_dir
'xpwr' columns::
  xpwr_id, xpwr_cfg_id,
  year, doy, utc, epoch, tsys, az, el, ha, dec, offset
'xpwr_cfg' columns::
  xpwr_cfg_id, source_id,
  axis, chan
'xscan' columns::
  xscan_id, xpwr_cfg_id,
  year, doy, utc, epoch, tsrc, stdev, bl_stdev, az, az_offset, el, el_offset,
  ha, dec, offset, bw, corr
"""
import cPickle as pickle
import datetime
import ephem
import logging
import MySQLdb
import os
import stat
import time

from math import cos, pi
from matplotlib.mlab import griddata
from matplotlib.pylab import date2num, num2date
from numpy import arange

from Data_Reduction.DSN.GAVRT.Mysql import BaseDB, MysqlException
from support.lists import unique

logger = logging.getLogger(__name__)

import Astronomy as A

from MonitorControl.Configurations.coordinates import DSS
from Radio_Astronomy import HPBW

_host,_user,_pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p",
                                     "rb" ))
dss28 = DSS(28)
longitude = -dss28.long*180/pi
latitude = dss28.lat*180/pi
f_max = 16. # GHz
wl_min = f_max/300
taper = 12 # dB
hpbw = HPBW(taper, wl_min, 34)*180/pi # deg
default_step = hpbw/3.

class Map(object):
  """
  Class for all the data and methods associated with a raster scan map
  
  Public attributes::
    cfg         - raster configuration
    cfg_id      - entry in the raster configuration table
    channels    - list of channels which took tlog data
    map_data    - dict of data from tlog table; 'tsrc' is dict keyed on channel
    logger      - logging.Logger object
    name        - map identifier
    raster_data - data from the raster table
    regrid      - computes map data onto a rectangular grid
    rss_cfg     - receiver configuration
    session     - observing session to which this map belongs
    start       - UNIX time at start of map
    end         - UNIX time at end of map
    
  Public methods::
    center_map          - converts map coordinates to be relative to Sun
    get_active_channels - returns channels which took tlog data during this map
    get_map_config      - returns a dict with the raster map configuration
    get_raster_data     - gets the data for a raster scan map used for Zplot
    get_raster_keys     - returns rasters associated with a given configuration
    maps_from_tlogs     - re-organizes tlog table data into map form
  """
  def __init__(self, parent, raster_cfg_id, name=None):
    """
    initialize a Map object
    """
    self.logger = logging.getLogger(parent.logger.name+".Map")
    self.session = parent
    self.cfg_id = raster_cfg_id
    if name:
      self.name = name
    else:
      self.name = "map %d" % self.cfg_id
    self.logger.debug("__init__: map for configuration %s is %s",
                      self.session, self.name)
    # this applies to the default channel.
    self.cfg = self.get_map_config()
    # this is used by get_raster_data
    self.get_raster_keys()
    # this defines 'start' and 'end'
    self.get_raster_data()
    # gets from 'tlog' the channels used between 'start' and 'end'
    self.get_active_channels()
         
  def get_map_config(self):
    """
    returns a dict with the raster map configuration
    
    This is the map generated for display during the observing session.  The
    observer selects a channel.  Other channel data are recorded in the t-logs
    for the channels enabled.
    
    Example::
      In [9]: map56.get_map_config()
      Out[9]: 
      {'chan': 4.0,
       'doy': 233.0,
       'epoch': 1503329279.0,
       'freq': 8000.0,
       'raster_cfg_id': 56.0,
       'rate': 0.035000000000000003,
       'rss_cfg_id': 109708.0,
       'source_id': 66.0,
       'step': 0.035000000000000003,
       'utc': datetime.timedelta(0, 55679),
       'year': 2017.0}

    """
    cfg = self.session.db.get_as_dict(
        "select * from raster_cfg where raster_cfg_id = " + \
                                                         str(self.cfg_id) +";")
    self.cfg = {}
    for key in cfg.keys():
      self.cfg[key] = cfg[key][0]
    return self.cfg
        
  def get_raster_keys(self):
    """
    Returns the rasters associated with a given configuration
    """
    rasterkeys = self.session.db.get(
                        "select raster_id from raster where raster_cfg_id = " + 
                                                          str(self.cfg_id)+";")
    rasterkeys.sort()
    self.raster_keys = rasterkeys[:,0]
    return self.raster_keys
  
  def get_raster_data(self):
    """
    gets the data for a raster scan map extracted for Zplot
    """
    data = self.session.db.get(
          "select epoch, xdecoff,decoff,tsrc from raster where raster_id >= " +
          str(self.raster_keys[0]) + " and raster_id <= " + 
          str(self.raster_keys[-1]) + ";")
    self.logger.debug("get_raster_data: data shape is %s", data.shape)
    self.raster_data = {}
    self.raster_data['UNIXtime'] = data[:,0].astype(float)
    self.raster_data['xdec']     = data[:,1].astype(float)
    self.raster_data['dec']      = data[:,2].astype(float)
    self.raster_data['tsrc']     = data[:,3].astype(float)
    self.start = self.raster_data['UNIXtime'][0]
    self.end   = self.raster_data['UNIXtime'][-1]
    return self.raster_data

  def get_active_channels(self):
    """
    returns the receiver channels that were active between 'start' and 'end'
    
    Example::
      In [5]: map56.get_active_channels()
      Out[5]: [2, 4]
    """
    response = self.session.db.get("select chan from tlog where epoch >= " +
                       str(self.start) + " and epoch <=" + str(self.end) + ";")
    self.channels = unique(list(response[:].flatten()))
    self.rss_cfg = {}
    for chan in self.channels:
      response = self.session.db.get(
                       "select rss_cfg_id from tlog where chan = " +str(chan)
                       + " and epoch >= " + str(self.start) + " and epoch <="
                       + str(self.end) + ";")
      rss_cfg_id = unique(list(response[:,0].flatten()))[0]
      self.rss_cfg[chan] = self.session.db.get_as_dict(
                                  "select * from rss_cfg where rss_cfg_id = " +
                                  str(rss_cfg_id) + ";")
    return self.channels
    
  def maps_from_tlogs(self):
    """
    Gets the data for the specified channel and polarization for this map
    
    NOTE: It occurs to me 
    """
    try:
      chan_list = self.channels
    except:
      self.get_active_channels()
    self.map_data = {}
    if self.channels:
      pass
    else:
      self.logger.warning("maps_from_tlogs: this map has no active channels")
      return None
    for channel in self.channels:
      ch_index = self.channels.index(channel)
      query = "select tlog.epoch, tlog.az, tlog.el, tlog.top, tlog.rss_cfg_id" +               \
            " from tlog, rss_cfg where tlog.rss_cfg_id = rss_cfg.rss_cfg_id" +\
            " and tlog.epoch >= "  + str(self.start) +                        \
            " and tlog.epoch <= "  + str(self.end  ) +                        \
            " and rss_cfg.chan = " + str(channel)      +";"
      self.logger.debug("maps_from_tlogs: query: %s", query)
      data = self.session.db.get(query)
      # get receiver configuration
      self.rss_cfg[channel] = self.session.db.get_as_dict(
                                  "select * from rss_cfg where rss_cfg_id = " +
                                  str(int(data[:,4][0])) + ";")
      if ch_index == 0:
        # first channel only: these are common to all channels
        self.map_data['UNIXtime']    = data[:,0].astype(float)
        self.map_data['azimuth']     = data[:,1].astype(float)
        self.map_data['elevation']   = data[:,2].astype(float)
        self.map_data['RA']          = []
        self.map_data['declination'] = []
        self.map_data['MPL_datenum'] = []
        self.map_data['VFC_counts']  = {}
        self.map_data['freq']        = {}
        self.map_data['pol']         = {}
        for index in range(len(self.map_data['UNIXtime'])):
          dt = datetime.datetime.utcfromtimestamp(
                                     self.map_data['UNIXtime'][index])
          time_tuple = (dt.year,
                        A.day_of_year(dt.year,dt.month,dt.day)
                        + (  dt.hour
                           + dt.minute/60.
                           + dt.second/3600.
                           + dt.microsecond/3600./1e6)/24.)
          ra, dec = A.AzEl_to_RaDec(
                             float(self.map_data['azimuth'][index]),
                             float(self.map_data['elevation'][index]),
                             latitude,
                             longitude,
                             time_tuple)
          self.map_data['RA'].append(ra)
          self.map_data['declination'].append(dec)
          self.map_data['MPL_datenum'].append(date2num(dt))
      # only the VFC counts differ between channels
      self.map_data['VFC_counts'][channel]  = data[:,3].astype(float)
      self.map_data['freq'][channel] = self.rss_cfg[channel]['sky_freq']
      self.map_data['pol'][channel] = self.rss_cfg[channel]['pol'][0].upper()
    return self.map_data

  def center_map(self, xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source

    @param body : source at map center
    @type  body : ephem source instance

    @return: (dxdecs,ddecs) in degrees
    """
    try:
      self.map_data.keys()
    except AttributeError:
      self.maps_from_tlogs()
    if self.map_data.has_key('MPL_datenum'):
      pass
    else:
      # 'tlog' did not have good data
      return None
    sun = ephem.Sun()
    self.map_data['dec_offset'] = []
    self.map_data['xdec_offset'] = []
    for count in range(len(self.map_data['MPL_datenum'])):
      dt = num2date(self.map_data['MPL_datenum'][count])
      sun.compute(dt)
      ra_center = sun.ra*12/pi    # hours
      dec_center = sun.dec*180/pi # degrees
      decrad = sun.dec
      # right ascension increases to the left, cross-dec to the right
      self.map_data['xdec_offset'].append(xdec_ofst - 
                       (self.map_data['RA'][count] - ra_center)*15*cos(decrad) )
      self.map_data['dec_offset'].append(  dec_ofst + 
                        self.map_data['declination'][count] - dec_center)
    return self.map_data['xdec_offset'], self.map_data['dec_offset']
  
  def regrid(self, width=1.0, height=1.0, step=None):
    """
    converts a map from observed coordinates to map coordinates
    
    @param width : map width in deg
    @type  width : float
    
    @param height : map height in deg
    @type  height : float
    
    @param step : map step size in deg
    @type  step : float
    """
    if step == None:
      # use the original step size
      step = self.cfg['step']
    nx = int(round(width/step))
    ny = int(round(height/step))
    self.map_data['grid_x'] = arange( -width/2,  width/2 + step, step/2)
    self.map_data['grid_y'] = arange(-height/2, height/2 + step, step/2)
    self.map_data['grid_z'] = {}
    for channel in self.channels:
      self.map_data['grid_z'][channel] = griddata(self.map_data['xdec_offset'],
                                                  self.map_data['dec_offset'],
                                          self.map_data['VFC_counts'][channel],
                                          self.map_data['grid_x'], 
                                          self.map_data['grid_y'])
    return self.map_data['grid_x'], self.map_data['grid_y'], \
           self.map_data['grid_z']


class Session(object):
  """
  Class for an observing session on a given year and DOY
  """
  def __init__(self, parent, year, doy, plotter=False):
    """
    """
    self.logger = logging.getLogger(parent.logger.name+".Session")
    self.db = parent
    self.year = year
    self.doy = doy
    if plotter == False:
      self.maps = self.get_maps()
   
  def get_maps(self):
    """
    Returns the raster configuration IDs for the specified date
    """
    map_cfg_ids = self.db.get(
                         "select raster_cfg_id from raster_cfg where year = " +
                         str(self.year) + " and doy = " + str(self.doy) + 
                         ";")
    self.logger.debug("get_maps: map IDs: %s", map_cfg_ids)
    self.maps = {}
    for map_id in map_cfg_ids[:,0]:
      self.logger.debug("get_maps: getting %d", map_id)
      self.maps[map_id] = Map(self, map_id)
    return self.maps
    
  def make_directory(self):
    """
    show maps in the database for this session
    """
    map_keys = self.maps.keys()
    map_keys.sort()
    response = ["Map YYYY/DDD -start-- --end--- chan -freq pol"]
    txt_format = "%3d %4d/%03d %02d:%02d:%02d %02d:%02d:%02d  %2d  %5d %s"
    for key in map_keys:
      start = time.gmtime(self.maps[key].start)
      end   = time.gmtime(self.maps[key].end)
      for channel in self.maps[key].channels:
        freq = self.maps[key].rss_cfg[channel]['sky_freq']
        pol = self.maps[key].rss_cfg[channel]['pol']
        txt = txt_format % (key, start.tm_year, start.tm_yday,
                            start.tm_hour, start.tm_min, start.tm_sec,
                            end.tm_hour, end.tm_min, end.tm_sec,
                            channel, freq, pol)
        response.append(txt)
    return response

  def save_map_data(self, mapkeys=None):
    """
    create a dict with the map data from the designated images
    
    @param mapkeys : numbers of the maps (default: all)
    @type  mapkeys : list of int
    """
    if mapkeys:
      self.logger.info("show_images:")
    else:
      mapkeys = self.maps.keys()
      mapkeys.sort()
    for key in mapkeys:
      try:
        self.maps[key].map_data.keys()
        self.logger.debug("save_map_data: mapdata[%d] exists", key)
      except AttributeError:
        self.maps[key].maps_from_tlogs()
        self.logger.debug("save_map_data: loaded mapdata[%d]", key)
      if self.maps[key].map_data.has_key('dec_offset'):
        self.logger.debug("save_map_data: mapdata[%d] is centered", key)
      else:
        self.maps[key].center_map()
        self.logger.debug("save_map_data: mapdata[%d] has been centered", key)
      if self.maps[key].map_data.has_key('grid_x'):
        self.logger.debug("save_map_data: mapdata[%d] is regridded", key)
      else:
        self.maps[key].regrid()
        self.logger.debug("save_map_data: mapdata[%d] has been regridded", key)
    export = {}
    for key in mapkeys:
      export[key] = self.maps[key].map_data
    filename = "maps-%4d-%03d.pkl" % (self.year, self.doy)
    exportfile = open(filename, "w")
    pickle.dump(export, exportfile)
    exportfile.close()
    return export
    
    

class DSS28db(BaseDB):
  """
  subclass for the DSS-28 EAC database
  
  provides methods for handling tables
  
  Attributes::
    logger   - logging.Logger object
    receiver - receivers which provide data
    sessions - dict of sessions obtained with 'get_session'
  """
  def __init__(self, host=_host, user=_user, pw=_pw,
                     name='dss28_eac', port=3306):
    """
    create an instance BaseDB subclass for the DSS-28 EAC database
    
    The defaults for BaseDB are for the DSS-28 EAC database
    """
    mylogger = logging.getLogger(logger.name+".DSS28db")
    BaseDB.__init__(self, host=host, user=user, pw=pw, name=name, port=port)
    self.logger = mylogger
    self.sessions = {}

  def insertRecord(self, table, rec):
    """
    not allowed for subclass
    """
    self.logger.warning("insertRecord: not allowed for %s", self.name)

  def updateValues(self, vald, table): 
    """
    not allowed for subclass
    """
    self.logger.warning("updateValues: not allowed for %s", self.name)
        
  def extract_boresight_data(self, year, doy):
    """
    Get the metadata for the boresights on the designated day.

    The boresights are extracted from table 'xscan'.  Missing 'el' data are
    obtained from table 'xpwr'.  The source, scan axis and channel are obtained
    from table 'xpwr_cfg'.  The receiver data are obtained from table 'rss_cfg'.

    Returns a dictionary like this::
      {'utc':        list of datetime.timedelta,
       'epoch':      list of float,
       'az':         list of float,
       'el':         list of value,
       'chan':       list of int,
       'tsrc':       list of float,
       'axis':       list of str,
       'source':     list of str,
       'xpwr_cfg_id: list of int',
       'xscan_id':   list of int,
       'source_id':  list of int,
       'rx':         list of dict}
     
    An 'rx' dict looks like this::
      { 2: {'if_bw':    float,
            'if_mode':  str,
            'pol':      str,
            'sky_freq': float,
            'utc':      datetime.timedelta},
        4: { ... },
       ....
       16: { ... }}

    @param year : year of observation
    @type  year : int

    @param doy : day of year
    @type  doy : int

    @return: dict
    """
    # Get the boresight data from xscan
    columns = "utc, epoch, tsrc, az, el, xscan_id, xpwr_cfg_id"
    boresight_data = self.get_rows_by_date("xscan", columns, year, doy)

    # Get the missing elevation data from xpwr
    times = boresight_data['utc']
    power_data = self.get_rows_by_time('xpwr',['utc','el','tsys'],
                                     year,doy,times)
    # Fix the missing elevation data
    boresight_data['el'] = power_data['el']

    # Get the source information from gavrt_sources.source
    columns = "source_id, axis, chan"
    for column in columns.split(','):
      boresight_data[column.strip()] = []
    for cfg_id in boresight_data['xpwr_cfg_id']:
      response = self.get_as_dict("select "
                        + columns
                        + " from xpwr_cfg where xpwr_cfg_id="+str(cfg_id)+";")
    for key in response.keys():
      boresight_data[key].append(response[key][0])
    boresight_data['source'] = []
    for source_id in boresight_data['source_id']:
      response = self.get_as_dict("select name from gavrt_sources.source where source_id="
                        +str(source_id)+";")
      boresight_data['source'].append(response['name'][0])

    # Get the receiver information from rss_cfg
    columns = "utc,sky_freq,pol,if_mode,if_bw"
    boresight_data['rx'] = []
    for time in times:
      boresight_data['rx'].append(self.get_receiver_data(year,doy,time,columns))

    return boresight_data

  def get_receiver_data(self, year, doy, time, columns):
    """
    Get the receiver state at a given time

    This creates a dictionary keyed with channel number and returns a dictionary
    of the receiver configuration, keyed with specified in the columns, that was
    in effect at the given time.

    Notes
    =====
    Logic
    -----
    The challenge here is to get the latest configuration data for each channel
    at or prior to the specified time.  That channel may have been configured on
    the same day or a prior day. The method we'll use is to find the ID of last
    configuration change and assume that the IDs are sequential in date/time.

    @param db : database
    @type  db : Mysql.BaseDB instance

    @param year : year of observation
    @type  year : int

    @param doy : day of year
    @type  doy : int

    @param time : UTC for the requested receiver state
    @type  time : datetime.timedelta

    @param columns : data items to be returned
    @type  columns : list of str

    @return: dict
    """
    columns = columns.replace(" ","")
    column_keys = columns.split(',')
    self.logger.debug("get_receiver_data: requesting %s", column_keys)
    latest_data = self.get_as_dict("select rss_cfg_id,year,doy,utc from rss_cfg"
                        +" where year <= "+str(year)
                        +" and doy <= "+str(doy)
                        +" and utc <= '"+str(time)
                       +"' order by year desc, doy desc, utc desc limit 1;")
    self.logger.debug("get_receiver_data: query response: %s",latest_data)
    cfg_ID = latest_data['rss_cfg_id'][0]
    self.receiver = {}
    for key in column_keys:
      self.receiver[key] = {}
      for chan in [2,4,6,8,10,12,14,16]:
        rx_data = self.get_as_dict("select "+columns
                       +" from rss_cfg where rss_cfg_id <= "+str(cfg_ID)
                       +" and chan = "+str(chan)
                       +" order by rss_cfg_id desc limit 1;")
        index = column_keys.index(key)
        self.receiver[key][chan] = rx_data[key][0]
    return self.receiver

  def get_Tsys(self, chan, start, stop):
    """
    Get system temperatures from tlog
    
    @param start : UNIXtime at start of selection
    @type  start : float
    
    @param stop : UNIXtime at end of selection
    @type  stop :float
    """
    query = \
        'select epoch, top from tlog where chan = %d and epoch >= %f and epoch <= %f' \
                                                                % (chan, start, stop)
    try:
      response = self.get_as_dict(query)
      return response
    except Exception, details:
      self.logger.error("get_Tsys: error: %s", str(details))
      return None
 
  def get_session(self, year, doy):
    """
    get IDs for an observing session
    """
    if self.sessions.has_key(year) == False:
      self.sessions[year] = {}
    self.sessions[year][doy] = Session(self, year, doy)
    return self.sessions[year][doy]

