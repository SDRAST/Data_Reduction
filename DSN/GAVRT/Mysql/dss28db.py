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
import numpy as NP
import os
import stat
import sys
import time

from math import cos, pi
from matplotlib.mlab import griddata
from matplotlib.pylab import date2num, num2date
from numpy import arange, unique, where
from time import gmtime, strftime

from local_dirs import projects_dir
import Astronomy
from Data_Reduction.DSN.GAVRT.Mysql import BaseDB, MysqlException
#from support.lists import unique

logger = logging.getLogger(__name__)

import Astronomy as A

from Astronomy.DSN_coordinates import DSS
from Radio_Astronomy import HPBW

_host,_user,_pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p", "rb" ))
dss28 = DSS(28)
longitude = -dss28.long*180/pi
latitude = dss28.lat*180/pi
f_max = 16. # GHz
wl_min = f_max/300
taper = 12 # dB
hpbw = HPBW(taper, wl_min, 34)*180/pi # deg
default_step = hpbw/3.

class Observation(object):
  """
  Class for any group of data for a single purpose.
  
  Attributes::
    parent - a collection or group of observations
  """
  def __init__(self, parent):
    """
    """
    self.logger = logging.getLogger(parent.logger.name+".Observation")
    self.parent = parent

  def get_conv_config(self, time, converter):
    """
    get last configuration change for this converter
    """
    # search the last 11.6 days
    self.conv_cfg = self.session.db.get_as_dict(
       "select conv_cfg_id,converter,lock_status,atten_a,atten_b from conv_cfg"
                       +" where epoch <= "+str(time)
                       +" and epoch >= "+str(float(time)-1e6)
                       +" and converter = "+str(converter)
                       +" order by epoch desc limit 1;")
  
  def get_channel_attenuation(self, time, channel):
    """
    get last configuration change for this channel
    """
    converter = 1+(channel-1)//4
    self.get_conv_config(time, converter)
    conv_chl = ((channel-1)%4)//2 # 0 means 'a'; 1 means 'b'
    side = chr(ord('a')+conv_chl)
    self.logger.debug("get_channel_attenuation: channel %d is converter %d%s",
                      channel, converter, side)
    attenuator = "atten_"+side
    self.logger.debug("get_channel_attenuation: using %s", attenuator)
    self.conv_cfg['atten'] = self.conv_cfg[attenuator]
    

class Map(Observation):
  """
  Class for all the data and methods associated with a raster scan map
  
  Public attributes::
    cfg         - raster configuration
    cfg_id      - entry in the raster configuration tableshape
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
    if self.raster_keys != None:
      # this defines 'start' and 'end'
      self.get_raster_data()
      # gets from 'tlog' the channels used between 'start' and 'end'
      self.get_active_channels()
    self.source = self.session.db.get_source_names([self.cfg['source_id']])['source'][0]
         
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
         "select rss_cfg_id, raster_cfg_id, step, source_id from raster_cfg " +
                              "where raster_cfg_id = " + str(self.cfg_id) +";")

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
    if rasterkeys.any():
      rasterkeys.sort()
      self.raster_keys = rasterkeys[:,0]
    else:
      self.raster_keys = None
    return self.raster_keys
    
      
  def get_raster_data(self):
    """
    gets the data for a raster scan map extracted for Zplot
    """
    data = self.session.db.get(
          "select epoch, xdecoff,decoff,tsrc from raster where raster_id >= " +
          str(self.raster_keys[0]) + " and raster_id <= " + 
          str(self.raster_keys[-1]) + ";")
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
    # find out which channels were active
    response = self.session.db.get("select chan from tlog where epoch >= " +
                       str(self.start) + " and epoch <=" + str(self.end) + ";")
    self.channels = NP.unique(response[:].flatten()) # unique(list(response[:].flatten()))
    # for the active channels get the rss_cfg data
    self.rss_cfg = {}
    for chan in self.channels:
      # get the RSS configuration for that channel
      response = self.session.db.get(
                       "select rss_cfg_id from tlog where chan = " +str(chan)
                       + " and epoch >= " + str(self.start) + " and epoch <="
                       + str(self.end) + ";")
      rss_cfg_id = NP.unique(response[:].flatten())[0] # unique(list(response[:,0].flatten()))[0]
      self.rss_cfg[chan] = self.session.db.get_as_dict(
                                  "select * from rss_cfg where rss_cfg_id = " +
                                  str(rss_cfg_id) + ";")
      # get the converter configuration for that channel
      self.get_channel_attenuation(self.start, chan)
      self.rss_cfg[chan].update(self.conv_cfg)
    return self.channels
    
  def maps_from_tlogs(self):
    """
    Gets the data for the specified channel and polarization for this map
    """
    try:
      chan_list = self.channels
    except:
      self.channels = self.get_active_channels()
    if self.channels.any():
      pass
    else:
      self.logger.warning("maps_from_tlogs: this map has no active channels")
      return None
    self.map_data = {}
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

  def center_map(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
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
    if source.lower() == "sun":
      src = ephem.Sun()
    else:
      pass
    self.map_data['dec_offset'] = []
    self.map_data['xdec_offset'] = []
    for count in range(len(self.map_data['MPL_datenum'])):
      dt = num2date(self.map_data['MPL_datenum'][count])
      if source.lower() == "sun":
        src.compute(dt)
      else:
        pass
      ra_center = src.ra*12/pi    # hours
      dec_center = src.dec*180/pi # degrees
      decrad = src.dec
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
    self.map_data['grid_x'] = NP.arange( -width/2,  width/2 + step, step/2)
    self.map_data['grid_y'] = NP.arange(-height/2, height/2 + step, step/2)
    self.map_data['grid_z'] = {}
    for channel in self.channels:
      cx = self.map_data['xdec_offset']
      cy = self.map_data['dec_offset']
      cz = self.map_data['VFC_counts'][channel]
      xi = self.map_data['grid_x']
      yi = self.map_data['grid_y']
      try:
        self.map_data['grid_z'][channel] = griddata(cx,cy,cz, xi, yi, interp=u'nn')
      except ValueError, details:
        self.logger.error("regrid: gridding failed: %s", str(details))
        self.logger.debug("regrid: channel %d length of cx is %d", channel, len(cx))
        self.logger.debug("regrid: channel %d length of cy is %d", channel, len(cy))
        self.logger.debug("regrid: channel %d length of cz is %d", channel, len(cz))
        continue
    return self.map_data['grid_x'], self.map_data['grid_y'], \
           self.map_data['grid_z']


class BoresightScan(Observation):
  """
  class for a single scan during a boresight
  
  Attributes::
    axis     - direction of the scan (footnote 1)
    cal_flux - flux of the calibrator source
    cal_src  - source used to set the flux scale
    chan     - channel used for the boresight
    data     - scan data (time, positions and Tsys)
    diode    - state of the noise diode
    epoch    - UNIX time
    freq     - channel frequency in MHz
    IFbw     - IF band width in MHz
    IFmode   - IF phasing
    logger   - logging.Logger object
    name     - identifier string based on xpwr_cfg_id
    pol      - channel polarization
    session  - Session object
    source   - name of scanned source
  """
  def __init__(self, parent, xpwr_cfg_id):
    """
    initialize the class
    
    @param parent : 'self' of the calling method
    @type  parent : Session object
    @param xpwr_cfg_id : row identifier in table 'xpwr_cfg'
    @type  xpwr_cfg_id : int
    """
    self.logger = logging.getLogger(parent.logger.name+".BoresightScan")
    self.session = parent
    self.name = "boresight"+str(xpwr_cfg_id)
    # get the metadata
    columns = "rss_cfg_id, chan, cal_flux, epoch, source_id, cal_src_id, axis"
    meta = self.session.db.get_as_dict(
              "select "+columns
              + " from xpwr_cfg where xpwr_cfg_id="+str(xpwr_cfg_id)+";")
    rss_cfg_id = int(meta['rss_cfg_id'])
    self.chan = int(meta['chan'])
    self.cal_flux = float(meta['cal_flux'])
    self.epoch = float(meta['epoch'])
    source_id = int(meta['source_id'])
    if source_id:  # source ID 0 is no source
      self.source = self.session.db.get_source_names([source_id])['source'][0]
      cal_src_id = int(meta['cal_src_id'])
      self.cal_src = self.session.db.get_source_names([cal_src_id])['source'][0]
      self.axis = meta['axis'][0]
      # receiver metadata
      rss_cols = "sky_freq, feed, pol, nd, if_mode, if_bw"
      rss_cfg = self.session.db.get_as_dict(
              "select "+rss_cols
              + " from rss_cfg where rss_cfg_id="+str(rss_cfg_id)+";")
      self.freq = float(rss_cfg['sky_freq'])
      self.pol = rss_cfg['pol'][0]
      self.diode = rss_cfg['nd'][0]
      self.IFmode = rss_cfg['if_mode'][0]
      self.IFbw = float(rss_cfg['if_bw'])
      # scan data
      scan_cols = "epoch, tsys, az, el, ha, `dec`"
      self.data = self.session.db.get_as_dict(
                "select "+scan_cols
                +" from xpwr where xpwr_cfg_id="+str(xpwr_cfg_id)+";")
      # tlog data
      if self.data.has_key('epoch'):
        start = self.data['epoch'][0]
        end   = self.data['epoch'][-1]
        self.channels = NP.unique(
                    self.session.db.get("select chan from tlog where epoch >= "
                                    +str(start)+" and epoch <= "+str(end)+";"))
        logcols = "epoch, chan, top, integ, az, el, diode, rss_cfg_id"
        self.logdata = {}
        for chan in self.channels:
          response = self.session.db.get("select "+logcols
                                       +" from tlog where epoch >= "
                                     +str(start)+" and epoch <= "+str(end)+";")
          self.logdata[chan] = {}
          self.logdata[chan]['epoch'] = response[:,0][NP.where(response[:,1] == chan)]
          self.logdata[chan]['top']   = response[:,2][NP.where(response[:,1] == chan)]
          self.logdata[chan]['integ'] = response[:,3][NP.where(response[:,1] == chan)]
          self.logdata[chan]['az']    = response[:,4][NP.where(response[:,1] == chan)]
          self.logdata[chan]['el']    = response[:,5][NP.where(response[:,1] == chan)]
          rss_cfg_id = NP.unique(response[:,7][NP.where(response[:,1] == chan)])[0]
          response = self.session.db.get_as_dict("select "
              +rss_cols +" from rss_cfg where rss_cfg_id="+str(rss_cfg_id)+";")
          self.logdata[chan]['feed']    = response['feed'][0]
          self.logdata[chan]['if_bw']   = response['if_bw'][0]
          self.logdata[chan]['if_mode'] = response['if_mode'][0]
          self.logdata[chan]['nd']      = response['nd'][0]
          self.logdata[chan]['pol']     = response['pol'][0]
          self.logdata[chan]['freq']    = response['sky_freq'][0]
          self.logdata.atten = get_channel_attenuation(self, start, chan)
    else:
      delete(self)

class Session(object):
  """
  Class for an observing session on a given year and DOY
  
  Public Attributes::
    boresights    - dict keyed on 'xpwr_cfg_id' with 2D arrays for scan metadata
    bs_channels   - dict keyed on 'xpwr_cfg_id' with lists of active channels
    bs_data       - dict keyed on 'xpwr_cfg_id' with 2D rrays for 'tlog' data
    db            - database
    doy           - day of year for session
    logger        - logging.Logger object
    maps          - maps in this session
    session_dir   - path to results from this session
    xpwr_metadata - 2D array with data for each 'xpwr' configuration
    year          - year for session
  
  Notes on Data Arrays::
    * 'boresights' 2D-arrays have a row for each scan of the boresight and
      columns for::
        0 - 'xscan_id',
        1 - 'xpwr_cfg_id', and
        2 - 'epoch'.
    * 'bs_data' 2D-arrays have a row for each 'tlog' row and columns for:: 
        0 - UNIX time,
        1 - counts,
        2 - integration time,
        3 - azimuth,
        4 - elevation,
        5 - noise diode state, and
        6 - channel [if argument chan=None; see get_boresight_data()]
    * 'xpwr_metadata' is a 2D-array with a row for each configuration and columns::
        0 - 'xpwr_cfg_id'
        1 - UNIX time,
        2 - rss_cfg_id,
        3 - source_id,
        4 - axis, and
        5 - chan
  """
  def __init__(self, parent, year, doy, plotter=False):
    """
    """
    self.logger = logging.getLogger(parent.logger.name+".Session")
    self.db = parent
    self.year = year
    self.doy = doy
    if plotter == False:
      # instantiating map plotters also gets the maps
      self.get_maps()
    self.get_boresights()
    self.get_session_dir()

  def get_session_dir(self):
      obs_dir = projects_dir+"SolarPatrol/Observations/dss28/"
      self.session_dir = obs_dir + "%4d" % self.year +"/"+ "%03d" % self.doy +"/"
      if not os.path.exists(self.session_dir):
        os.makedirs(self.session_dir, mode=0775)

  # ------------------------------ maps ---------------------------------------
  
  def get_maps(self):
    """
    Returns maps from the raster configuration IDs for the specified date
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
    self.logger.info("%4d/%03d found %d maps", self.year, self.doy,
                     len(self.maps.keys()))

  def get_boresights(self):
    """
    Returns boresights from the xpwr configurations
    """
    xpwr_cfg_ids = self.db.get("select xpwr_cfg_id from xpwr_cfg where year = "
                          +str(self.year)+" and doy = "+str(self.doy)+";")[:,0]
    xpwr_cfg_ids.sort()
    self.boresights = {}
    for xpwr_cfg_id in xpwr_cfg_ids:
      self.boresights[xpwr_cfg_id] = BoresightScan(self, xpwr_cfg_id)
    self.logger.info("%4d/%03d found %d boresights", self.year, self.doy,
                     len(self.boresights.keys()))
          
  def list_maps(self, save=False):
    """
    """
    if save:
      fileobj = open(self.session_dir+"maps.txt", "w")
    else:
      fileobj = sys.stdout
    print >> fileobj, "----------------- Session Maps for %4d/%03d -------------------" %\
          (self.year, self.doy)
    print >> fileobj, " ID start-stop ch  freq.  pol.  b.w. IFmode attn.        source"
    print >> fileobj, "--- ---------- -- ------ ----- ----- ------ ----- -------------"
    mapkeys = self.maps.keys()
    mapkeys.sort()
    if mapkeys == []:
      print >> fileobj, "no valid maps with tlog data found"
      return False
    for mapno in self.maps.keys():
      try:
        channels = self.maps[mapno].channels
        for chno in channels:
          print >> fileobj, " %3d %4s-%4s %2d %6.0f %4s %4.2f %4s %4.1d %16s" % (
            mapno,
            strftime("%H%M", gmtime(self.maps[mapno].start)),
            strftime("%H%M", gmtime(self.maps[mapno].end)),
            chno,
            self.maps[mapno].rss_cfg[chno]["sky_freq"],
            self.maps[mapno].rss_cfg[chno]['pol'][0],
            self.maps[mapno].rss_cfg[chno]["if_bw"],
            self.maps[mapno].rss_cfg[chno]["if_mode"][0],
            self.maps[mapno].rss_cfg[chno]["atten"],
            self.maps[mapno].source)
      except AttributeError:
        print "map", mapno, "has no channels"
    return True

  def save_map_data(self, mapkeys=None):
    """
    create a dict with the map data from the designated images
    
    This speeds up retrieval of images
    
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

  # ---------------------------- receiver data --------------------------------
  
  def get_receiver_data(self, time, columns):
    """
    Get the receiver state at a given time

    Columns is a string with column names separated by commas.
    
    This creates a dictionary keyed with channel number and returns a dictionary
    of the receiver configuration, keyed with specified in the columns, that was
    in effect at the given time.
    
    The columns in the 'rss_cfg' table are::
      rss_cfg_id - primary key
      year       -
      doy        -
      utc        -
      epoch      - UNIX time
      chan       -
      sky_freq   -
      feed       -
      pol        -
      nd         -
      if_mode    -
      if_bw      -
      bb_bw      -
      fiber_chan -
    
    Returns a dict of dicts keyed on column name, where the sub-dicts are keyed
    on channel number.

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
    latest_data = self.db.get_as_dict("select rss_cfg_id,year,doy,utc from rss_cfg"
                        +" where epoch <= '"+str(time)
                        +"' order by epoch desc limit 1;")
    # remove any spaces between names and commas
    columns = columns.replace(" ","")
    # split string into a list
    column_keys = columns.split(',')
    cfg_ID = latest_data['rss_cfg_id'][0]
    receiver = {}
    for key in column_keys:
      receiver[key] = {}
      for chan in [2,4,6,8,10,12,14,16]:
        rx_data = self.db.get_as_dict("select "+columns
                       +" from rss_cfg where rss_cfg_id <= "+str(cfg_ID)
                       +" and chan = "+str(chan)
                       +" order by rss_cfg_id desc limit 1;")
                        
        index = column_keys.index(key)
        receiver[key][chan] = rx_data[key][0]
    return receiver

  # --------------------------- method for boresights -------------------------
  
  def get_good_boresights(self):
    """
    Retrieves data from 'tlog' table for boresights with a given channel
    
    Returns a numpy array with columns containing::
      0 - UNIX time
      1 - counts
      2 - integration time
      3 - azimuth
      4 - elevation
      5 - noise diode state
      6 - chan (if chan=None)
    
    """
    keys = self.boresights.keys()
    keys.sort()
    self.good_boresights = {}
    for key in keys:
      self.good_boresights[key] = []
      try:
        channels = list(self.boresights[key].channels)
      except AttributeError:
        self.logger.warning("get_good_boresights: %d has no channels", key)
      else:
        if bool(channels):
          for ch in channels:
            try:
              start = self.boresights[key].logdata[ch]['epoch'][0]
              end   = self.boresights[key].logdata[ch]['epoch'][-1]
            except:
              continue
            self.good_boresights[key].append(ch)
      if self.good_boresights[key] == []:
        self.good_boresights.pop(key)
    return self.good_boresights
  
  def make_bs_dir(self, good_only=False, save=False):
    """
    Notes
    =====
    Each good boresight consists of two scans
    """
    if save:
      fileobj = open(self.session_dir+"xscans.txt", "w")
    else:
      fileobj = sys.stdout
    if good_only:
      bs_keys = self.get_good_boresights().keys()
    else:
      # these are the keys for all boresights, good or bad
      bs_keys = self.boresights.keys()
    bs_keys.sort()
    num_scans = len(bs_keys)
    if num_scans == 0:
      # no data
      print >> fileobj, " Boresight Summary for %4d/%03d" % (self.year, self.doy)
      print >> fileobj, "\nNo valid boresights with tlog data found"
      return False
    print >> fileobj, " Boresight Summary for %4d/%03d" % (self.year, self.doy)
    print >> fileobj, "  ID   date          ch axis  freq.  pol IF bw   source         Top   diode   az    el"
    print >> fileobj, "------ ------------- -- ---- ------ ---- ---- ---------------- ------ ------ ----- ----"
    for bs in bs_keys:
      source =  self.boresights[bs].source
      try:
        bs_channels = self.boresights[bs].channels
      except AttributeError:
        print >> fileobj, "%6d has no channels" % bs
      else:
        bs_channels.sort()
        #self.logger.debug("make_bs_dir: boresight %d channels: %s", bs, bs_channels)
        #self.logger.debug("make_bs_dir: boresight %d channels is %s", bs, bool(bs_channels))
        if bool(bs_channels.any()):
          for ch in bs_channels:
            UNIXtime = self.boresights[bs].epoch
            top =      self.boresights[bs].data['tsys'][0]
            axis =     self.boresights[bs].axis
            az =       self.boresights[bs].data['az'][0]
            el =       self.boresights[bs].data['el'][0]  
            print >> fileobj, "%6d %13s %2s %4s %6.0f %4s %4.0f %16s %6.2f %6s %5.1f %4.1f" % (
                            bs, 
                            strftime("%Y/%j %H%M", gmtime(UNIXtime)),
                            ch, axis,
                            self.boresights[bs].freq,
                            self.boresights[bs].pol,
                            self.boresights[bs].IFbw,
                            source, top, 
                            self.boresights[bs].diode, az, el)
        else:
          print >> fileobj, "%6d has no channels" % bs
    return True


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
    #self.logger.debug("get_receiver_data: requesting %s", column_keys)
    latest_data = self.get_as_dict("select rss_cfg_id,year,doy,utc from rss_cfg"
                        +" where year <= "+str(year)
                        +" and doy <= "+str(doy)
                        +" and utc <= '"+str(time)
                        +"' order by year desc, doy desc, utc desc limit 1;")
    #self.logger.debug("get_receiver_data: query response: %s",latest_data)
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

  def get_source_names(self, source_IDs):
    """
    Get the source information from gavrt_sources.source
    
    Returns a dict with source names for the source IDs provided
    """
    names = {'source': []}
    self.logger.debug("get_source_names: for %s", source_IDs)
    for source_id in source_IDs:
      if source_id:  # no source_id = 0
        response = self.get_as_dict(
                       "select name from gavrt_sources.source where source_id="
                       +str(source_id)+";")
        names['source'].append(response['name'][0])
      else:
        names['source'].append([None])
    return names

