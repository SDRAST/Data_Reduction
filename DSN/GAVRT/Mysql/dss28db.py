# -*- coding: utf-8 -*-
"""
Classes for GAVRT mysql database

Example
=======
DBPlotter (from Data_Reduction.DSN.GAVRT.Mysql.plotter) is used to reduce data
stored in the LCER GAVRT MySQL database.  The following example gets the 
coordinate data for a map made during a given session::
  In [1]: from Data_Reduction.DSN.GAVRT.Mysql.plotter import DBPlotter
  In [2]: pl = DBPlotter()
  In [3]: sp = pl.get_session_plotter(2017,233)
  In [4]: map69data = sp.maps[69].maps_from_tlogs()
  In [5]: xdec,dec = sp.maps[69].center_map()

Databases
=========
The databases and their schemas are described in
http://gsc.lewiscenter.org/data_info/dss28_eac.php.

The server has these databases::
 'dss28_eac'
 'dss28_spec'
 'gavrt_sources'.


Database 'dss28_eac'
--------------------
has these tables::
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

Database 'gavrt_sources'
------------------------
has these tables::
 'catalog',
 'class',
 'source'

Table columns
-------------
'angles' columns::
  angles_id,
  year, doy, utc, epoch, az, el, status
'catalog' columns::
  catalog_id, name
'chan_cfg' columns::
  chan_cfg_id,
  year, doy, utc, epoch, chan, center_freq, tdiode
'class' columns::
  class_id, name, description
'conv_cfg' columns::
  conv_cfg_id,
  year, doy, utc, epoch, converter, mode_a, ifbw_a, bbbw_a, atten_a, mode_b,
  ifbw_b, bbbw_b, atten_b, lock_status
'five_point' columns::
  five_point_id,
  xpwr_cfg_id, year, doy, utc, epoch, source_id, chan, tsrc, az, el, ha, dec, 
  xdec_off, dec_off
'pointing_cfg' columns::
  pointing_cfg_id,
  year, doy, utc, epoch, man, plx, semod, refrctn, delut, model
'raster' columns::
  raster_id,
  raster_cfg_id, year, doy, utc, epoch, xdecoff, decoff, ha, dec, tsrc
'raster_cfg' columns::
  raster_cfg_id,
  rss_cfg_id, year, doy, utc, epoch, source_id, chan, freq, rate, step
'rf_cfg' columns::
  rf_cfg_id,
  year, doy, utc, epoch, feed, diodex, diodey, pol, transfer
'rss_cfg' columns::
  rss_cfg_id,
  year, doy, utc, chan, sky_freq, feed, pol, nd, if_mode, if_bw, bb_bw, fiber_chan
'source' columns::
  source_id, catalog_id, class_id,
  name, RA, Dec, size_dec, size_xdec, reference, aka
'tlog' columns::
  tlog_id,
  rss_cfg_id, year, doy, utc, epoch, chan, top, integ, az, el, diode, level, cryo
'weather' columns::
  weather_id,
  datetime, pressure, temp, humidity, wind_speed, wind_dir
'xpwr' columns::
  xpwr_id,
  xpwr_cfg_id, year, doy, utc, epoch, tsys, az, el, ha, dec, offset
'xpwr_cfg' columns::
  xpwr_cfg_id, 
  rss_cfg_id, source_id, cal_src_id, year, doy, utc, epoch, axis, chan, cal_flux
'xscan' columns::
  xscan_id,
  xpwr_cfg_id, year, doy, utc, epoch, tsrc, stdev, bl_stdev, az, az_offset, el,
  el_offset, ha, dec, offset, bw, corr
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
from numpy import arange, array, unique, where
from time import gmtime, strftime

from local_dirs import projects_dir
import Astronomy as A
from Astronomy.Ephem import calibrator, DSS, Quasar
from Data_Reduction.DSN.GAVRT.Mysql import BaseDB, MysqlException
from Math.least_squares import fit_gaussian, gaussian_error_function, st_dev
from Radio_Astronomy import HPBW
from support import nearest_index

logger = logging.getLogger(__name__)

_host,_user,_pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p", "rb" ))
dss28        = DSS(28)
longitude    = -dss28.long*180/pi
latitude     = dss28.lat*180/pi
f_max        = 16. # GHz
wl_min       = f_max/300
taper        = 12 # dB
hpbw         = HPBW(taper, wl_min, 34)*180/pi # deg
default_step = hpbw/3.

def DSS28_beamtaper(freq):
  """
  ad hoc fit to beamwidth vs frequency plot
  """
  if freq < 7:
    taper=0
  else:
    taper = 50*(log10(freq)-log10(7))
  return taper
  
def DSS28_beamwidth(freq):
  """
  beamwidth in deg. with edge taper
  """
  return HPBW(DSS28_beamtaper(freq), 0.3/float(freq), 34)*180/pi

class Observation(object):
  """
  Class for any group of data for a single purpose.
  
  Attributes::
    channels - active channels
    conv_cfg - converter configuration
    data     - result of get_data_from_tlogs()
    end      - provided by the subclasses
    logger   - logging.Logger instance
    parent   - a collection or group of observations
    rss_cfg  - receiver configuration
    start    - provided by the subclasses
  Methods::
    get_active_channels - returns channels which took tlog data during this map    
    get_data_from_tlogs - re-organizes tlog table data into map form
    get_offsets
    get_channel_attenuation
  """
  def __init__(self, parent):
    """
    """
    self.logger = logging.getLogger(parent.logger.name+".Observation")
    self.session = parent

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

  def get_active_channels(self):
    """
    returns the receiver channels that were active between 'start' and 'end'
    
    This requires attributes 'start' and 'end' to be defined which happens
    during 'BoresightScan' or 'Map' initialization
    
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
      # get the attenuation for that channel
      self.rss_cfg[chan]['atten'] = self.get_channel_attenuation(self.start, chan)
    return self.channels
    
  def get_data_from_tlogs(self):
    """
    Gets the data for the specified channel and polarization for this observation
    """
    try:
      chan_list = self.channels
    except:
      self.channels = self.get_active_channels()
    if self.channels.any():
      pass
    else:
      self.logger.warning("get_data_from_tlogs: this map has no active channels")
      return None
    self.data = {}
    for channel in self.channels:
      ch_index = list(self.channels).index(channel)
      query = "select tlog.epoch, tlog.az, tlog.el, tlog.top, tlog.rss_cfg_id" +               \
            " from tlog, rss_cfg where tlog.rss_cfg_id = rss_cfg.rss_cfg_id" +\
            " and tlog.epoch >= "  + str(self.start) +                        \
            " and tlog.epoch <= "  + str(self.end  ) +                        \
            " and rss_cfg.chan = " + str(channel)      +";"
      self.logger.debug("get_data_from_tlogs: query: %s", query)
      data = self.session.db.get(query)
      # get receiver configuration
      self.rss_cfg[channel] = self.session.db.get_as_dict(
                                  "select * from rss_cfg where rss_cfg_id = " +
                                  str(int(data[:,4][0])) + ";")
      if ch_index == 0:
        # first channel only: these are common to all channels
        self.data['UNIXtime']    = data[:,0].astype(float)
        self.data['azimuth']     = data[:,1].astype(float)
        self.data['elevation']   = data[:,2].astype(float)
        self.data['RA']          = []
        self.data['declination'] = []
        self.data['MPL_datenum'] = []
        self.data['VFC_counts']  = {}
        self.data['freq']        = {}
        self.data['pol']         = {}
        for index in range(len(self.data['UNIXtime'])):
          dt = datetime.datetime.utcfromtimestamp(
                                     self.data['UNIXtime'][index])
          time_tuple = (dt.year,
                        A.day_of_year(dt.year,dt.month,dt.day)
                        + (  dt.hour
                           + dt.minute/60.
                           + dt.second/3600.
                           + dt.microsecond/3600./1e6)/24.)
          ra, dec = A.AzEl_to_RaDec(
                             float(self.data['azimuth'][index]),
                             float(self.data['elevation'][index]),
                             latitude,
                             longitude,
                             time_tuple)
          self.data['RA'].append(ra)
          self.data['declination'].append(dec)
          self.data['MPL_datenum'].append(date2num(dt))
      # only the VFC counts differ between channels
      self.data['VFC_counts'][channel]  = data[:,3].astype(float)
      self.data['freq'][channel] = self.rss_cfg[channel]['sky_freq']
      self.data['pol'][channel] = self.rss_cfg[channel]['pol'][0].upper()
    return self.data

  def get_offsets(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source
    
    If the source is the default, the position of the Sun will be computed for
    the time of each sample. IT SEEMS LIKE A GOOD IDEA TO DO THIS FOR PLANETS
    ALSO.
    
    This adds elements with keys 'xdec_offset' and 'dec_offset' to the
    attribute 'data'.  If this attribute does not exist then
    'get_data_from_tlogs' is called first

    @param source : source at map center
    @type  source : ephem source instance

    @param xdec_ofst : relative X-dec position of sample
    @type  xdec_ofst : float

    @param dec_ofst : relative dec position of sample
    @type  dec_ofst : float
    
    @return: (dxdecs,ddecs) in degrees
    """
    try:
      self.data.keys()
    except AttributeError:
      self.get_data_from_tlogs()
    if self.data.has_key('MPL_datenum'):
      pass
    else:
      # 'tlog' did not have good data
      return None
    if source.lower() == "sun":
      src = ephem.Sun()
    else:
      src = self.calibrator
    self.data['dec_offset'] = []
    self.data['xdec_offset'] = []
    for count in range(len(self.data['MPL_datenum'])):
      dt = num2date(self.data['MPL_datenum'][count])
      if type(src) == Quasar:
        pass
      else:
        src.compute(dt)
      ra_center = src.ra*12/pi    # hours
      dec_center = src.dec*180/pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['RA'][count] - ra_center)*15*cos(decrad) )
      self.data['dec_offset'].append(  dec_ofst + 
                        self.data['declination'][count] - dec_center)
    return self.data['xdec_offset'], self.data['dec_offset']
  
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
    atten = self.conv_cfg[attenuator]
    self.rss_cfg[channel]['atten'] = atten
    return atten
    

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
    get_map_config      - returns a dict with the raster map configuration
    get_raster_data     - gets the data for a raster scan map used for Zplot
    get_raster_keys     - returns rasters associated with a given configuration
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
    
  def maps_from_tlogs(self):
    """
    Gets the data for the specified channel and polarization for this map
    
    For backward compatibility.  map_data is not used anywhere in this module.
    It's used by plot_centered_offsets().
    """
    self.map_data = self.get_data_from_tlogs()
    return self.map_data

  def center_map(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source

    @param body : source at map center
    @type  body : ephem source instance

    @return: (dxdecs,ddecs) in degrees
    """
    try:
      self.data.keys()
    except AttributeError:
      self.maps_from_tlogs()
    if self.data.has_key('MPL_datenum'):
      pass
    else:
      # 'tlog' did not have good data
      return None
    if source.lower() == "sun":
      src = ephem.Sun()
    else:
      pass
    self.data['dec_offset'] = []
    self.data['xdec_offset'] = []
    for count in range(len(self.data['MPL_datenum'])):
      dt = num2date(self.data['MPL_datenum'][count])
      if source.lower() == "sun":
        src.compute(dt)
      else:
        pass
      ra_center = src.ra*12/pi    # hours
      dec_center = src.dec*180/pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['RA'][count] - ra_center)*15*cos(decrad) )
      self.data['dec_offset'].append(  dec_ofst + 
                        self.data['declination'][count] - dec_center)
    return self.data['xdec_offset'], self.data['dec_offset']
  
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
    self.data['grid_x'] = NP.arange( -width/2,  width/2 + step, step/2)
    self.data['grid_y'] = NP.arange(-height/2, height/2 + step, step/2)
    self.data['grid_z'] = {}
    for channel in self.channels:
      cx = self.data['xdec_offset']
      cy = self.data['dec_offset']
      cz = self.data['VFC_counts'][channel]
      xi = self.data['grid_x']
      yi = self.data['grid_y']
      try:
        self.data['grid_z'][channel] = griddata(cx,cy,cz, xi, yi, interp=u'nn')
      except ValueError, details:
        self.logger.error("regrid: gridding failed: %s", str(details))
        self.logger.debug("regrid: channel %d length of cx is %d", channel, len(cx))
        self.logger.debug("regrid: channel %d length of cy is %d", channel, len(cy))
        self.logger.debug("regrid: channel %d length of cz is %d", channel, len(cz))
        continue
    return self.data['grid_x'], self.data['grid_y'], \
           self.data['grid_z']


class BoresightScan(Observation):
  """
  class for a single scan during a boresight
  
  Attributes::
    axis     - direction of the scan
    bs_data  - scan data (time, positions and Tsys) from 'xpwr' table
    cal_flux - flux of the calibrator source
    cal_src  - source used to set the flux scale
    chan     - channel used for the boresight by EAC program 'xant'
    data     - scan data (time, positions and Tsys) from 'tlog' table
    diode    - state of the noise diode
    epoch    - UNIX time start of scan
    freq     - frequency xant used to fit 'xpwr' data, in MHz
    IFbw     - IF band width in MHz
    IFmode   - IF phasing
    logger   - logging.Logger object
    log_data - data from 'tlog' table
    name     - identifier string based on xpwr_cfg_id
    pol      - channel polarization
    session  - parent Session object
    source   - name of scanned source
  """
  def __init__(self, parent, xpwr_cfg_id):
    """
    initialize the class
    
    @param parent : 'self' of the calling method
    @type  parent : Session object
    @param xpwr_cfg_id : row identifier in table 'xpwr_cfg'
    @type  xpwr_cfg_id : int
    
    From examination of the data it was concluded that the correct values of
    'axis' and 'source-Id' are those for row 'xpwr_cfg_id" + 1.
    """
    Observation.__init__(self, parent)
    self.logger = logging.getLogger(parent.logger.name+".BoresightScan")
    self.name = "boresight"+str(xpwr_cfg_id)
    # get the metadata
    #   see docstring for this ad hoc fix
    columns = "rss_cfg_id, chan, cal_flux, epoch, source_id, cal_src_id, axis"
    #good_columns = "rss_cfg_id, chan, cal_flux, epoch, source_id, cal_src_id, axis"
    meta = self.session.db.get_as_dict(
              "select "+columns
              + " from xpwr_cfg where xpwr_cfg_id="+str(xpwr_cfg_id)+";")
    rss_cfg_id = int(meta['rss_cfg_id'])
    self.cal_flux = float(meta['cal_flux'])
    self.epoch = float(meta['epoch'])
    source_id = int(meta['source_id'])
    if source_id:  # source ID 0 is no source
      # this is the source at the center of the scan
      self.source = self.session.db.get_source_names([source_id])['source'][0]
      self.logger.debug("__init__: central source is %s", self.source)
      self.calibrator = calibrator(self.source)
      # this is the source whose flux is used for calibration
      cal_src_id = int(meta['cal_src_id'])
      self.cal_src = self.session.db.get_source_names([cal_src_id])['source'][0]
      self.axis = meta['axis'][0]
      # receiver metadata
      #self.rss_cfg = self.session.db.get_as_dict(
      #        "select "+rss_cols
      #        + " from rss_cfg where rss_cfg_id="+str(rss_cfg_id)+";")
      # scan data
      scan_cols = "epoch, tsys, az, el, ha, `dec`"
      self.bs_data = self.session.db.get_as_dict(
                "select "+scan_cols
                +" from xpwr where xpwr_cfg_id="+str(xpwr_cfg_id)+";")
      # tlog data
      if self.bs_data.has_key('epoch'):
        self.start = self.bs_data['epoch'][0]
        self.end   = self.bs_data['epoch'][-1]
        self.channels = NP.unique(
                    self.session.db.get("select chan from tlog where epoch >= "
                          +str(self.start)+" and epoch <= "+str(self.end)+";"))
        logcols = "epoch, chan, top, integ, az, el, diode, rss_cfg_id"
        rss_cols = "sky_freq, feed, pol, nd, if_mode, if_bw"
        self.logdata = {}
        self.channels = self.get_active_channels()
        for chan in self.channels:
          response = self.session.db.get("select "+logcols
                                       +" from tlog where epoch >= "
                           +str(self.start)+" and epoch <= "+str(self.end)+";")
          channel_selector = NP.where(response[:,1] == chan)
          self.logdata[chan] = {}
          self.logdata[chan]['epoch'] = response[:,0][channel_selector]
          self.logdata[chan]['top']   = response[:,2][channel_selector]
          self.logdata[chan]['integ'] = response[:,3][channel_selector]
          self.logdata[chan]['az']    = response[:,4][channel_selector]
          self.logdata[chan]['el']    = response[:,5][channel_selector]
          rss_cfg_id = NP.unique(response[:,7][channel_selector])[0]
          response = self.session.db.get_as_dict("select "
              +rss_cols +" from rss_cfg where rss_cfg_id="+str(rss_cfg_id)+";")
          self.logdata[chan]['feed']    = response['feed'][0]
          self.logdata[chan]['if_bw']   = response['if_bw'][0]
          self.logdata[chan]['if_mode'] = response['if_mode'][0]
          self.logdata[chan]['nd']      = response['nd'][0]
          self.logdata[chan]['pol']     = response['pol'][0]
          self.logdata[chan]['freq']    = response['sky_freq'][0]
          self.logdata[chan]['atten']   = self.get_channel_attenuation(
                                                              self.start, chan)
    else:
      delete(self)
    
  def fit_gaussian(self, channel, beam_limit=2.5):
    """
    Fit the scan to a Gaussian and a baseline
    
    Extract the appropriate data::
        For raster scans, 'xdec' means that 'xdec' stays fixed while the
        antenna moves up and down; 'dec' means that 'dec' stays fixed while the
        left and right.
    The Gaussian is assumed to fit the inner five beamwidths of the data,
    though that limit can be adjusted.  The baseline is the rest of the data,
    although the lower baseline includes at least data[:5] and the upper
    baseline includes data[-5:]
    
    @param channel : channel whose data will be fit (required)
    @type  channel : int
    
    @param beam_limit : distance from the center included in Gaussian fit
    @type  beam_limit : float
    """
    self.logger.debug("fit_gaussian: direction is %s", self.axis)
    # get offsets if necessary
    if self.__dict__.has_key('data') == False:
      self.get_offsets()
    # remember that GAVRT nomenclature seems backwards
    if self.axis.lower() == 'xdec':
      x = NP.array(self.data['dec_offset'])  # NP.array(self.ddecs)
    else:
      x = NP.array(self.data['xdec_offset']) # NP.array(self.dxdecs)
    self.logger.debug("fit_gaussian: selected x: %s", x)
    tsys = self.data['VFC_counts'][channel]
    # define the domain of the Gaussian fit:
    beam_index = tsys.argmax()                  # NP.array(self.data).argmax()
    self.logger.debug("fit_gaussian: peak at index %d", beam_index)
    beam_center = x[beam_index]
    self.logger.debug("fit_gaussian: peak at x = %f", beam_center)
    beamwidth = DSS28_beamwidth(self.data['freq'][channel]/1000)
    self.logger.debug("fit_gaussian: beamwidth = %f deg", beamwidth)
    lower_limit = beam_center - beam_limit*beamwidth # source scan starts here
    upper_limit = beam_center + beam_limit*beamwidth # source scan ends here
    self.logger.debug("fit_gaussian: scan lower limit: %f", lower_limit)
    self.logger.debug("fit_gaussian: scan upper limit: %f", upper_limit)
    # Define baseline ranges for the lower end and the upper end of the spectrum
    #  * 'lower_baseline' and 'upper_baseline' are 2-item lists
    #  * assume that there are at least 5 data points for each baseline section
    if x[0] < x[-1]: # increasing X-coordinate
      # scans go from low sample to high sample
      if lower_limit < x[5]: # source scan starts inside lower baseline segment
        lower_baseline = [0,5] # force 5 baseline points
      else:
        lower_baseline = [0, nearest_index(x, lower_limit)]
      if upper_limit > x[-5]: # source scan starts inside upper baseline segment
        upper_baseline = [-6,-1] # force 5 baseline points
      else:
        upper_baseline = [nearest_index(x, upper_limit), -1]
    else:
      # scans go from high sample to low sample
      if upper_limit > x[5]: 
        upper_baseline = [0, nearest_index(x,upper_limit)]
      else:
        upper_baseline = [0,5]
      if lower_limit < x[-5]:
        lower_baseline = [-6,-1]
      else:
        lower_baseline = [nearest_index(x,lower_limit), -1]
    self.logger.debug("fit_gaussian: lower baseline: %s", lower_baseline)
    self.logger.debug("fit_gaussian: upper baseline: %s", upper_baseline)
    # define the baseline data
    xdata = NP.append(x[lower_baseline[0]:lower_baseline[1]],
                      x[upper_baseline[0]:upper_baseline[1]]).astype(float)
    ydata = NP.append(tsys[lower_baseline[0]:lower_baseline[1]],
                      tsys[upper_baseline[0]:upper_baseline[1]]).astype(float)
    #   Fit baseline
    self.baseline_pars = NP.polyfit(xdata,ydata,1)
    self.logger.debug("fit_gaussian: baseline parameters: %s", self.baseline_pars)
    #   Fit the beam
    zdata = NP.array(tsys).astype(float)
    self.logger.debug("fit_gaussian: zdata: %s", zdata)
    height = zdata[beam_index] - NP.polyval(self.baseline_pars, x[beam_index])
    self.logger.debug("fit_gaussian: height: %s", height)
    sigma = st_dev(beamwidth)
    initial_guess = [height, beam_center, sigma]
    # in this case we only fit out to one beamwidth
    if x[0] < x[-1]:
      xfit =  x[nearest_index(x,beam_center-beamwidth):nearest_index(x,beam_center+beamwidth)]
      y = zdata[nearest_index(x,beam_center-beamwidth):nearest_index(x,beam_center+beamwidth)]
    else:
      xfit =  x[nearest_index(x,beam_center+beamwidth):nearest_index(x,beam_center-beamwidth)]
      y = zdata[nearest_index(x,beam_center+beamwidth):nearest_index(x,beam_center-beamwidth)]
    self.pars, err = fit_gaussian(gaussian_error_function,
                            initial_guess,
                            xfit,
                            y-NP.polyval(self.baseline_pars,xfit))
    return self.baseline_pars, self.pars, err

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
    if parent:
      self.logger = logging.getLogger(parent.logger.name+".Session")
    else:
      self.logger = logging.getLogger(logger.name+".Session")
    if parent:
      self.db = parent
    else:
      self.db = DSS28db() # default is GAVRT
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
  
  def summary(self, save=False):
    if not self.list_maps(save=save):
      print "no usable maps found"
    else:
      self.show_images()
    if not self.make_bs_dir(save=save):
      print "no usable boresights found"
    else:
      self.show_boresights()

  # ------------------------------ maps ---------------------------------------
  def get_map_IDs(self):
    """
    """
    map_cfg_ids = self.db.get(
                         "select raster_cfg_id from raster_cfg where year = " +
                         str(self.year) + " and doy = " + str(self.doy) + 
                         ";")
    self.logger.debug("get_maps: map IDs: %s", map_cfg_ids)
    return map_cfg_ids

  def get_maps(self, map_IDs=[]):
    """
    Returns maps from the raster configuration IDs for the specified date
    """
    if map_IDs == []:
      map_cfg_ids = self.get_map_IDs()
    else:
      map_cfg_ids = array(map_IDs)
    self.logger.debug("get_maps: map IDs: %s", map_cfg_ids)
    if map_cfg_ids.any():
      self.maps = {}
      for map_id in map_cfg_ids[:,0]:
        self.logger.debug("get_maps: getting %d", map_id)
        self.maps[map_id] = Map(self, map_id)
      self.logger.info("%4d/%03d found %d maps", self.year, self.doy,
                     len(self.maps.keys()))
    else:
      self.logger.info("No maps found for %4d/%03d", self.year, self.doy)
      
  def get_boresights(self):
    """
    Returns boresights from the xpwr configurations
    """
    try:
      xpwr_cfg_ids = self.db.get("select xpwr_cfg_id from xpwr_cfg where year = "
                          +str(self.year)+" and doy = "+str(self.doy)+";")[:,0]
    except IndexError:
      # 'too many indices for array' means no data were returned
      xpwr_cfg_ids =  []
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
            top =      self.boresights[bs].bs_data['tsys'][0]
            axis =     self.boresights[bs].axis
            az =       self.boresights[bs].bs_data['az'][0]
            el =       self.boresights[bs].bs_data['el'][0]  
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

