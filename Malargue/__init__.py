# -*- coding: utf-8 -*-
"""
Classes for reducing observations

This is barely started by extracting from GAVRT.  UNDER DEVELOPMENT

An ``Observation`` object comprises a data structure in this form:

  UNIXtime     NP.float (N,)           seconds since 1970.0
  azimuth      NP.float (N,)           horizon system longitude (deg)
  elevation    NP.float (N,)           horizon system latitude (deg)
  RA           NP.float (N,)           apparent topocentric right ascension
  dec          NP.float (N,)           apparent topocentric declination
  MPL_datenum  NP.float (N,)           number of days since 0001-01-01 UTC, plus 1
                                       ({*e.g.* 0001-01-01, 06:00 is 1.25)
  power       NP.float  {ch:(N,),... } or equivalent, like detector volts,
                                       for each channel
  freq        float     {ch: ?}        frequency of the channel (MHz)
  pol         str       {ch: ?}        polarization of the channel

Notes
=====

  * If ``azimuth`` and ``elevation`` are given, then apparent ``RA`` and ``dec`` are computed.
  * IF apparent ``RA`` and ``dec`` are given, then ``azimuth`` and ``elevation`` are computed.
  * Other items may be defined for the dict (*e.g* ``RA2000``, ``dec2000``) but might not be
    used by the ``Observation`` object.
  * If the mean astrometic geocentric position is given, then the apparent 
    topocentric position is computed, and then ``azimuth`` and ``elevation``.
  * ``MPL_datenum`` is computed from ``UNIXtime``.
  * The telescope position and other attributes are determined from its DSN
    number desigbation, *e.g.* ``DSS(28)``.  
"""
import pickle as pickle
import datetime
import ephem
import logging
import math
import matplotlib.dates as mpd
import numpy as NP
import os
import stat
import sys
import time

# from matplotlib.mlab import griddata  NEED NEW GRIDDING 

import Astronomy as A
import Astronomy.Ephem as Aeph
import local_dirs
import Math.least_squares as lsq
import Radio_Astronomy as RA
import support

logger = logging.getLogger(__name__)

Malargue        = Aeph.DSS(84)
longitude    = -Malargue.long*180/math.pi
latitude     = Malargue.lat*180/math.pi
f_max        = 32. # GHz for default step size
wl_min       = f_max/300
taper        = 12 # dB
hpbw         = RA.HPBW(taper, wl_min, 34)*180/math.pi # deg
default_step = hpbw/3.

def DSS84_beamtaper(freq):
  """
  ad hoc fit to beamwidth vs frequency plot
  """
  if freq < 7:
    taper=0
  else:
    taper = 50*(log10(freq)-log10(7))
  return taper
  
def DSS84_beamwidth(freq):
  """
  beamwidth in deg. with edge taper
  """
  return RA.HPBW(DSS84_beamtaper(freq), 0.3/float(freq), 34)*180/math.pi

class Observation(object):
  """
  Class for any group of data for a single purpose.
  
  **Attributes**
  
  channels : list of str
    active channels
  conv_cfg : dict
    converter configuration
  data : dict
    result of ``get_data_from_logs()``
  end : float
    observation end UNIX time, provided by the subclasses
  logger
    ``logging.Logger`` instance
  parent : ``Session`` instance
    a collection or group of observations
  rss_cfg : dict
    receiver configuration
  start : float
    observation start UNIX time,  provided by the subclasses
  """
  def __init__(self, parent):
    """
    Initialize an Observation object
    """
    self.logger = logging.getLogger(parent.logger.name+".Observation")
    self.session = parent

  def get_active_channels(self):
    """
    Returns IDs of channels which took data during this observation
    """
    self.channels = None # some or all of ["XL", "XR", "KaL", "KaR"]
    return self.channels
    
  def get_data_from_logs(self):
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
      self.logger.warning("get_data_from_logs: this map has no active channels")
      return None
    self.data = {}
    data = None # replace with numpy.loadtxt with a column for each channel
    for channel in self.channels:
      ch_index = list(self.channels).index(channel)
      if ch_index == 0:
        # first channel only: these are common to all channels
        # actual columns depend on way file is structured
        self.data['UNIXtime']    = data[:,0].astype(float)
        self.data['azimuth']     = data[:,1].astype(float)
        self.data['elevation']   = data[:,2].astype(float)
        self.data['RA']          = []
        self.data['declination'] = []
        self.data['MPL_datenum'] = [] # needed for time series plots
        self.data['power']  = {}
        self.data['freq']        = {}
        self.data['pol']         = {}
        # this only needed if coords are az/el
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
      # only the power differs between channels
      self.data['power'][channel]  = data[:,3].astype(float)
      self.data['freq'][channel] = self.rss_cfg[channel]['sky_freq']
      self.data['pol'][channel] = self.rss_cfg[channel]['pol'][0].upper()
    return self.data

  def get_offsets(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source
    
    If the source is the default, the position of the Sun will be computed for
    the time of each sample. IT SEEMS LIKE A GOOD IDEA TO DO THIS FOR PLANETS
    ALSO.
    
    This adds elements with keys ``xdec_offset`` and ``dec_offset`` to the
    attribute ``data``.  If this attribute does not exist then
    ``get_data_from_logs`` is called first.

    @param source : source at map center
    @type  source : ephem source instance

    @param xdec_ofst : relative X-dec position of sample
    @type  xdec_ofst : float

    @param dec_ofst : relative dec position of sample
    @type  dec_ofst : float
    
    @return: (dxdecs,ddecs) in degrees
    """
    # get the data if not yet read ib
    try:
      list(self.data.keys())
    except AttributeError:
      self.get_data_from_logs()
    if 'MPL_datenum' in self.data:
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
      dt = mpd.num2date(self.data['MPL_datenum'][count])
      if type(src) == Aeph.Quasar:
        pass
      else:
        src.compute(dt)
      ra_center = src.ra*12/math.pi    # hours
      dec_center = src.dec*180/math.pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['RA'][count] - ra_center)*15*math.cos(decrad) )
      self.data['dec_offset'].append(  dec_ofst + 
                        self.data['declination'][count] - dec_center)
    return self.data['xdec_offset'], self.data['dec_offset']
      

class Map(Observation):
  """
  Class for all the data and methods associated with a raster scan map
  
  ** Attributes**
  
  cfg : dict
    raster configuration
  cfg_id : int
    entry in the raster configuration tableshape
    
  map_data : dict
    data from log table; ``tsrc`` is dict keyed on channel
    
  logger : ``logging.Logger`` object
  
  name : str
    map identifier

  rss_cfg
    receiver configuration
  session
    observing session to which this map belongs
  """
  def __init__(self, parent, name=None, source="Sun"):
    """
    initialize a Map object
    """
    self.logger = logging.getLogger(parent.logger.name+".Map")
    self.session = parent
    if name:
      self.name = name
    else:
      self.name = "map %d" % self.cfg_id
    self.logger.debug("__init__: map for configuration %s is %s",
                      self.session, self.name)
    self.source = source
         
  def center_map(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source

    @param body : source at map center
    @type  body : ephem source instance

    @return: (dxdecs,ddecs) in degrees
    """
    try:
      list(self.data.keys())
    except AttributeError:
      self.maps_from_tlogs()
    if 'MPL_datenum' in self.data:
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
      dt = mpd.num2date(self.data['MPL_datenum'][count])
      if source.lower() == "sun":
        src.compute(dt)
      else:
        pass
      ra_center = src.ra*12/math.pi    # hours
      dec_center = src.dec*180/math.pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['RA'][count] - ra_center)*15*math.cos(decrad) )
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
        self.data['grid_z'][channel] = griddata(cx,cy,cz, xi, yi, interp='nn')
      except ValueError as details:
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
  def __init__(self, parent, xpwr_cfg_id, source=None, axis='dec'):
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
    if source:  
      # this is the source at the center of the scan
      self.source = source
      self.logger.debug("__init__: central source is %s", self.source)
      # we assume this source is a recognized calibrator
      self.calibrator = Aeph.calibrator(self.source)
      self.axis = axis
      # code to get data
      
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
    if ('data' in self.__dict__) == False:
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
    beamwidth = DSS84_beamwidth(self.data['freq'][channel]/1000)
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
        lower_baseline = [0, support.nearest_index(x, lower_limit)]
      if upper_limit > x[-5]: # source scan starts inside upper baseline segment
        upper_baseline = [-6,-1] # force 5 baseline points
      else:
        upper_baseline = [support.nearest_index(x, upper_limit), -1]
    else:
      # scans go from high sample to low sample
      if upper_limit > x[5]: 
        upper_baseline = [0, support.nearest_index(x,upper_limit)]
      else:
        upper_baseline = [0,5]
      if lower_limit < x[-5]:
        lower_baseline = [-6,-1]
      else:
        lower_baseline = [support.nearest_index(x,lower_limit), -1]
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
    sigma = lsq.st_dev(beamwidth)
    initial_guess = [height, beam_center, sigma]
    # in this case we only fit out to one beamwidth
    if x[0] < x[-1]:
      xfit =  x[support.nearest_index(x,beam_center-beamwidth):support.nearest_index(x,beam_center+beamwidth)]
      y = zdata[support.nearest_index(x,beam_center-beamwidth):support.nearest_index(x,beam_center+beamwidth)]
    else:
      xfit =  x[support.nearest_index(x,beam_center+beamwidth):support.nearest_index(x,beam_center-beamwidth)]
      y = zdata[support.nearest_index(x,beam_center+beamwidth):support.nearest_index(x,beam_center-beamwidth)]
    self.pars, err = lsq.fit_gaussian(lsq.gaussian_error_function,
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
    * 'boresights' 2D-arrays have a row for each scan of the boresight and columns for::
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
      obs_dir = local_dirs.projects_dir+"SolarPatrol/Observations/Malargue/"
      self.session_dir = obs_dir + "%4d" % self.year +"/"+ "%03d" % self.doy +"/"
      if not os.path.exists(self.session_dir):
        os.makedirs(self.session_dir, mode=0o775)
  
  def summary(self, save=False):
    if not self.list_maps(save=save):
      print("no usable maps found")
    else:
      self.show_images()
    if not self.make_bs_dir(save=save):
      print("no usable boresights found")
    else:
      self.show_boresights()

  # ------------------------------ maps ---------------------------------------
  def get_map_IDs(self):
    """
    """
    map_cfg_ids = []
    self.logger.debug("get_maps: map IDs: %s", map_cfg_ids)
    return map_cfg_ids

  def get_maps(self, map_IDs=[]):
    """
    Returns maps from the raster configuration IDs for the specified date
    """
    if map_IDs == []:
      map_cfg_ids = self.get_map_IDs()
    else:
      map_cfg_ids = NP.array(map_IDs)
    self.logger.debug("get_maps: map IDs: %s", map_cfg_ids)
    if map_cfg_ids.any():
      self.maps = {}
      for map_id in map_cfg_ids[:,0]:
        self.logger.debug("get_maps: getting %d", map_id)
        self.maps[map_id] = Map(self, map_id)
      self.logger.info("%4d/%03d found %d maps", self.year, self.doy,
                     len(list(self.maps.keys())))
    else:
      self.logger.info("No maps found for %4d/%03d", self.year, self.doy)
      
  def get_boresights(self):
    """
    Returns boresights from the xpwr configurations
    """
    self.boresights = {}
    pass
    self.logger.info("%4d/%03d found %d boresights", self.year, self.doy,
                     len(list(self.boresights.keys())))
          
  def list_maps(self, save=False):
    """
    """
    if save:
      fileobj = open(self.session_dir+"maps.txt", "w")
    else:
      fileobj = sys.stdout
    print("----------------- Session Maps for %4d/%03d -------------------" %\
          (self.year, self.doy), file=fileobj)
    print(" ID start-stop ch  freq.  pol.  b.w. IFmode attn.        source", file=fileobj)
    print("--- ---------- -- ------ ----- ----- ------ ----- -------------", file=fileobj)
    mapkeys = list(self.maps.keys())
    mapkeys.sort()
    if mapkeys == []:
      print("no valid maps with tlog data found", file=fileobj)
      return False
    for mapno in list(self.maps.keys()):
      try:
        channels = self.maps[mapno].channels
        for chno in channels:
          print(" %3d %4s-%4s %2d %6.0f %4s %4.2f %4s %4.1d %16s" % (
            mapno,
            time.strftime("%H%M", time.gmtime(self.maps[mapno].start)),
            time.strftime("%H%M", time.gmtime(self.maps[mapno].end)),
            chno,
            self.maps[mapno].rss_cfg[chno]["sky_freq"],
            self.maps[mapno].rss_cfg[chno]['pol'][0],
            self.maps[mapno].rss_cfg[chno]["if_bw"],
            self.maps[mapno].rss_cfg[chno]["if_mode"][0],
            self.maps[mapno].rss_cfg[chno]["atten"],
            self.maps[mapno].source), file=fileobj)
      except AttributeError:
        print("map", mapno, "has no channels")
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
      mapkeys = list(self.maps.keys())
      mapkeys.sort()
    for key in mapkeys:
      try:
        list(self.maps[key].map_data.keys())
        self.logger.debug("save_map_data: mapdata[%d] exists", key)
      except AttributeError:
        self.maps[key].maps_from_tlogs()
        self.logger.debug("save_map_data: loaded mapdata[%d]", key)
      if 'dec_offset' in self.maps[key].map_data:
        self.logger.debug("save_map_data: mapdata[%d] is centered", key)
      else:
        self.maps[key].center_map()
        self.logger.debug("save_map_data: mapdata[%d] has been centered", key)
      if 'grid_x' in self.maps[key].map_data:
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
    **Logic**
    
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
    keys = list(self.boresights.keys())
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
      bs_keys = list(self.get_good_boresights().keys())
    else:
      # these are the keys for all boresights, good or bad
      bs_keys = list(self.boresights.keys())
    bs_keys.sort()
    num_scans = len(bs_keys)
    if num_scans == 0:
      # no data
      print(" Boresight Summary for %4d/%03d" % (self.year, self.doy), file=fileobj)
      print("\nNo valid boresights with tlog data found", file=fileobj)
      return False
    print(" Boresight Summary for %4d/%03d" % (self.year, self.doy), file=fileobj)
    print("  ID   date          ch axis  freq.  pol IF bw   source         Top   diode   az    el", file=fileobj)
    print("------ ------------- -- ---- ------ ---- ---- ---------------- ------ ------ ----- ----", file=fileobj)
    for bs in bs_keys:
      source =  self.boresights[bs].source
      try:
        bs_channels = self.boresights[bs].channels
      except AttributeError:
        print("%6d has no channels" % bs, file=fileobj)
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
            print("%6d %13s %2s %4s %6.0f %4s %4.0f %16s %6.2f %6s %5.1f %4.1f" % (
                            bs, 
                            time.strftime("%Y/%j %H%M", time.gmtime(UNIXtime)),
                            ch, axis,
                            self.boresights[bs].freq,
                            self.boresights[bs].pol,
                            self.boresights[bs].IFbw,
                            source, top, 
                            self.boresights[bs].diode, az, el), file=fileobj)
        else:
          print("%6d has no channels" % bs, file=fileobj)
    return True

