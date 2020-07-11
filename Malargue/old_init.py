# -*- coding: utf-8 -*-
"""
Classes for reducing observations

This was designed to support Malargue observations but can be used for other
observatories (identified by a DSS number) that produce a single text file in
row/column format.

An ``Observation`` object comprises a data structure in this form::

  UNIXtime     NP.float (N,)            seconds since 1970.0
  azimuth      NP.float (N,)            horizon system longitude (deg)
  elevation    NP.float (N,)            horizon system latitude (deg)
  RA           NP.float (N,)            apparent topocentric right ascension
  dec          NP.float (N,)            apparent topocentric declination
  MPL_datenum  NP.float (N,)            number of days since 0001-01-01 UTC +1
                                        ({*e.g.* 0001-01-01, 06:00 is 1.25)
  power        NP.float {ch:(N,),... }  or equivalent, like detector volts,
                                        for each channel
  freq         float    {ch: ?}         frequency of the channel (MHz)
  pol          str      {ch: ?}         polarization of the channel

Notes
=====

Data from a text file::

  * The original data are provided as a 2D text array. Each provided parameter
    must be in a column.
  * If the column names in the text file are not suitable, then there must be
    a dict to map parameter name to column number.
  * The telescope position and other attributes are determined from its DSN
    number designation, *e.g.* ``DSS(28)``.
  * If the mean astrometic geocentric position (RA, dec from a catalog at some
    epoch) is given, then the topocentric apparent celestial position
    (RA, dec for the date of observation) is computed, and then the
    topocentric horizontal position (az, el) is computed.
  * If a topocentric position is given, then the topocentric apparent celestial
    position is computed.

Data from a ``numpy`` (structured) ``ndarray``::

  * If the data are provided as a structured NP array, then the ``names`` item
    of the array must use those defined above.

Additional items::

  * Other items may be defined for the dict but might not be used by the
    ``Observation`` object.
  * ``MPL_datenum`` is computed from ``UNIXtime`` for convenience in plotting
    time series.
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

import Astronomy as A
import Astronomy.Ephem as Aeph
import Data_Reduction as DR
import local_dirs
import Math.least_squares as lsq
import Radio_Astronomy as RA
import support

logger = logging.getLogger(__name__)

Malargue     = Aeph.DSS(84)
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

class Observation(DR.Observation, DR.DataGetterMixin):
  """
  Class for any group of data for a single purpose.
  
  **Attributes**
  
  channels : list of str
    active channels
  conv_cfg : dict
    converter configuration
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
  def __init__(self, parent=None, name=None, # for the superclass
                     dss=None, date=None, project='Malargue',
                     datafile=None, delimiter=" ", # for the data file
                     skip_header=0, names=True, aliases=None, props=None):
    """
    Initialize an Observation object
    
    Args:
      parent (Session):  optional
      name (str):        optional name for this observation
      dss (int):         required station number
      date (str):        YEAR/DOY, required to find session files
      project (str):     required used to find session files
      datafile (str):    name of file in session directory
      delimiter (str):   as defined for ``genfromtxt``
      skip_header (int): rows at start to ignore
      names (bool or list of str): first row has column names
      aliases (list of str): map text column names to data names
      props (dict of dicts): signal properties
    """
    logger.debug("Observation.__init__: initializing...")
    if parent:
        loggername = parent.logger.name+".Observation"
    else:
        loggername = "ObservationLogger"
    self.logger = logging.getLogger(loggername)
    self.session = parent
    DR.Observation.__init__(self, parent=parent, name=name, dss=dss, date=date,
                            project=project)
    # map text file column names to ndarray data names
    self.aliases = {}
    if aliases:
        self.aliases.update(aliases)
    # get the data
    self.datafile = datafile
    if datafile:
      data = self.open_datafile(datafile, delimiter=delimiter, 
                                     names=names, 
                                     skip_header=skip_header)
      numdata = len(data)
      self.logger.debug("__init__: %d samples", numdata)
      names = data.dtype.names
      self.logger.info("__init__: column names = %s", names)
      metadata, signals = self.get_data_channels(data)
      self.make_data_struct(data, metadata, signals) # this makes self.data
      self.make_channels(signals, props=props)
    else:
      self.logger.warning("__init__: must call method 'open_datafile()' next")

  def get_active_channels(self, filename):
    """
    Returns IDs of channels which took data during this observation
    
    This is intended to be Malargue-specific in that the allowed signal names
    are predefined. 
    """
    allowed = ["XL", "XR", "KaL", "KaR", "power"]
    try:
        self.logger.debug("get_active_channels: %s", self.channels)
    except:
        self.channels = []
        fd = open(filename,"rt")
        self.names = fd.readline().strip().split()
        self.logger.debug("get_active_channels: columns=%s", self.names)
        for name in self.names:
            if name in allowed:
                self.channels.append(name)
        fd.close()
    return self.channels

class Map(Observation, DR.Map):
  """
  """
  def __init__(self, parent=None, name=None, dss=84, date=None, 
                     project="SolarPatrol", datafile=None, source="Sun",
                     step=None, props=None):
    """
    put Malargue-specific initialization here
    """
    logger.debug("Map.__init__: initializing...")
    Observation.__init__(self, parent=parent, name=name, dss=dss, date=date,
                         project=project, datafile=datafile, step=step,
                         props=props)
    self.get_offsets(source="Venus")
    
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
        self.maps[key].get_offsets()
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

