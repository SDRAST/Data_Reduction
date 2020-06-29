# -*- coding: utf-8 -*-
"""
Modules to support data reduction in Python.

The base module is intended for reading data from a text file as a ``numpy``
structured array.  The function ``examine_text_data_file()`` reveals the
structure of the file.

Examples
========
Using the original header line::

  obs = Observation(dss=28, datafile=datadir+'t12127.10', 
                    delimiter=[17,16,3,11,7,9,8,2,6], skip_header=1,
                    names="UTC Epoch Chan Tsys Int Az El Diode Level".split())
                      
In this case the column "UTC" comprises the first three data columns.

Using a new set of headers::

  obs2 = Observation(dss=28, datafile=datadir+'t12127.10', 
                     skip_header=1,
          names="Year DOY UTC Epoch Chan Tsys Integr Az El Diode Level".split())

This adds names for the first two columns and replaces "Int" with "Integr".

Note
====

Still IN PROGRESS!
"""
# standard Python modules
import datetime
import glob
import logging
import math
import matplotlib.dates as MPLd
import numpy as NP
import os
import re
import readline
import scipy.fftpack
from scipy.interpolate import griddata

import Astronomy as A
import Astronomy.Ephem as AE
import DatesTimes as DT
from local_dirs import fits_dir, hdf5_dir, projects_dir, wvsr_dir
import support

# enable raw_input Tab completion
readline.parse_and_bind("tab: complete")

logger = logging.getLogger(__name__) # module logger

def examine_text_data_file(filename):
  """
  Examine a file to guide ``genfromtxt()``
  
  Things to look for::
  
    * Is there a header line with column names? If not, use argument ``names``.
    * Is the number of names equal to the number of columns? If not::
      - use argument ``names`` and ``skip_header=1``, or
      - use argument ``delimiter`` with a list of column widths 
        and ``skip_header=1``.
  """
  print(examine_text_data_file.__doc__)
  fd = open(filename, "r")
  lines = fd.readlines()
  fd.close()
  topline = lines[0].strip().split()
  print("          1         2         3         4         5         6         7")
  print("01234567890123456789012345678901234567890123456789012345678901234567890123456789")
  print(lines[0].strip())
  print(lines[1].strip())
  print(" ...")
  print(lines[-1].strip())
  data = NP.genfromtxt(filename, dtype=None, names=None, skip_header=1, encoding=None)
  print("%d datatypes:" % len(data.dtype.fields))
  for item in data.dtype.fields:
    print(item, data.dtype.fields[item])
  
class Observation(object):
  """
  superclass for a data structure and methods
  
  Attributes
  ==========
    aliases      - (dict) data keys to replace column names
    channel      - (dict) signal paths
    data         - (dict) table contents
    end          - (float) UNIX time at the end
    latitude     - (float) from obs
    logger       - (logging.Logger)
    longitude    - (float) from obs
    name         - (str) user assigned, defaults to YEAR/DOY
    numdata      - (int) number of data samples
    obs          - (AE.DSS) observatory
    session      - (Session) set of observations
    session_path - (str) directory for session files
  
  **Reserved Column Names**
  
  These column names are recognized.
  
  These quantities must be present in some form::
  
    unixtime   (float) UNIX time in sec
    chan_name  (str)   channel name
    integr     (float) integration (exposure) in sec
    azel       (float,float) azimuth and elevation in decimal deg
    
  Alternative for ``unixtime``::
  
    year       (int) year of observation 
    doy        (int) day of year
    utc        (str) HH:MM:SS
  
  Alternative for ``unixtime``::
  
    timestr    (str) something like 2020/06/14/14:22:21.00
  
  Alternative for ``chan_name``::
  
    chan       (int) index in receiver channel names

  Alternative for ``azel``::
  
    radec      (float,float) precessed right ascension in decimal hours and
               precessed declination in decimal deg
    radec1950  (float,float) mean right ascension in decimal hours and
               mean declination in decimal deg at epoch
    radec2000  (float,float) mean right ascension in decimal hours and
               mean declination at epoch in decimal deg
    az         (float) azimuth in decimal deg
    el         (float) elevation in decimal deg
    ra         (float) precessed right ascension in decimal hours
    dec        (float) precessed declination in decimal deg
    ra1950     (float) mean right ascension in decimal hours at epoch
    dec1950    (float) mean declination in decimal deg at epoch
    ra2000     (float) mean right ascension in decimal hours at epoch
    dec2000    (float) mean declination in decimal deg at epoch
    
  Optional::
  
    tsys       (float) power level if only a single channel
    top        (float) alternative for ``tsys``
    diode      (float) 0 or power in K (integers OK)
    level      (float) (unidentified -- in ``tlog`` table)
    cryotemp   (float) cryostat temp in K
    windspeed  (float) km/hr
    winddir    (float) deg
    ambtemp    (float) deg C
    pressure   (float) mbar
  
  Columns to be computed::
  
    mpldatenum (float) matplotlib ``datenum``
  
  Notes
  =====
  * The ``data`` structure is a dict.
  * The value of a ``data`` item is either a numpy array or a object 
    like ``float``, ``int``, or ``str``.
  * The keys have reserved words defined above and will be lowercase.
  * Items with other keys may be added, typically by a child class.
  * Coordinates shall be in pairs, `e.g. ``azel``, ``radec``. (This way you 
    never get one without the other.)
  """
  reserved = ['unixtime','chan_name','integr','az','el','year','doy','utc',
              'timestr','chan','tsys','top','diode','level','cryotemp',
              'windspeed','winddir','ambtemp','pressure',
              'ra','dec','ra1950','dec1950','ra2000','dec2000']
              
  def __init__(self, parent=None, name=None, dss=None, date=None, project=None):
    """
    Initialize the base Observation class.
    
    This is not meant to be initialized by itself.  The subclass determines how
    data are read in.
    
    Args:
      parent (Session):   session to which thos ibservation belongs
      name (str):         an identifier; default is station ID + "obs"
      dss (int):          station number
      date (str):         "YEAR/DOY"
      project (str):      directory under /usr/local/projects
    """
    logger.debug("Observation.__init__: initializing...")
    self.session = parent
    if name:
      self.name = name
    else:
      self.name = "DSS"+str(dss)+"obs"
    self.logger = logging.getLogger(logger.name+".Observation")
    if dss:
      self.obs = AE.DSS(dss)
      self.longitude = self.obs.long
      self.latitude = self.obs.lat
    else:
      self.logger.error("__init__: requires observatory location")
      raise Exception("Where were the data taken?")
    if date:
      y,d = date.split('/')
      self.year = int(y); self.DOY = int(d)
      projdatapath, self.sessionpath, rawdatapath = \
                              get_obs_dirs(project, dss, self.year, self.DOY,
                                           datafmt=None)
      self.logger.debug("__init__: session path: %s", self.sessionpath)
    else:
      self.logger.error("__init__: requires a date")
      raise Exception("When were the date taken?")
    if project:
      self.project = project
    else:
      self.logger.error("__init__: requires a project")
      raise Exception("Where are the session's working files?")
    # accomodate subclass arguments
    self.aliases = {}
      
  def open_datafile(self, filename, delimiter=" ", names=True, skip_header=0):
    """
    Opens and reads a data file
    
    This is used by ``Malargue`` (one data files) and ``GAVRT`` (one data file
    for each signal).
    
    Args:
      filename (str):    text data file name
      delimiter (str):   separator between columns (default: whitespace)
      names (bool):      file row has column names (default: True)
      skip_header (int): number of rows to skip at beginning of file
    
    Returns:
      ndarray: 
    """
    data = NP.genfromtxt(filename, 
                         delimiter=delimiter,
                         dtype=None,
                         names=names,
                         case_sensitive='lower', 
                         skip_header=skip_header,
                         encoding=None)
    return data
  
  def get_data_channels(self, data, ignore=None):
    """
    Gets or sets the names of the signal columns
    
    Column names are separated into metadata and signals. Names in
    ``ignore`` re ignored. Names in ``aliases`` are replaced.
    Args:
      data (ndarray): data read from text file
      ignore (list of str): columns to ignore; default None
    
    Returns:
      (list of str, list of str): metadata, signals
    """
    names = data.dtype.names
    metadata = []
    signals = []
    for name in names:
      if ignore:
        if name in ignore:
          pass
      if name.casefold() in map(str.casefold, self.aliases):
        key = self.aliases[name].lower() # we use only lower case names
      else:
        key = name.lower()
      self.logger.debug("get_data_channels: doing %s for %s", key, name)
      if key in map(str.casefold, Observation.reserved):
        if key.casefold() in ['top', 'tsys']:
          signals.append(key)
        else:
          metadata.append(key)
      else:
        signals.append(key)
    self.logger.debug("get_data_channels: signals: %s", signals)
    self.logger.debug("get_data_channels: metadata: %s", metadata)
    return metadata, signals  
      
  def make_channels(self, signals, props=None):
    """
    Assign properties to the channels.  
    
    The prop keys are "freq", "pol", and "IFtype".
    
    Args:
      props (dict of dicts): signal channel properties.
    """
    self.channel = {}
    for ch in signals:
      chindex = signals.index(ch)
      if props:
        self.channel[ch] = self.Channel(self, ch, 
                                        freq=props['freq'][ch],
                                        bw=props['bw'][ch],
                                        pol=props['pol'][ch],
                                        IFtype=props['IFtype'][ch])
      else:
        self.channel[ch] = self.Channel(self, ch)
      if props:
        for prop in props:
          self.channel[ch][prop] = props[prop]

  def make_data_struct(self, data, metadata, signals):
    """
    Takes a text table with headers and converts it into a numpy ``ndarray``.
    
    That means that a column can be extracted using `data[label]`.

    Ideally, there must be a header word for each column and the separator must
    not be  ambiguous, for example, a space if the column headers also include
    spaces.  However, the ``numpy`` function ``genfromtxt()`` is very flexible.
    Here is an example of reading DSS-28 text files from solar observations in
    2012, 9 columns but 7 header words, when database access was not yet 
    practical::
    
      In [1]: from Data_Reduction import Observation
      In [2]: from Astronomy.DSN_coordinates import DSS
      In [3]: import logging
      In [4]: logger = logging.getLogger()
      In [5]: logger.setLevel(logging.INFO)
      In [6]: obs = Observation(dss=28)

    The ``names`` string was just copied by hand from the first line in the
    file. It could also be done with ``file.readline()``.
    
    However, the subclass ``GAVRT.Observation`` will read all the t-files in a
    directory and instead of ``Tsys`` has ``counts`` which is itself a dict with
    an item for each receiver channel.
    
    By default, columns with unknown names are assumed to be power from the
    channel whose name is the column name. Columns that have unknown names
    (i.e. not in 'reserved') must have their names in 'ignore'.

    @param filename : name of the data file
    @type  filename : str

    @param delimiter : column separator
    @type  delimiter : str

    @param aliases : alternative names for keywords
    @type  aliases : dict
    
    @return: column headers (str), data array(str)
    """
    # get the known columns:
    self.data = {}
    self.numdata = len(data)
    self.logger.debug("make_data_struct: using aliases: %s", self.aliases)
    # get columns that are better computed once instead of as needed
    for signal in signals:
      self.logger.debug("make_data_struct: for signal: %s", signal)
      if signal in self.aliases.items():
        self.data[signal] = data[next(key for key, value in self.aliases.items()
                                                       if value == signal)][idx]
      else:
        self.data[signal] = data[signal]
    if 'unixtime' in metadata:
      if 'unixtime' in data.dtype.names:
        self.data['unixtime'] = data['unixtime']
      else:
        self.data['unixtime'] = data[next(key 
                                          for key, value in self.aliases.items()
                                          if value == 'unixtime')]
      self.data['datetime'] = []
      self.data['date_num'] = []
      for idx in list(range(self.numdata)):
        if 'unixtime' in data.dtype.names:
          tm = data['unixtime'][idx]
        else:
          tm = data[next(key for key, value in self.aliases.items()
                                                   if value == 'unixtime')][idx]          
        dt = datetime.datetime.utcfromtimestamp(tm)
        self.data['datetime'].append(dt)
        self.data['date_num'].append(MPLd.date2num(dt))
    if self.check_for(data, 'azel'):
      # azel exists; compute radec if needed; then radec2000 if needed
      if self.check_for(data, 'radec'):
        pass
      else:
        self.radec_from_azel()
        if self.check_for(data, 'radec2000'):
          # ra2000 and dec2000 already exist
          pass
        else:
          self.radec2000_from_radec()
    elif self.check_for(data, 'radec2000'):
      # coordinates exist; compute back to azimuth and elevation
      if self.check_for(data, 'radec'):
        pass
      else:
        # compute observed RA and dec
        self.radec_from_radec2000()
        if self.check_for(data, 'azel'):
          pass
        else:
          self.azel_from_radec()
    # in here check for 'radec'
    else:
      self.logger.error("no coordinates found in data")
      raise Exception("check INFO logging for columns found")
    self.start = self.data['unixtime'].min()
    self.end   = self.data['unixtime'].max()
      
  def splitkey(self, longlat):
    """
    Checks for presence of coordinates in pairs or singles
    
    @param longlat : "azel", or "radec", or "radecEPOC"
    @type  longlat : str
    """
    longitude = longlat[:2] # 'az' or 'ra'
    if len(longlat) > 5: # has epoch
        epoch = longlat[-4:]
        longitude += epoch
        latitude = longlat[2:-4]+epoch
    else: # date of observation
        latitude = longlat[2:]
        epoch = None
    return longitude, latitude, epoch
  
  def check_for(self, data, longlat):
    """
    Checks for separate coordinates and splits if coord pairs
    
    Args:
      data (dict):    attribute ``data``
      longlat (str): "azel", or "radec", or "radecEPOC"
    """
    longitude, latitude, epoch = self.splitkey(longlat)
    if longitude in data.dtype.names and \
       latitude  in data.dtype.names:
      self.logger.debug("check_for: data has %s and %s", longitude, latitude)
      self.data[longitude] = data[longitude]
      self.data[latitude]  = data[latitude]
      return True
    elif longlat in data.dtype.names:
      self.logger.debug("check_for: data has %s", longlat)
      self.data[longitude],self.data[latitude] = map(None, *data[longlat])
      self.logger.debug("check_for: added %s and %s to data",
                        longitude, latitude)
      return True
    else:
      # coords need to be computed from other coords
      
      return False
  
  def radec_from_azel(self):
    """
    compute RA and dec from az and el
    """
    RA = []; decs = []; RAdecs = []
    for idx in list(range(self.numdata)):
      # setup
      dt = self.data['datetime'][idx]
      time_tuple = (dt.year,
                    DT.day_of_year(dt.year,dt.month,dt.day)
                      + (  dt.hour
                         + dt.minute/60.
                         + dt.second/3600.
                         + dt.microsecond/3600./1e6)/24.)
      azimuth = self.data['az'][idx]
      elevation = self.data['el'][idx]
      # compute
      ra,dec = A.AzEl_to_RaDec(azimuth, elevation,
                               self.latitude, self.longitude,
                               time_tuple)
      RA.append(ra)
      decs.append(dec)
      RAdecs.append((RA,decs))
    self.data['ra'] = RA
    self.data['dec'] = dec
    self.data['radec'] = RAdecs
  
  def radec2000_from_radec(self):
    """
    compute RA2000 and dec2000 from observed RA and dec
    """
    RA2000 = []; decs2000 = []; RAdec2000 = []
    for idx in list(range(self.numdata)):
      # setup
      tm = self.data['unixtime'][idx]
      mjd = DT.UnixTime_to_MJD(tm)
      MJD = int(mjd)
      UT = 24*(mjd-MJD)
      ra = self.data['ra']
      dec = self.data['dec']
      # compute
      ra2000,dec2000 = A.apparent_to_J2000(MJD,UT, 
                                           ra, dec,
                                           self.longitude, self.latitude)
      RA2000.append(ra2000)
      decs2000.append(dec2000)
      RAdec2000.append((ra2000,dec2000))
    self.data['ra2000'] = RA2000
    self.data['dec2000'] = dec2000
    self.data['radec2000'] = RAdec2000

  def radec_from_radec2000(self):
    """
    compute apparent RA and dec. from J2000 RA and dec
    """
    RA = []; decs = []; RAdecs = []
    for idx in list(range(self.numdata)):
      # setup
      tm = self.data['unixtime'][idx]
      mjd = DT.UnixTime_to_MJD(tm)
      MJD = int(mjd)
      UT = 24*(mjd-MJD)
      ra2000 = self.data['ra2000'][idx]
      dec2000 = self.data['dec2000'][idx]
      # compute
      ra, dec = A.J2000_to_apparent(MJD, UT, 
                                    ra2000*math.pi/12, dec2000*math.pi/180)
      RA.append(ra)
      decs.append(dec)
      RAdecs.append((ra,dec))
    self.data['ra'] = RA
    self.data['dec'] = decs
    self.data['radec'] = RAdecs

  def azel_from_radec(self):
    """
    compute azimuth and elevation from apparent right ascension and declination
    """
    azs = []; els = []; azels = []
    for idx in list(range(self.numdata)):
      # setup
      ra = self.data['ra'][idx]
      dec = self.data['dec'][idx]
      timetuple = self.data['datetime'][idx].timetuple()
      year = timetuple.tm_year
      doy = timetuple.tm_yday + (timetuple.tm_hour
                              +(timetuple.tm_min+timetuple.tm_sec/60)/60)/24
      # compute
      az, el = A.RaDec_to_AzEl(ra, dec, 
                                     self.latitude, self.longitude, (year,doy))
      azs.append(az)
      els.append(el)
      azels.append((az,el))
    self.data['az'] = azs
    self.data['el'] = els
    self.data['azel'] = azels
    
  def get_offsets(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source
    
    If the source is the default, the position of the Sun will be computed for
    the time of each sample. IT SEEMS LIKE A GOOD IDEA TO DO THIS FOR PLANETS
    ALSO.
    
    This adds elements with keys ``xdec_offset`` and ``dec_offset`` to the
    attribute ``data``.

    @param source : source at map center
    @type  source : ephem source instance

    @param xdec_ofst : relative X-dec position of sample
    @type  xdec_ofst : float

    @param dec_ofst : relative dec position of sample
    @type  dec_ofst : float
    
    @return: (dxdecs,ddecs) in degrees
    """
    if 'date_num' in self.data:
      pass
    else:
      self.logger.debug("get_offsets: this does not look like a good dataset")
      raise Exception("check the datafile")
    if source.lower() == "sun":
      src = AE.ephem.Sun()
    else:
      src = AE.calibrator(source)
    self.data['dec_offset'] = []
    self.data['xdec_offset'] = []
    for count in range(len(self.data['date_num'])):
      dt = MPLd.num2date(self.data['date_num'][count])
      if type(src) == AE.Quasar:
        pass
      else:
        src.compute(dt)
      ra_center = src.ra*12/math.pi    # hours
      dec_center = src.dec*180/math.pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['ra'][count] - ra_center)*15*math.cos(decrad) )
      self.data['dec_offset'].append(  dec_ofst + 
                        self.data['dec'][count] - dec_center)
    # change list to NP.array
    self.data['xdec_offset'] = NP.array(self.data['xdec_offset'])
    self.data['dec_offset'] = NP.array(self.data['dec_offset'])
    return self.data['xdec_offset'], self.data['dec_offset']
    
  def unpack_to_complex(self, rawdata):
    """
    Converts a sequence of alternating real/imag samples to complex

    @param rawdata : alternating real and imaginary bytes
    @type  rawdata : numpy array of signed int8

    @return: numpy array of complex
    """
    datalen = len(rawdata)
    real = rawdata[0:datalen:2]
    imag = rawdata[1:datalen:2]
    data = real + 1j*imag
    return data

  def sideband_separate(self, data):
    """
    Converts a complex spectrum array and returns two reals with USB and LSB

    This applies a Hilbert transform to the complex data.
    """
    usb = (data.real + scipy.fftpack.hilbert(data).imag)
    lsb = (scipy.fftpack.hilbert(data).real + data.imag)
    return lsb,usb
  
  class Channel(support.PropertiedClass):
    """
    Class for a signal path
    """
    def __init__(self, parent, name, freq=None, bw=None, pol=None, IFtype=None):
      """
      Notes
      =====
      The properties can be accessed as if the class were a dict.
      
      Arguments
      =========
      freq:float or int: center frequency in MHz
      bw:float or int:   bandwidth in MHz
      pol:str:           polarization code
      """
      support.PropertiedClass.__init__(self)
      self.parent = parent
      self.logger = logging.getLogger(self.parent.name+".Channel")
      self.name = name
      self.data['freq'] = freq
      self.data['bw'] = bw
      self.data['pol'] = pol 

class Map(object):
  """
  Class for all the data and methods associated with a raster scan map
  
  It is expected that the parent class is a subclass of ``Observation`` already
  by virtue of it being a superclass of subclass which inherits these methods.
  
  Attrs:
    cfg (dict): 
    data (numpy array):      from ``Observation``
    logger (logging.Logger): replaces ``Observation`` logger
    name (str):              replaces ``Observation`` name
    session (Session):
    source (str):
    step (float):            map step size
    
  """  
  def regrid(self, width=1.0, height=1.0, step=None):
    """
    converts a map from observed coordinates to map coordinates
    
    @param width : map width in deg
    @type  width : float
    
    @param height : map height in deg
    @type  height : float
    
    @param step : map step size in X and Y in deg
    @type  step : (float, float)
    """
    if step == None:
      # use the original step size
      self.xstep = abs(self.data['xdec_offset'][1:]
                      -self.data['xdec_offset'][:-1]).mean().round(6)
      self.ystep = abs(self.data['dec_offset'][1:]
                      -self.data['dec_offset'][:-1]).mean().round(6)
    else:
      self.xstep, self.ystep = step
    nx = int(round(width/self.xstep))
    ny = int(round(height/self.ystep))
    self.data['grid_x'] = NP.arange( -width/2,  width/2 + self.xstep, self.xstep/2)
    self.data['grid_y'] = NP.arange(-height/2, height/2 + self.ystep, self.ystep/2)
    self.data['grid_z'] = {}
    for chnl in self.channel:
      cx = self.data['xdec_offset']
      cy = self.data['dec_offset']
      cz = self.data[chnl]
      xi = self.data['grid_x']
      yi = self.data['grid_y']
      try:
        self.data['grid_z'] = griddata(cx,cy,cz, xi, yi, method='nearest')
      except ValueError as details:
        self.logger.error("regrid: gridding failed: %s", str(details))
        self.logger.debug("regrid: channel %d length of cx is %d", chnl, len(cx))
        self.logger.debug("regrid: channel %d length of cy is %d", chnl, len(cy))
        self.logger.debug("regrid: channel %d length of cz is %d", chnl, len(cz))
        continue
    return self.data['grid_x'], self.data['grid_y'], \
           self.data['grid_z']

# ------------------------ module functions -------------------------------

def get_obs_session(project=None, dss=None, date=None, path='proj'):
  """
  Provides project, station, year and DOY, asking as needed.
  
  It follows one of several possible paths to get to the session::
  
    proj - path through /usr/local/projects/project
    hdf5 - path through /usr/local/RA_data/HDF5
    fits - path through /usr/local/RA_data/FITS
    wvsr - path through /data
  
  @param project : optional name as defined in /usr/local/projects
  @type  project : str
  
  @param dss : optional station number
  @type  dss : int
  
  @param date : optional YYYY/DDD
  @type  date : str
  
  @return: project, DSS, year, DOY.
  """
  def get_directory(path):
    """
    """
    # only one trailing /
    path = path.rstrip('/')+"/*"
    logger.debug("get_obs_session:get_directory: from %s", path)
    names = glob.glob(path)
    if names:
      dirs = []
      for name in names:
        if os.path.isdir(name):
          dirs.append(os.path.basename(name))
      dirs.sort()
      for name in dirs:
        print((name), end=' ')
      return input('\n>')
    else:
      return []
      
  def from_wvsr_dir():
    """
    this needs to be completed and tested on crab14 or an auto host
    """
    session = get_directory(wvsr_dir)
    return session
    
  cwd = os.getcwd()
  # get the project
  if project:
    pass
  else:
    os.chdir(projects_dir)
    project = get_directory(projects_dir)
  projectpath = projects_dir+project
  # get the station
  if path[:4].lower() == 'wvsr':
    # special call
    print("from_wvsr_dir()")
  if path[:4].lower() == 'proj':
    os.chdir(projectpath+"/Observations/")
  elif path[:4].lower() == 'hdf5':
    os.chdir(hdf5_dir)
  elif path[:4].lower() == 'fits':
    os.chdir(fits_dir)
  
  # get the station
  if dss:
    pass
  else:
    # This seems odd but get_directory() needs '/' and int does not
    station = get_directory(os.getcwd()+"/").rstrip('/')
    dss = int(station[-2:])
  stationpath = os.getcwd()+"/dss"+str(dss)
  # get the date
  if date:
    items = date.split('/')
    year = int(items[0])
    DOY = int(items[1])
  else:
    year = int(get_directory(stationpath))
    yearpath = stationpath+"/"+str(year)
    DOY = int(get_directory(yearpath))
  os.chdir(cwd)
  return project, dss, year, DOY  

def get_obs_dirs(project, station, year, DOY, datafmt=None):
  """
  Returns the directories where data and working files are kept
  
  @param project : project code string, e.g., RRL
  @type  project : str
  
  @param station : DSN station number
  @type  station : int
  
  @param year : year of observation
  @type  year : int
  
  @param DOY : day of year of observations
  @type  DOY : int
  
  @param datafmt : raw data format
  @type  datafmt : str
  """
  logger.debug("get_obs_dirs: type %s for %s, DSS%d, %4d/%03d",
               datafmt, project, station, year, DOY)
  obspath = "dss%2d/%4d/%03d/" %  (station,year,DOY)
  if project:
    projdatapath = "/usr/local/project_data/"+project+"/"+obspath
    projworkpath = "/usr/local/projects/"+project+"/Observations/"+obspath
  else:
    projdatapath = ""
    projworkpath = ""
  if datafmt:
    rawdatapath = "/usr/local/RA_data/"+datafmt+"/"+obspath
  else:
    rawdatapath = ""
  return projdatapath, projworkpath, rawdatapath

def select_data_files(datapath, name_pattern="", load_hdf=False):
  """
  Provide the user with menu to select data files.
  
  Finding the right data store is complicated. As originally coded::
  
    * If the input files are ``.h5`` the ``load_hdf=True`` and the directory area is
      RA_data/HDF5/.
    * If the input files are .pkl (obsolete) then load_hdf=False and the
      directory area is project_data/<project>/.
      
  Now we need to add the possibility of getting datafiles from RA_data/FITS/.
  Now the location of the data are implicit in ``datapath``::
  
    * If datapath is ...RA_data/HDF5/... then the files could be .h5 (Ashish) or .hdf5 (Dean).
    * If datapath is ...RA_data/FITS/... then the extent is .fits.
    * If datapath is ...project_data... then the extent is .pkl
  
  @param datapath : path to top of the tree where the DSS subdirectories are
  @type  datapath : str
  
  @param name_pattern : pattern for selecting file names, e.g. source
  @type  name_pattern : str
  
  @param load_hdf : use RA_data/HDF5 directory if True
  @type  load_hdf : bool
  
  @return: list of str
  """
  # Get the data files to be processed
  logger.debug("select_data_files: looking in %s", datapath)
  if name_pattern:
    name,extent = os.path.splitext(name_pattern)
    if extent.isalpha(): # a proper extent with no wildcards
      # take name pattern as is
      pass
    else:
      # only one * at front and back of pattern
      name_pattern = "*"+name_pattern.rstrip('*')+"*"
  else:
    # no pattern specified.  All files.
    name_pattern = "*"
  logger.debug("select_data_files: for pattern %s", name_pattern)
  if re.search('HDF5', datapath):
    load_hdf = True
  elif re.search('project_data', datapath):
    load_hdf = False
    datafiles = support.text.select_files(datapath+name_pattern+"[0-9].pkl")
  elif re.search('FITS', datapath):
    datafiles = support.text.select_files(datapath+name_pattern+".fits")
  if load_hdf:
    full = datapath+name_pattern+".h*5"
  else:
    full = datapath+name_pattern
  logger.debug("select_data_files: from: %s", full)
  datafiles = support.text.select_files(full)

  logger.debug("select_data_files: found %s", datafiles)
  if datafiles == []:
    logger.error("select_data_files: None found. Is the data directory mounted?")
    raise RuntimeError('No data files found.')  
  if type(datafiles) == str:
    datafiles = [datafiles]
  logger.info("select_data_files: to be processed: %s", datafiles)
  return datafiles

def get_num_chans(linefreq, bandwidth, max_vel_width):
  """
  compute the base 2 number of output channels for the specified resolution
  """
  kmpspMHz = 300000./linefreq
  BW_kmps = bandwidth*kmpspMHz
  est_num_chan_out = BW_kmps/max_vel_width
  logger.debug("get_num_chans: estimated num chans out = %d",
               est_num_chan_out)
  return 2**int(math.ceil(math.log(est_num_chan_out,2)))
    
def reduce_spectrum_channels(spectrum, refval, refpix, delta,
                             num_chan=1024, axis=0):
  """
  Reduce the number of channels in the spectrum.
  
  The default option is to reduce the spectrum to a specified number of
  channels with a default of 1024. The input spectrum is presumed to have
  2**N channels so that num_chan/num_chan_in is an integer.
  
  If 'spectrum' is an N-D array, then the spectrum axis is given by 'axis'
  which defaults to 0.
  
  'delta' is negative for lower sideband or reversed double sideband spectra.
    
  @param spectrum : spectrum values
  @type  spectrum : list or nparray
  
  @param refval : X-axis value at the reference pixel of 'spectrum'
  @type  refval : float
  
  @param refpix : reference pixel for 'spectrum'
  @type  refpix : int
  
  @param delta : interval between pixels on the X-axis
  @type  delta : float
  
  @param num_chan : optional number of channels to be returned (default: 2^10)
  @type  num_chan : int
  
  @return: numpy.array
  """
  if math.log(num_chan,2) % 1:
    raise RuntimeError("num_chan = %d is not a power of 2", num_chan)
  if type(spectrum) == NP.ndarray:
    num_chans_in = spectrum.shape[axis]
  else:
    num_chans_in = len(spectrum)
  if math.log(num_chans_in,2) % 1:
    raise RuntimeError("input spectrum length = %d is not a power of 2",
                                                                  num_chans_in)
  logger.debug("reduce_spectrum_channels: %d channels in", num_chans_in)
  
  num_chan_avg = num_chans_in/num_chan
  newrefpix = refpix/num_chan_avg
  logger.debug("reduce_spectrum_channels: refpix from %d to %d",
               refpix, newrefpix)

  newdelta = delta*num_chan_avg
  logger.debug("reduce_spectrum_channels: delta from %.3f to %.3f",
               delta, newdelta)
  newrefval = refval + delta*(num_chan_avg/2 - 1)
  logger.debug("reduce_spectrum_channels: refval from %.3f to %.3f",
               refval, newrefval)
  logger.debug("reduce_spectrum_channels: averaging %d channels", num_chan_avg)
  
  specout = NP.array([spectrum[index*num_chan_avg:(index+1)*num_chan_avg].mean()
                                                 for index in range(num_chan)])
  logger.debug("reduce_spectrum_channels: %d channels out", num_chan)
  return specout, newrefval, newrefpix, newdelta

def get_freq_array(bandwidth, n_chans):
  """
  Create an array of frequencies for the channels of a backend

  @param bandwidth : bandwidth
  @type  bandwidth : float

  @param n_chans : number of channels
  @type  n_chans : int

  @return: frequency of each channel in same units as bandwidth
  """
  return NP.arange(n_chans)*float(bandwidth)/n_chans

def freq_to_chan(frequency,bandwidth,n_chans):
  """
  Returns the channel number where a given frequency is to be found.

  @param frequency : frequency of channel in sane units as bandwidth.
  @type  frequency : float

  @param bandwidth : upper limit of spectrometer passband
  @type  bandwidth : float

  @param n_chans : number of channels in the spectrometer
  @type  n_chans : int

  @return: channel number (int)
  """
  if frequency < 0:
    frequency = bandwidth + frequency
  if frequency > bandwidth:
    raise RuntimeError("that frequency is too high.")
  return round(float(frequency)/bandwidth*n_chans) % n_chans

def get_smoothed_bandshape(spectrum, degree = None, poly_order=15):
  """
  Do a Gaussian smoothing of the spectrum and then fit a polynomial.
  Optionally, the raw and smoothed data and the fitted polynomial can be
  plotted.

  Note
  ====
  ``numpy.polyfit(x, y, deg, rcond=None, full=False, w=None, cov=False)``
  Least squares polynomial fit.
  Fit a polynomial::
  
     p(x) = p[0] * x**deg + ... + p[deg]
  
  of degree deg to points (x, y).
  Returns a vector of coefficients p that minimises the squared error.

  @param spectrum : input data
  @type  spectrum : list of float

  @param degree : number of samples to smoothed (Gaussian FWHM)
  @type  degree : int

  @param poly_order : order of the polynomial
  @type  poly_order : int

  @param plot : plotting option
  @type  plot : boolean

  @return: (polynomial_coefficient, smoothed_spectrum)
  """
  if degree == None:
    degree = len(spectrum)/100
  # normalize the spectrum so max is 1 and convert to dB.
  max_lev = NP.max(spectrum)
  norm_spec = NP.array(spectrum)/float(max_lev)
  norm_spec_db = 10*NP.log10(norm_spec)
  # optionally plot normalized spectrum
  if plot:
    pylab.plot(norm_spec_db)
  # do a Gaussian smoothing
  norm_spec_db_smoothed = smoothListGaussian(norm_spec_db, degree=degree)
  # deal with the edges by making them equal to the smoothed end points
  norm_spec_db_smoothed_resized = NP.ones(len(spectrum))
  # left end
  norm_spec_db_smoothed_resized[0:degree] = norm_spec_db_smoothed[0]
  # middle
  norm_spec_db_smoothed_resized[degree:degree+len(norm_spec_db_smoothed)] = \
      norm_spec_db_smoothed
  # right end
  norm_spec_db_smoothed_resized[degree+len(norm_spec_db_smoothed):] = \
      norm_spec_db_smoothed[-1]
  if plot:
    pylab.plot(norm_spec_db_smoothed_resized)
    poly = NP.polyfit(list(range(len(norm_spec_db_smoothed))),
                         norm_spec_db_smoothed,poly_order)
    pylab.plot(NP.polyval(poly, list(range(len(norm_spec_db_smoothed)))))
    pylab.show()
  return poly, norm_spec_db_smoothed_resized
