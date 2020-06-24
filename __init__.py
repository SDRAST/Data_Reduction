# -*- coding: utf-8 -*-
"""
Modules to support data reduction in Python.

Note
====
These provides the base class ``Observation``. IN PROGRESS!  The definition
of the keywords is done.
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

import Astronomy as A
import DatesTimes as DT
from local_dirs import fits_dir, hdf5_dir, projects_dir, wvsr_dir
import support

# enable raw_input Tab completion
readline.parse_and_bind("tab: complete")

module_logger = logging.getLogger(__name__) # module logger

class Observation(object):
  """
  superclass for a data structure and methods
  
  Attributes
  ==========
    channel - (dict) a signal path
    data    - (dict) table contents
  
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
    Top        (float) alternative for ``tsys``
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
              'timestr','chan','Tsys','top','diode','level','cryotemp',
              'windspeed','winddir','ambtemp','pressure',
              'ra','dec','ra1950','dec1950','ra2000','dec2000']
  def __init__(self, channel_names=None, 
                     longitude=None,
                     latitude=None):
    """
    Arguments
    =========
    channel_names:string: like "XL"
    """
    self.logger = logging.getLogger(module_logger.name+".Observation")
    if longitude and latitude:
      self.longitude = longitude
      self.latitude = latitude
    else:
      self.logger.error("__init__: this requires observatory location")
      raise Exception("Where were the data taken?")
    self.channel = {}
    if channel_names:
      for ch in channel_names:
        self.channel[ch] = Observation.Channel(ch)

  def load_text_with_header(self, filename, delimiter=" ",
                            channels=None, ignore=[], skip_header=0,
                            names=True, aliases={"epoch": "unixtime"}):
    """
    Takes a text table with headers and converts it into a structured numpy
    array. That means that a column can be extracted using `data[label]`.

    Ideally, there must be a header word for each column and the separator must
    not be  ambiguous, for example, a space if the column headers also include
    spaces.  However, the ``numpy`` function ``genfromtxt()`` is very flexible.
    Here is an example of reading DSS-28 text files from solar observations in
    2012, 9 columns but 7 header words, when database access was not yet 
    practical::
    
      In [1]: from Data_Reduction import Observation
      In [2]: from Astronomy.DSN_coordinates import DSS
      In [3]: import logging
      In [4]: tel = DSS(28)
      In [5]: obs = Observation(longitude=tel.long, latitude=tel.lat)
      In [6]: obs.logger.setLevel(logging.INFO)
      In [7]: obs.load_text_with_header('t12127.10',
                      delimiter=[17,16,3,11,7,9,8,2,6],
                      names="UTC Epoch Chan Tsys Int Az El Diode Level".split(),
                      skip_header=1)

    The ``names`` string was just copied by hand from the first line in the
    file. It could also be done with ``file.readline()``.
    
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
    data = NP.genfromtxt(filename, delimiter=delimiter, dtype=None, names=names,
                         case_sensitive='lower', skip_header=skip_header)
    numdata = len(data)
    self.logger.debug("load_text_with_header: %d samples", numdata)
    names = data.dtype.names
    self.logger.debug("load_text_with_header: names type is %s", type(names))
    self.logger.info("load_text_with_header: names = %s", names)
    # get the known columns:
    self.data = {}
    self.unknown_cols = []
    for name in names:
        if name.casefold() in map(str.casefold, aliases):
            key = aliases[name].lower() # we use only lower case names
        else:
            key = name.lower()
        if key in map(str.casefold, Observation.reserved):
            self.data[key] = data[name]
        elif name in ignore:
            pass
        else:
            self.unknown_cols.append(name)
    for name in self.unknown_cols:
        self.data[name] = data[name]
    # get columns that are better computed once instead of as needed
    colnames = list(self.data.keys())
    if 'unixtime' in colnames:
        self.data['datetime'] = []
        self.data['date_num'] = []
        for idx in list(range(numdata)):
            tm = data[next(key for key, value in aliases.items() if value == 'unixtime')][idx]
            dt = datetime.datetime.utcfromtimestamp(tm)
            self.data['datetime'].append(dt)
            self.data['date_num'].append(MPLd.date2num(dt))
    if self.check_for(data, 'azel'):
      # azel exists; compute radec if needed; then radec2000 if needed
      if self.check_for(data, 'radec'):
        pass
      else:
        # compute RA and dec
        RA = []; decs = []; RAdecs = []
        for idx in list(range(numdata)):
          # setup
          dt = self.data['datetime'][idx]
          time_tuple = (dt.year,
                        DT.day_of_year(dt.year,dt.month,dt.day)
                        + (  dt.hour
                           + dt.minute/60.
                           + dt.second/3600.
                           + dt.microsecond/3600./1e6)/24.)
          azimuth = data['az'][idx]
          elevation = data['el'][idx]
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
        if self.check_for(data, 'radec2000'):
          # ra2000 and dec2000 already exist
          pass
        else:
          # compute RA2000 and dec2000 from observed RA and dec
          RA2000 = []; decs2000 = []; RAdec2000 = []
          for idx in list(range(numdata)):
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
      # finished with azel -> radec -> radec2000
    elif self.check_for(data, 'radec2000'):
      # coordinates exist; compute back to azimuth and elevation
      if self.check_for(data, 'radec'):
        pass
      else:
        # compute observed RA and dec
        RA = []; decs = []; RAdecs = []
        for idx in list(range(numdata)):
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
        if self.check_for(data, 'azel'):
          pass
        else:
          # compute az and el
          azs = []; els = []; azels = []
          for idx in list(range(numdata)):
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
    # in here check for 'radec'
    else:
      self.logger.error("no coordinates found in data")
      raise Exception("check INFO logging for columns found")
      
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
    """
    longitude, latitude, epoch = self.splitkey(longlat)
    if longitude in data.dtype.names and \
       latitude  in data.dtype.names:
      self.logger.debug("check_for: data has %s and %s", longitude, latitude)
      return True
    elif longlat in data.dtype.names:
      self.logger.debug("check_for: data has %s", longlat)
      data[longitude],data[latitude] = map(None, *data[longlat])
      self.logger.debug("check_for: added %s and %s to data", longitude, latitude)
      return True
    else:
      # coords need to be computed from other coords
      
      return False
    
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
    def __init__(self, name, freq=None, bw=None, pol=None):
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
      self.name = name
      self.data['freq'] = freq
      self.data['bw'] = bw
      self.data['pol'] = pol 
  
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
    module_logger.debug("get_obs_session:get_directory: from %s", path)
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
  
  
def get_obs_session_old(project=None, dss=None, date=None):
  """
  Asks user for parameters to locate observation session paths
  
  This expects one of two directory trees to exist.  If dtype is given::
    /usr/local/RA_data/dtype/
      dssXX/
        YEAR/
          DOY/
  
  or if project is given::
    /usr/local/projects/project/Observations/
      dssXX/
        YEAR/
          DOY/
          
  If neither is given it will prompt for a project.
  
  @param project : optional name as defined in /usr/local/projects
  @type  project : str
  
  @param dss : optional station number
  @type  dss : int
  
  @param date : optional YYYY/DDD
  @type  date : str
  
  @return: project, DSS, year, DOY.
  """
  # get the path to the session directory
  if project:
    projectpath = projects_dir+project+"/"
  else:
    projectpath = support.text.select_files(projects_dir+"*", ftype="dir",
                               text="Select a project by index: ", single=True)
    project = os.path.basename(projectpath)
    if projectpath[-1] != "/":
      projectpath += "/"
  module_logger.debug("get_obs_session: project path: %s", projectpath)

  # get the path to the project DSS sub-directory
  project_obs_path = projectpath + "Observations/"
  if dss:
    dsspath = project_obs_path+"dss"+str(dss)+"/"
  else:
    dsspath = support.text.select_files(project_obs_path+"dss*", ftype="dir",
                           text="Select a station by index: ", single=True)
    module_logger.debug("get_obs_session: selected: %s", dsspath)
    dss = int(os.path.basename(dsspath)[-2:])
  module_logger.debug("get_obs_session: DSS path: %s", dsspath)
  if date:
    items = date.split('/')
    yr = int(items[0])
    doy = int(items[1])
  else:
    yrpath = support.text.select_files(dsspath+"*", ftype="dir",
                                  text="Select a year by index: ", single=True)
    if yrpath:
      module_logger.debug("get_obs_session: year path: %s", yrpath)
      yr = int(os.path.basename(yrpath))
      yrpath += "/"
      doypath = support.text.select_files(yrpath+"/*", ftype="dir",
                                   text="Select a day BY INDEX: ", single=True)
      doy = int(os.path.basename(doypath))
      doypath += '/'
      module_logger.debug("get_obs_session: DOY path: %s", doypath)
    else:
      module_logger.warning("get_obs_session: no data for dss%2d", dss)
      return project, None, 0, 0
  module_logger.debug("get_obs_session: for %s, DSS%d, %4d/%03d",
                    project, dss, yr, doy)
  return project, dss, yr, doy

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
  module_logger.debug("get_obs_dirs: type %s for %s, DSS%d, %4d/%03d",
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
  module_logger.debug("select_data_files: looking in %s", datapath)
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
  module_logger.debug("select_data_files: for pattern %s", name_pattern)
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
  module_logger.debug("select_data_files: from: %s", full)
  datafiles = support.text.select_files(full)

  module_logger.debug("select_data_files: found %s", datafiles)
  if datafiles == []:
    module_logger.error("select_data_files: None found. Is the data directory mounted?")
    raise RuntimeError('No data files found.')  
  if type(datafiles) == str:
    datafiles = [datafiles]
  module_logger.info("select_data_files: to be processed: %s", datafiles)
  return datafiles

def get_num_chans(linefreq, bandwidth, max_vel_width):
  """
  compute the base 2 number of output channels for the specified resolution
  """
  kmpspMHz = 300000./linefreq
  BW_kmps = bandwidth*kmpspMHz
  est_num_chan_out = BW_kmps/max_vel_width
  module_logger.debug("get_num_chans: estimated num chans out = %d",
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
  module_logger.debug("reduce_spectrum_channels: %d channels in", num_chans_in)
  
  num_chan_avg = num_chans_in/num_chan
  newrefpix = refpix/num_chan_avg
  module_logger.debug("reduce_spectrum_channels: refpix from %d to %d",
               refpix, newrefpix)

  newdelta = delta*num_chan_avg
  module_logger.debug("reduce_spectrum_channels: delta from %.3f to %.3f",
               delta, newdelta)
  newrefval = refval + delta*(num_chan_avg/2 - 1)
  module_logger.debug("reduce_spectrum_channels: refval from %.3f to %.3f",
               refval, newrefval)
  module_logger.debug("reduce_spectrum_channels: averaging %d channels", num_chan_avg)
  
  specout = NP.array([spectrum[index*num_chan_avg:(index+1)*num_chan_avg].mean()
                                                 for index in range(num_chan)])
  module_logger.debug("reduce_spectrum_channels: %d channels out", num_chan)
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
