# -*- coding: utf-8 -*-
"""
Modules to support data reduction in Python.

A lot of the modules in DSN need work because of the change from Observatory
to MonitorControl.
"""
# standard Python modules
import datetime
import logging
import math
import re
import scipy.fftpack

from numpy import array, loadtxt, ndarray
from os.path import basename, splitext

from support import nearest_index # for older code
from support.text import select_files

logger = logging.getLogger(__name__)

def get_obs_session(datafmt=None, project=None, dss=None, date=None):
  """
  Asks user for parameters to locate observation session paths
  
  This expects two directory trees to exist.  For raw data::
    /usr/local/RA_data/HDF5/
      dssXX/
        YEAR/
          DOY/
  and for pickled data::
    /usr/local/project_data/
      dssXX/
        YEAR/
          DOY/
  
  @param raw : if True data will come from /usr/local/RA_data/
  @type  raw : bool
  
  @param datafmt : optional "HDF5" or "SDFITS" or ... for sub-directory
  @type  datafmt : str
  
  @param project : optional name as defined in /usr/local/projects
  @type  project : str
  
  @param dss : optional station number
  @type  dss : int
  
  @param date : optional YYYY/DDD
  @type  date : str
  
  @return: project, DSS, year, DOY.
  """
  # get the path to the project directory
  if project:
    projectpath = "/usr/local/projects/"+project+"/"
  else:
    projectpath = select_files("/usr/local/projects/*",
                               text="Select a project by index: ", single=True)
    project = basename(projectpath)
    if projectpath[-1] != "/":
      projectpath += "/"
  logger.debug("get_obs_session: project path: %s", projectpath)
  # get the path to the raw data
  rawdatapath = "/usr/local/RA_data/"
  if datafmt:
    rawfmtpath = "/usr/local/RA_data/"+datafmt
  else:
    rawfmtpath = select_files(rawdatapath+"[A-Z]*",
                           text="Select a data format by index: ", single=True)
  rawfmt = basename(rawfmtpath)
  rawfmtpath = rawfmtpath+"/"
  logger.debug("get_obs_session: raw data path: %s", rawfmtpath)
  # get the path to the project DSS sub-directory
  if dss:
    dsspath = projectpath+"Observations/dss"+str(dss)+"/"
  else:
    dsspath = select_files(projectpath+"Observations/dss*",
                           text="Select a station by index: ", single=True)
    logger.debug("get_obs_session: selected: %s", dsspath)
    dss = int(basename(dsspath)[-2:])
  logger.debug("get_obs_session: DSS path: %s", dsspath)
  if date:
    items = date.split('/')
    yr = int(items[0])
    doy = int(items[1])
  else:
    yrpath = select_files(dsspath+"/20*",
                                  text="Select a year by index: ", single=True)
    if yrpath:
      logger.debug("get_obs_session: year path: %s", yrpath)
      yr = int(basename(yrpath))
      yrpath += "/"
      doypath = select_files(yrpath+"/*",
                                   text="Select a day BY INDEX: ", single=True)
      doy = int(basename(doypath))
      doypath += '/'
      logger.debug("get_obs_session: DOY path: %s", doypath)
    else:
      logger.warning("get_obs_session: no data for dss%2d", dss)
      return project, None, 0, 0, None
  logger.debug("get_obs_session: for %s, DSS%d, %4d/%03d, raw data is %s",
                    project, dss, yr, doy, rawfmt)
  return project, dss, yr, doy, datafmt

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
  projdatapath = "/usr/local/project_data/"+project+"/"+obspath
  projworkpath = "/usr/local/projects/"+project+"/Observations/"+obspath
  if datafmt:
    rawdatapath = "/usr/local/RA_data/"+datafmt+"/"+obspath
  else:
    rawdatapath = ""
  return projdatapath, projworkpath, rawdatapath

def select_data_files(datapath, name_pattern="", load_hdf=False):
  """
  Provide the user with menu to select data files.
  
  Finding the right data store is complicated. As originally coded::
    * If the input files are .h5 the load_hdf=True and the directory area is
      RA_data/HDF5/.
    * If the input files are .pkl (obsolete) then load_hdf=False and the
      directory area is project_data/<project>/.
  Now we need to add the possibility of getting datafiles from RA_data/FITS/.
  Now the location of the data are implicit in 'datapath'::
    * If datapath is ...RA_data/HDF5/... then the files could be .h5 (Ashish)
      or .hdf5 (Dean).
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
    name,extent = splitext(name_pattern)
    if extent.isalpha(): # a proper extent with no wildcards
      # take name pattern as is
      pass
    else:
      # only one * at front and back of pattern
      name_pattern = "*"+name_pattern.strip('*')+"*"
  else:
    # no pattern specified.  All files.
    name_pattern = "*"
  logger.debug("select_data_files: for pattern %s", name_pattern)
  if re.search('HDF5', datapath):
    load_hdf = True
  elif re.search('project_data', datapath):
    load_hdf = False
    datafiles = select_files(datapath+name_pattern+"[0-9].pkl")
  elif re.search('FITS', datapath):
    datafiles = select_files(datapath+name_pattern+".fits")
  if load_hdf:
    full = datapath+name_pattern+".h*5"
  else:
    full = datapath+name_pattern
  logger.debug("select_data_files: from: %s", full)
  datafiles = select_files(full)

  logger.debug("select_data_files: found %s", datafiles)
  if datafiles == []:
    logger.error("select_data_files: None found. Is the data directory mounted?")
    raise RuntimeError('No data files found.')  
  if type(datafiles) == str:
    datafiles = [datafiles]
  logger.info("select_data_files: to be processed: %s", datafiles)
  return datafiles

def load_csv_with_header(filename,delimiter=""):
  """
  Takes a text table with headers and converts it into an ASCII array

  There must be a header for each column and the separator must not be
  ambiguous, for example, a space if the column headers also include spaces.

  The data are ASCII and the program or programmer will have to convert
  them to an appropriate form, e.g. int or float or datetime, if necessary.

  A column can be extracted using data[label].

  @param filename : name of the data file
  @type  filename : str

  @param delimiter : column separator
  @type  delimiter : str

  @return: column headers (str), data array(str)
  """
  datafile = open(filename,"r")
  labels = datafile.readline().strip().split(',')
  num_cols = len(labels)
  fmts = ('S20',)*num_cols
  datafile.close()
  data = loadtxt(filename,skiprows=1,delimiter=',',
                 dtype = {'names': tuple(labels),
                          'formats':fmts})
  return labels,data

def unpack_to_complex(rawdata):
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

def sideband_separate(data):
  """
  Converts a complex array time series and returns two reals with USB and LSB

  This applies a Hilbert transform to the complex data.
  """
  usb = (data.real + scipy.fftpack.hilbert(data).imag)
  lsb = (scipy.fftpack.hilbert(data).real + data.imag)
  return lsb,usb

def get_num_chans(linefreq, bandwidth, max_vel_width):
  """
  compute number of output channels for the specified resolution
  """
  kmpspMHz = 300000./linefreq
  BW_kmps = bandwidth*kmpspMHz
  est_num_chan_out = BW_kmps/max_vel_width
  logger.debug("get_num_chans: estimated num chans out = %d",
               est_num_chan_out)
  return 2**int(math.ceil(math.log(est_num_chan_out,2)))
    

def reduce_spectrum_channels(spectrum, num_chan=1024,
                            linefreq=None, bandwidth=None, max_vel_width=None):
  """
  Reduce the number of channels in the spectrum.
  
  The default option is to reduce the spectrum to a specified number of
  channels with a default of 1024. The input spectrum is presumed to have
  2**N channels so that num_chan/num_chan_in is an integer. If linefreq,
  bandwidth and max_vel_width are given, the number of channels is computed
  and overrides the argument num_chan.
  
  @param spectrum : spectrum values
  @type  spectrum : list or nparray
  
  @param num_chan : optional number of channels to be returned (default: 2^10)
  @type  num_chan : int
  
  @param linefreq : optional line frequency in MHz
  @type  linefreq : float or int
  
  @param bandwidth : optional width of the spectrum in MHz
  @type  bandwidth : float or int
  
  @param max_vel_width : optional maximum channel width in km/s
  @type  max_vel_width : float
  """
  if type(spectrum) == ndarray:
    num_chans_in = spectrum.shape[0]
  else:
    num_chans_in = len(spectrum)
  logger.debug("reduce_spectrum_channels: %d channels in", num_chans_in)
  logger.debug("reduce_spectrum_channels: input: %s", spectrum)
  if num_chan:
    num_chans = num_chan
  elif linefreq and bandwidth and max_vel_width:
    num_chans = get_num_chans(linefreq, bandwidth, max_vel_width)
  else:
    logger.error("reduce_spectrum_channels: could not get number of channels")
    return None
  logger.debug("reduce_spectrum_channels: %d channels out", num_chan)
  num_chan_avg = num_chans_in/num_chan
  logger.debug("reduce_spectrum_channels: averaging %d channels", num_chan_avg)
  if num_chan_avg:
    specout = array([spectrum[index*num_chan_avg:(index+1)*num_chan_avg].mean()
                                                 for index in range(num_chan)])
    logger.debug("reduce_spectrum_channels: new values: %s", specout)
    return specout
  else:
    return None
