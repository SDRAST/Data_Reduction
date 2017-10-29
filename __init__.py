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

from scipy import array, argmin, argmax, mean, sqrt, var
from support import nearest_index # for older code
from support.text import select_files

logger = logging.getLogger(__name__)

def get_obs_session(project=None, dss=None, date=None):
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
  compute the base 2 number of output channels for the specified resolution
  """
  kmpspMHz = 300000./linefreq
  BW_kmps = bandwidth*kmpspMHz
  est_num_chan_out = BW_kmps/max_vel_width
  logger.debug("get_num_chans: estimated num chans out = %d",
               est_num_chan_out)
  return 2**int(math.ceil(math.log(est_num_chan_out,2)))
    
def reduce_spectrum_channels(spectrum, refval, refpix, delta, num_chan=1024):
  """
  Reduce the number of channels in the spectrum.
  
  The default option is to reduce the spectrum to a specified number of
  channels with a default of 1024. The input spectrum is presumed to have
  2**N channels so that num_chan/num_chan_in is an integer.
    
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
  if type(spectrum) == ndarray:
    num_chans_in = spectrum.shape[0]
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
  newrefval = refval+newdelta/2
  logger.debug("reduce_spectrum_channels: refval from %.3f to %.3f",
               refval, newrefval)
  logger.debug("reduce_spectrum_channels: averaging %d channels", num_chan_avg)
  
  specout = array([spectrum[index*num_chan_avg:(index+1)*num_chan_avg].mean()
                                                 for index in range(num_chan)])
  logger.debug("reduce_spectrum_channels: %d channels out", num_chan)
  return specout, newrefval, newrefpix, newdelta
  
def trim_extremes(data):
  """
  Remove extreme values from a data array.

  Extreme values are those greater than 10x the standard deviation
  and are 'skinny'.: numpy array

  @param data : numpy array
  """
  data_array = array(data)
  amean   = mean(data_array)
  avar    = var(data_array)
  astdev  = sqrt(avar)
  amax = data_array.max()
  amin = data_array.min()
  # Check the maximum
  if abs(amax-amean) > 10*astdev:
    index = argmax(data_array)
    if is_skinny(data_array,index):
      data_array = clobber(data_array,index)
  # check the minimum
  if abs(amin-amean) > 10*astdev:
    index = argmin(data_array)
    if is_skinny(data_array,index):
      data_array = clobber(data_array,index)
  return data_array

def is_skinny(data_array,index):
  """
  Test whether a data value is an isolated outlier

  Returns True if the data values adjacent to the test value are
  less that 1/10 of the test value, i.e., the data point is a spike

  @param data_array : numpy array

  @param index : int

  @return: boolean
  """
  amean   = mean(data_array)
  test_value = abs(data_array[index]-amean)
  if index == 0:
    ref_value = abs(data_array[index+1] - amean)
  elif index == len(data_array)-1:
    ref_value = abs(data_array[index-1] - amean)
  else:
    ref_value = (data_array[index-1] + data_array[index+1])/2. - amean
  if test_value > 10 * ref_value:
    return True
  else:
    return False

def clobber(data_array,index):
  """
  Replace a data array value with the adjacent value(s)

  @param data_array : numpy array

  @param index : int
  """
  if index == 0:
    data_array[index] = data_array[index+1]
  elif index == len(data_array)-1:
    data_array[index] = data_array[index-1]
  else:
    data_array[index] = (data_array[index-1] + data_array[index+1])/2.
  return data_array

