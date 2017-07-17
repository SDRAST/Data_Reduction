# -*- coding: utf-8 -*-
"""
Modules to support data reduction in Python.

A lot of the modules in DSN need work because of the change from Observatory
to MonitorControl.
"""
# standard Python modules
import logging
import datetime
import scipy.fftpack

# import from standard Python modules
from os.path import basename
from numpy import array, loadtxt

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

def get_obs_dirs(project, station, year, DOY, datafmt="SDFITS"):
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
  rawdatapath = "/usr/local/RA_data/"+datafmt+"/"+obspath
  return projdatapath, projworkpath, rawdatapath

def select_data_files(datapath, name_pattern="", load_hdf=False):
  """
  Provide the user with menu to select data files.
  
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
    # only one * at ont and back of pattern
    name_pattern = "*"+name_pattern.strip('*')+"*"
  else:
    name_pattern = "*"
  logger.debug("select_data_files: for pattern %s", name_pattern)
  if load_hdf:
    allfiles = datapath+name_pattern+".spec.h5" + datapath+name_pattern+".hdf5"
    datafiles = select_files(allfiles)
  else:
    datafiles = select_files(datapath+name_pattern+"[0-9].pkl")
  if datafiles == []:
    logger.error("select_data_files: Is the data directory mounted?")
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
