"""
module for reducing SAO spectrometer data
"""
import cPickle # import dill as cPickle
import logging
import re

from numpy import array, ndarray, zeros
from os.path import basename, exists

from DatesTimes import module_logger as dtl
dtl.setLevel(logging.WARNING)
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from MonitorControl.BackEnds.ROACH1.SAOspec import SAOhdf5
from support import mkdir_if_needed, nearest_index

logger = logging.getLogger(__name__)

class Data(object):
  """
  A subset of data extracted from an average difference spectrum
  
  Attributes::
    channels -
    frame    -
    logger   -
    x        -
    y        -
  """
  def __init__(self, x, y):
    """
    @param x : X-values
    @type  x : list or nparray
    """
    self.logger = logging.getLogger(logger.name+".Data")
    if type(x) == list:
      self.x = array(x)
    elif type(x) == ndarray:
      self.x = x
    else:
      self.logger.error("__init__: type %s is not valid for X", type(x))
    if type(y) == dict:
      if y.keys() == [0, 1]:
        self.y = y
      else:
        self.logger.error("__init__: pols keys must be [0,1]")
    elif type(y) == list:
      if type(y[0]) == list:
        self.y = array(y)
      else:
        self.logger.error("__init__: %s is not a dict of valid Y-values")
    self.frame = None
    self.channels = None

class SAOexaminer(object):
  """
  Class for examining SAOdataset objects
  """
  def __init__(self, project=None, dss=None, date=None):
    """
    Create an SAOexaminer object
    
    Attributes::
      avg_diff     - average of all scans for each pol
      datafiles    - list of pkl files in datapath
      datapath     - location of project data for this observation session
      dataset      - collection of scana and metadata taken from HDF5 file
      logger       - logging.Logger object
      projworkpath - obs. session sub-dir. in project directory
      rawdatapath  - local of HDF5 files for this obs. session
      tau          - total integration times for averaged spectra
    
    Methods::
      align_spectrum
      beam_switched_average - average spectra for all the scan
      get_datasets          - get datasets from cPickle files
      get_obs_dirs          - get the directories used in this obs. session
      load_HDF5_files       - get datasets from HDF5 files
      load_ephemeris        - get the ephemeris for the designated source
      plot_Tsys             -
      plot_spectra          - plots avg_diff
      save_average          -
      smoothed              -
    
    @param datafiles : optional list of cPickle data files
    @type  datafiles : list of str
    """
    self.logger = logging.getLogger(logger.name+'.SAOexaminer')
    self.avg_diff = {}
    self.datafiles = None
    self.dataset = {}
    self.tau = {}      # cumulative integration time
    try:
      self.logger.debug("__init__: getting pickle files for %s at DSS-%s on %s",
                        project, dss, date)
      self.get_datasets(project=project, dss=dss, date=date)
    except TypeError, details:
      if re.search("NoneType", str(details)):
        self.logger.debug("__init__: getting HDF5 files instead")
        self.load_HDF5_files(project=project, dss=dss, date=date)
      else:
        self.logger.error("__init__: failed with TypeError: %s", details)
        raise RuntimeError("__init__: failed with TypeError: %s", details)
    except RuntimeError, details:
      if re.search("No data files found", str(details)):
        self.logger.debug("__init__: getting HDF5 files")
        self.load_HDF5_files(project=project, dss=dss, date=date)
      else:
        self.logger.error("__init__: failed with RuntimeError: %s", details)
        raise RuntimeError("__init__: failed with RuntimeError: %s", details)
    self.beam_switched_average()
    self.frame = None
    
  def get_obs_dirs(self, project=None, dss=None, date=None, rawfmt="HDF5"):
    """
    Gets the directories used for an observing session.
    
    If the session details are not provided, a menu-based inquiry is made
    
    The attribute rawdatapath is also computed but not returned.
    
    @param project : the project for which the observations were made
    @type  project : str
    
    @param dss : DSN station number
    @type  dss : int
    
    @param date : year and day of year as YYYY/DDD
    @type  data : str
    
    @param raw : get HDF5 data from RA_data store if True
    @type  raw : bool
    """
    self.datapath, self.projworkpath, self.rawdatapath = get_obs_dirs(
         *get_obs_session(dss=dss, project=project, date=date, datafmt=rawfmt))
    self.logger.debug("get_obs_dirs: found %s, %s, %s",
                            self.datapath, self.projworkpath, self.rawdatapath)
    return self.datapath, self.projworkpath, self.rawdatapath
  
  def load_ephemeris(self, sourcename):
    """
    Loads the ephemeris and returns a spline for the radial velocity with date
    """
    # the source for which the ephemeris was calculated
    if sourcename == "67P":
      # only works in /home/kuiper/Projects/Comets/67P for now
      import imp
      ephemeris_67P = imp.load_source('ephemeris_67P',
                                    '/usr/local/projects/67P/ephemeris_67P.py')
      from ephemeris_67P import get_ephemeris
      return get_ephemeris()
    else:
      self.logger.error("load_ephemeris: only have ephemeris for Comet 67P")
      raise RuntimeError("No ephemeris for %s" % opts.sourcename)

  def get_datasets(self, project=None, dss=None, date=None):
    """
    Loads selected datasets and defines additional attributes.
    
    @param project : the project for which the observations were made
    @type  project : str
    
    @param dss : DSN station number
    @type  dss : int
    
    @param date : year and day of year as YYYY/DDD
    @type  data : str
    """
    # get the datafiles to be processed
    self.logger.debug("get_datasets: for %s at %s on %s", project, dss, date)
    dirs = self.get_obs_dirs(project=project, dss=dss, date=date)
    self.logger.debug("get_datasets: dirs are %s", dirs)
    self.datafiles = select_data_files(self.datapath)
    if self.datafiles:
      self.logger.debug("get_datasets: pickle files found")
      for datafile in self.datafiles:
        self.logger.info("get_datasets: processing %s", datafile)
        index = self.datafiles.index(datafile)
        self.savebasename = self.projworkpath+basename(datafile).split('.')[0]
        self.logger.info(" will save work as %s...", self.savebasename)
        # get the data
        fd = open(datafile, "rb")
        self.dataset[index] = cPickle.load(fd)
        fd.close()
        self.dataset[index].file    = datafile
        datapath_parts = datafile.split('/')
        self.dataset[index].project = datapath_parts[-5]
        self.dataset[index].dss     = int(datapath_parts[-4][-2:])
        self.dataset[index].year    = int(datapath_parts[-3])
        self.dataset[index].doy     = int(datapath_parts[-2])
        self.dataset[index].source  = self.dataset[index].header['OBJECT'][0]
      return len(self.dataset)
    else:
      self.logger.error("get_datasets: getting HDF5 files")
      return self.load_HDF5_files(project=project, dss=dss, date=date)

  def load_HDF5_files(self, project=None, dss=None, date=None):
    """
    """
    try:
      self.datafiles = select_data_files(self.rawdatapath, load_hdf=True)
    except AttributeError:
      self.get_obs_dirs(project=project, dss=dss, date=date)
      self.datafiles = select_data_files(self.rawdatapath, load_hdf=True)
    if self.datafiles:
      for datafile in self.datafiles:
        self.logger.info("load_HDF5_files: processing %s", datafile)
        index = self.datafiles.index(datafile)
        hdf = SAOhdf5(datafile)
        self.dataset[index] = hdf.to_dataset() # loads the data into the dataset
        mkdir_if_needed(self.datapath)
        self.dataset[index].save_pickle(self.datapath)
        self.dataset[index].file    = datafile
        datapath_parts = datafile.split('/')
        self.dataset[index].project = datapath_parts[-5]
        self.dataset[index].dss     = int(datapath_parts[-4][-2:])
        self.dataset[index].year    = int(datapath_parts[-3])
        self.dataset[index].doy     = int(datapath_parts[-2])
        self.dataset[index].source  = self.dataset[index].header['OBJECT'][0]
      return len(self.dataset)
    else:
      self.logger.error("load_HDF5_files: no files found in %s",
                        self.rawdatapath)
      raise RuntimeError("load_HDF5_files: no files found")
    
  def beam_switched_average(self, index=None, mode='LINE-PBSW'):
    """
    Returns pol average spectrum and total integration of given dataset(s)
    
    Also computes the average spectrum and total integration time for each pol
    and accumulates them in the SAOexaminer attributes avg_diff and tau.
    
    @param index : dataset index
    @type  index : int
    
    @param mode : observing mode (which aught to be in the header)
    @type  mode : str
    """
    if self.dataset == {}:
      self.get_datasets()
    if index == None:
      self.indices = self.dataset.keys()
    else:
      if type(index) == int:
        self.indices = [index]
      elif type(index) == list:
        self.indices = index
      else:
        self.logger.error("beam_switched_average: invalid indices: %s", index)
    first_scan = self.dataset[0].data.keys()[0]
    num_chan = self.dataset[0].data[first_scan].shape[0]
    for dindex in self.indices: # loop over datasets
      if mode == 'LINE-PBSW':
        spec, tau = self.dataset[dindex].average_calibrated_spectrum()
        if self.avg_diff == {}:
          # cumulative average
          # xan't do it in __init__ because num_chan is not known
          self.avg_diff = {0: zeros(num_chan), 1: zeros(num_chan)}
          self.tau = {0: 0.0, 1: 0.0}      # cumulative integration time
        for pol in [0,1]:
          self.avg_diff[pol] += spec[pol]
          self.tau[pol] += tau[pol]
      else:
        self.logger.error("get_datasets: mode %s not yet implemented", mode)
    for pol in [0,1]:
      self.avg_diff[pol] /= len(self.indices)
  
  def extract_window(self, xlimits, dsetidx=None, frame="RADI-LSR"):
    """
    Extracts a subset of a avg_diff spectrum
    """
    if self.avg_diff == {} or dsetidx:
      self.beam_switched_average(index=dsetidx, mode='LINE-PBSW')
    # we assume that all scans in all datasets have the same number of channels
    first_idx = self.indices[0]
    first_scan = self.dataset[first_idx].data.keys()[0]
    first_ds = self.dataset[first_idx]
    num_chan = first_ds.data[first_scan].shape[0]
    if frame:
      if frame == "RADI-LSR":
        x = ds.compute_X_axis("RADI-LSR",ds.dss)
      elif frame == "RADI-OBJ":
        if re.match('67P', first_ds.source):
          vspline = self.load_ephemeris('67P')
          x = first_ds.compute_X_axis("RADI-OBJ", first_ds.dss, vspline=vspline)
        else:
          self.logger.error("plot_spectra: no ephemeris for %s", self.object)
      else:
        x = first_ds.compute_X_axis(frame, first_ds.dss)
    else:
      frame = firs_ds.frame
    self.logger.debug("extract_window: x=%s", x)
    self.frame = frame
    # extract the designated range
    x1,x2 = xlimits
    if x[0] < x[-1]:
      ch1 = nearest_index(x, x1)
      ch2 = nearest_index(x, x2)
    else:
      ch1 = nearest_index(x, x2)
      ch2 = nearest_index(x, x1)
    extract = {}
    for pol in self.avg_diff.keys():
      extract[pol] = self.avg_diff[pol][ch1:ch2]
    d = Data(x[ch1:ch2], extract)
    d.frame = frame
    d.channels = (ch1,ch2)
    return d

#################################### Methods ##################################

def parse_filename(fname):
  """
  Parse an SAO .hf5 filename for date and source
  """
  parts = fname[7:-8].split('.')
  UNIXtime = float('.'.join(parts[:2]))
  sourcename = '.'.join(parts[2:])
  return UNIXtime, sourcename

