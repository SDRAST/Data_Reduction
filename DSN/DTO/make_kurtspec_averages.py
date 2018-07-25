"""
Convert raw kurtosis spectrometer data into 1 sec averages

A kurtosis spectrometer raw data file consists a number of groups containing
one second of data.  The datasets in a group are::
 hdr         - first eight bytes of a packet
 krt-I       - kurtosis spectrum of the I input
 krt-Q       - kurtosis spectrum of the Q input
 pkt cnt sec - 
 pwr-I       - power spectrum of the I input
 pwr-Q       - power spectrum of the Q input
 raw pkt cnt - number of the packet
 sec cnt     -
 time        - time when the packet was received
Each dataset has 10,000 packets (5,000 until a firmware issue is resolved).

The output file has the signal as an attribute. It has datasets 'actime',
'pkt_cnt', 'power' and 'kurtosis'. The first index in the array is a counter
from 0 to number of averages - 1.  The second index, if any, is spectrum
channel.
"""
import cPickle
import glob
import h5py
import logging
import numpy
import os
import sys
import time

from Data_Reduction import get_obs_dirs
from Data_Reduction.DSN import path_to_remote
from support.logs import get_loglevel, initiate_option_parser, init_logging

logger = logging.getLogger(__name__)

class DataAverager(object):
  """
  """
  def __init__(self, filename=None):
    """
    """
    self.logger = logging.getLogger(logger.name +
                                         ".{}".format(self.__class__.__name__))
    self.logger.debug("__init__: logger is %s", self.logger.name)
    if filename:
      self.open_data_file(filename)
    else:
      self.logger.warning("__init__: use 'open_data_file' to load data")
  
  def open_data_file(self, filename):
    """
    """
    self.logger.debug("open_data_file: %s", filename)
    self.data = h5py.File(filename)
    self.grp_keys = self.data.keys()
    # extract signals and stations from header
    header = self.data[self.grp_keys[0]]['hdr'][0]
    self.signal = {'I': "".join(header[:4]),
                   'Q': "".join(header[4:8])}
    self.logger.info("open_data_file: signals are %s", self.signal)
    dss = {"I": int(self.signal['I'][1:3]), "Q": int(self.signal['Q'][1:3])}
    band = {"I": self.signal['I'][0], "Q": self.signal['Q'][0]}
    self.logger.info("open_data_file: stations are %s", dss)
    # make destination paths and file names
    corename = os.path.splitext(os.path.basename(filename))[0]
    ftype,yearstr,doystr,timestr = corename.split('-')
    year = int(yearstr); doy = int(doystr)
    hr = int(timestr[:2]); mn = int(timestr[2:4]); sec = int(timestr[4:])
    self.logger.debug("open_data_file: session is %4d/3d %02d:%02d:%02d",
                      year, doy, hr, mn, sec)
    datapath = {}
    fname = {}
    #sessionpath = {}
    for inp in ["I", "Q"]:
      #   create the monitor data files
      fname[inp] = "mon-" + band[inp] + time.strftime("-%Y-%j-%H%M%S.hdf5",
                                          (year, 0, 0, hr, mn, sec, 0, doy, 0))
      session_data_path, projworkpath, ignore2 = \
                                      get_obs_dirs("PESD", dss[inp], year, doy)
      self.logger.debug("open_data_file: mon data for %s at %s",
                        inp, session_data_path)
      datapath[inp] = path_to_remote(args.workstation, session_data_path)
      #sessionpath[inp] = path_to_remote(args.workstation, projworkpath)
    self.create_monitor_datasets(datapath, fname)

  def create_monitor_datasets(self, session_data_path, fname):
    """
    opens files for 1-sec monitor data
    
    There is a file for the I input and one for the Q input
    
    @param session_data_path : destination of file in/usr/local/project_data
    @type  session_data_path : str
    
    @param fname : name of output file
    @type  fname : dict of str
    """
    self.monitor_file = {} # path and file name
    self.aqtime = {}       # UNIX time data was received
    self.pkt_cnt = {}      #
    self.power = {}        # power spectra
    self.kurtosis = {}     # kurtosis spectra
    for RF in "I", "Q":
      # replace if it exists?
      fullpath = session_data_path[RF] + fname[RF]
      if os.path.exists(fullpath):
        if raw_input("create_monitor_datasets: replace %s for %s?"
                     % (fname[RF], RF))[0] in "yY":
          os.remove(fullpath)
        else:
          raise RuntimeError("%s already exists", os.basename(fullpath))
      self.logger.debug("create_monitor_datasets: creating %s for %s",
                        fullpath, RF)
      try:
        self.monitor_file[RF] = h5py.File(session_data_path[RF] + fname[RF])
      except IOError:
        os.makedirs(session_data_path[RF])
        self.monitor_file[RF] = h5py.File(session_data_path[RF] + fname[RF])
      # create the datasets, initial length of 1
      self.monitor_file[RF].attrs['signal'] = self.signal[RF]
      self.monitor_file[RF].attrs['channel'] = RF
      self.aqtime[RF] = self.monitor_file[RF].create_dataset('time', (1,),
                                                        maxshape=(None,),
                                                        dtype=numpy.float64)
      self.pkt_cnt[RF] = self.monitor_file[RF].create_dataset('pkt cnt', (1,),
                                                         maxshape=(None,),
                                                         dtype=numpy.int32)
      self.power[RF] = self.monitor_file[RF].create_dataset('power', (1,1024),
                                                        maxshape=(None,1024))
      self.kurtosis[RF] = self.monitor_file[RF].create_dataset('kurtosis',
                                                              (1,1024),
                                                          maxshape=(None,1024))

  def average_one_second(self):
    """
    writes 1-sec averages to file
  
    This takes the four spectra power-I, kurtosis-I, power-Q and kurtosis-Q and
    arranged them in a 2D array which is written to disk
  
    @param inqueue : where the data come from
    """  
    for grp in self.grp_keys:
      if 'one_second' in grp:
        # to handle data pre 2018.200
        counter = int(grp.replace('one_second', ''))
      else:
        counter = int(grp)
      average = {}
      self.logger.debug("average_one_second: first packet counter: %d", 
                        self.data[grp]['raw pkt cnt'][0])
      for IF in ['I', 'Q']:
        self.logger.debug("average_one_second: doing record %d for %s",
                          counter, IF)
        self.aqtime[IF][counter] = self.data[grp]['time'].value       # float64
        self.pkt_cnt[IF][counter] = self.data[grp]['raw pkt cnt'][0]  # int32
        # convert list of 1D arrays into 2D array and average along the list
        # axis 0 is frequency; axis 1 is time
        self.power[IF][counter] = numpy.array(
                                           self.data[grp]['pwr-'+IF]).mean(axis=0) 
        self.kurtosis[IF][counter] = numpy.array(
                                           self.data[grp]['krt-'+IF]).mean(axis=0)
        self.monitor_file[IF].flush()
        self.logger.debug("average_one_second: added row %d to %s for %s",
                     counter, self.monitor_file[IF].file, IF)
      for IF in ['I', 'Q']:
        # size is one more than the current count, which is two more than the
        # current index
        self.aqtime[IF].resize(counter+2, axis=0)
        self.pkt_cnt[IF].resize(counter+2, axis=0)
        self.power[IF].resize(counter+2, axis=0)
        self.kurtosis[IF].resize(counter+2, axis=0)
    
    
if __name__ == "__main__":
  logging.basicConfig(level=logging.DEBUG)
  mylogger=logging.getLogger()

  p = initiate_option_parser("Kurtosis Data Averager","")
  p.usage = "python make_kurtspec_averages.py <kwargs>"
  # Add other options here
  p.add_argument('--date',
                 dest = 'date',
                 type = str,
                 default = None,
                 help = 'Date of observation as YEAR/DOY string')
  p.add_argument('-p', '--project',
                 dest = 'project',
                 type = str,
                 default = "PESD",
                 help = "Project code")
  p.add_argument('-w', '--workstation',
                 dest = 'workstation',
                 type = str,
                 default = None,
                 help = "Workstation for post-processing")
  
  args = p.parse_args()
  
  # This cannot be delegated to another module or class
  mylogger = init_logging(logging.getLogger(),
                          loglevel   = get_loglevel(args.file_loglevel),
                          consolevel = get_loglevel(args.console_loglevel),
                          logname    = args.logpath+__name__+".log")
  mylogger.debug("arguments: %s", args)
  mylogger.debug(" Handlers: %s", mylogger.handlers)

  if args.date:
    yearstr, doystr = args.date.split('/')
    year = int(yearstr)
    doy = int(doystr)
  else:
    mylogger.error(" 'date' is a required argument")
    sys.exit(1)
    
  datapath = "/data/kurtspec/%4d/%03d/" % (year, doy)
  mylogger.debug("data at %s", datapath)
    
  files = glob.glob(datapath+"*.hdf5")
  mylogger.debug("found: %s", files)
  
  averager = DataAverager()
  for filename in files:
    mylogger.debug("processing %s", filename)
    averager.open_data_file(filename)
    averager.average_one_second()
    averager.data.close()
    for inp in ['I', 'Q']:
      averager.monitor_file[inp].close()
  mylogger.debug("done")
