"""
"""
import cPickle
import glob
import h5py
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import os.path
import socket

from Data_Reduction import get_obs_dirs
from DatesTimes import UnixTime_to_MPL
from local_dirs import wvsr_dir
from support.logs import get_loglevel, initiate_option_parser, init_logging

logger = logging.getLogger(__name__)

def path_to_remote(workstation, local_path):
  """
  """
  thishost = socket.gethostname()
  if thishost == workstation:
    datapath = local_path
  else:
    workstationpath = "/home/kuiper/mnt/"+workstation
    if thishost == 'dto':
      datapath = workstationpath+local_path
    else:
      dtopath = "/home/kuiper/mnt/dto"
      if thishost == 'kuiper':
        datapath = dtopath +workstationpath+local_path
      else:
        return None
  return datapath
  
class KurtspecMonitor(object):
  """
  class to examine spectra of 1-sec kurtspec packet averages
  """
  array_index = {"power":    {"I": 0,
                              "Q": 2},
                 "kurtosis": {"I": 1,
                              "Q": 3}}

  def __init__(self, workstation=None, filename=None):
    """
    """
    self.logger = logging.getLogger(logger.name +
                                                ".{}".format(self.__class__.__name__))
    self.logger.debug("__init__: logger is %s", self.logger.name)
    self.workstation = workstation
    if filename:
      data_array = self.open_averages(filename)
    else:
      self.logger.warning("__init__: no file opened; use 'open_averages(name)'")
  
  def open_averages(self, filename):
    """
    """
    dfile = open(filename, "rb+")
    self.filename = filename
    data = cPickle.load(dfile)
    self.data_array = data.reshape(data.shape + (1,))
    count = 0
    reading = True
    while reading:
      try:
        data = cPickle.load(dfile)
        self.data_array = numpy.append(self.data_array, 
                                       data.reshape(data.shape + (1,)),
                                       axis=2)
        count += 1
        print "\r"+str(count),
      except EOFError:
        reading = False
    dfile.close()
    self.get_session_times()
    return self.data_array

  def get_session_times(self):
    """
    set start and end UNIXtime attributes and return matplotlibdates
    
    The time information is in the corresponding HDF5 raw data file.
    """
    # convert .pkl file name to .hdf5 file name
    h5name = os.path.basename(self.filename).replace('mon','kurt').replace('pkl','hdf5')
    self.logger.debug("get_session_times: HDF5 name is %s", h5name)
    # get the date/time information
    dummy, yearstr, doystr, rest = h5name.split('-')
    self.year = int(yearstr)
    self.doy = int(doystr)
    self.timestr = rest.split('.')[0]
    # open the .hdf5 file
    h5file = path_to_remote(self.workstation,
                            "/data/kurtspec/"+yearstr+"/"+doystr+"/"+h5name)
    self.logger.debug("get_session_times: HDF5 file is %s", h5file)
    rawdata = h5py.File(h5file)
    # get the times of each one second average
    seconds = rawdata.keys()
    seconds.sort()
    self.logger.debug("get_session_times: found %s", seconds)
    self.times = []
    for second in seconds:
      self.times.append(rawdata[second]['time'].value)
    self.logger.debug("get_session_times: %s", self.times)
    # get the signal sources
    self.signal = {}
    header = rawdata[seconds[0]]['hdr'][0] # 5000 headers, one per packet
    self.signal['I'] = "".join(list(header[:4]))
    self.signal['Q'] = "".join(list(header[4:]))
    rawdata.close()
  
  def get_MPL_times(self):
    """
    """
    return UnixTime_to_MPL(self.times)
    
  def extract_dynamic_spectra(self):
    """
    """
    self.dyn_spec = {}
    for moment in ["power","kurtosis"]:
      self.dyn_spec[moment] = {}
      for inpt in ["I","Q"]:
        index = KurtspecMonitor.array_index[moment][inpt]
        if moment == "power":
          self.dyn_spec[moment][inpt] = numpy.log10(self.data_array[:,index,:])
        else:
          self.dyn_spec[moment][inpt] = self.data_array[:,index,:]
    return self.dyn_spec
    
  def average_power(self, log=False):
    """
    """
    avg_pwr = {}
    for inpt in ["I","Q"]:
      index = KurtspecMonitor.array_index["power"][inpt]
      mean = self.data_array[:,index,:].mean(axis=1)
      if log:
        avg_pwr[inpt] = numpy.log10(mean)
      else:
        avg_pwr[inpt] = mean
    return avg_pwr

  def kurtosis_timeseries(self):
    """
    """
    kurtosis = {}
    for inpt in ["I","Q"]:
      index = KurtspecMonitor.array_index["kurtosis"][inpt]
      mean = self.data_array[50:-50,index,:].mean(axis=0)
      kurtosis[inpt] = mean
    return kurtosis
    
  def plot_spectra(self, log=True):
    """
    Note that this default differs from the one in 'average_power'
    """
    avg_pwr = self.average_power(log=log)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(numpy.arange(0,650,650./1024), avg_pwr['I'], label="I")
    ax.plot(numpy.arange(0,650,650./1024), avg_pwr['Q'], label="Q")
    ax.set_xlim(0, 650)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Channel')
    ax.legend()
    ax.grid()
    ax.set_title("Average Passband")
    fig.savefig(os.path.dirname(self.filename)
                + "/passband" + "-" + self.timestr + ".png")
  
  def plot_kurtosis(self):
    """
    """
    kurtosis = self.kurtosis_timeseries()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    times = self.get_MPL_times()
    self.logger.debug("plot_kurtosis: %d times", len(times))
    self.logger.debug("plot_kurtosis: times: %s", times)
    ax.plot(times, kurtosis['I'], label="I")
    ax.plot(times, kurtosis['Q'], label="Q")
    ax.grid(True)
    ax.set_ylim(2.8,3.0)
    ax.xaxis_date()
    fig.autofmt_xdate()
    fig.savefig(os.path.dirname(self.filename)
                +"/kurtosis" + "-" + self.timestr + ".png")

  def show_spectrogram(self, type='power'):
    """
    for the y-axis, about a pixel for every two seconds
    """
    times = self.get_MPL_times()
    self.extract_dynamic_spectra()
    height, width = self.dyn_spec['kurtosis']['Q'].shape
    figwidth = 12.*width/1024
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(figwidth,9))
    
    im1 = ax[0].imshow(self.dyn_spec[type]['I'],
                       extent=[times[0], times[-1],0,650],
                       aspect='auto', origin='lower')
    ax[0].grid(True)
    ax[0].get_xaxis().set_visible(False)
    ax[0].xaxis_date()
    ax[0].set_ylabel("Frequency (MHz)")
    ax[0].set_title(type+" I")
    plt.colorbar(im1, ax=ax[0])
    
    im2 = ax[1].imshow(self.dyn_spec[type]['Q'],
                       extent=[times[0],times[-1],0,650],
                       aspect='auto', origin='lower')
    ax[1].grid(True)
    ax[1].xaxis_date()
    ax[1].set_title(type+" Q")
    plt.colorbar(im2, ax=ax[1])
    fig.autofmt_xdate()
    fig.tight_layout()
    plt.savefig(os.path.dirname(self.filename)
                + "/"+type+"-"+"specgram" + "-" + self.timestr + ".png")
    
if __name__ == "__main__":
  logging.basicConfig(level=logging.DEBUG)
  mylogger=logging.getLogger()

  p = initiate_option_parser("Kurtosis monitor","")
  p.usage = "python kurtspec_monitor.py <kwargs>"
  # Add other options here
  p.add_argument('--date',
                 dest = 'date',
                 type = str,
                 default = None,
                 help = 'Date of observation as YEAR/DOY string')
  p.add_argument('-D', '--dss',
                 dest = 'dss',
                 type = int,
                 default = None,
                 help = 'DSN station number')
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
  if args.dss:
    pass
  else:
    mylogger.error(" 'dss' is a required argument")
    sys.exit(1)
  if args.workstation:
    pass
  else:
    mylogger.error(" 'workstation' is a required argument")
    sys.exit(1)
    
  projectdatapath, projworkpath, datapath = \
                    get_obs_dirs("PESD", args.dss, year, doy)
  mylogger.debug("project data path: %s", projectdatapath)
  mylogger.debug("project work path: %s", projworkpath)
  mylogger.debug("raw data path: %s", datapath)
  
  datapath = path_to_remote(args.workstation, projworkpath)
  mylogger.debug("path for data: %s", datapath)
  files = glob.glob(datapath+"*.pkl")
  mylogger.debug("found: %s", files)
  
  for filename in files:
    mon = KurtspecMonitor(workstation=args.workstation, filename=filename)
    mon.plot_spectra()
    mon.plot_kurtosis()
    mon.show_spectrogram(type="power")
    mon.show_spectrogram(type="kurtosis")
    
