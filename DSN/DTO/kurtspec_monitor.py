"""
Produce summary plots

Example:
  python kurtspec_monitor.py --date=2018/188 --dss=14 --workstation=gpu1
"""
import cPickle
import glob
import h5py
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import os
import socket
import sys
import time

from collections import OrderedDict
from os.path import basename, dirname, exists, splitext
from time import gmtime, strftime

from Data_Reduction import get_obs_dirs
from Data_Reduction.DSN import path_to_remote
from DatesTimes import UnixTime_to_MPL
from local_dirs import wvsr_dir
from support.logs import get_loglevel, initiate_option_parser, init_logging

logger = logging.getLogger(__name__)
  
class KurtspecSummarizer(object):
  """
  class to examine spectra of 1-sec kurtspec packet averages
  
  Attributes::
    data_array  - contents of pickle file
    doy         - from file name
    dyn_spec    - freq (Y) and time (X) 2D image indexed on type and signal
    filename    - name of picklefile
    logger      - logging.Logger instance
    signal      - signal name dict, keyed on 'I'/'Q' as BDDP (band, dss, pol)
    times       - UNIX time for each average
    timestr     - time as HHMMSS from file name
    workstation - local host
    year        - file name
  """
  array_index = {"power":    {"I": 0,
                              "Q": 2},
                 "kurtosis": {"I": 1,
                              "Q": 3}}

  def __init__(self, workstation=None, filename=None):
    """
    Initializes KurtspecSummarizer to produce summaries of a session.
    
    The workstation which has the averaged data file must be specified.
    
    If no file name is given, the user can later open one with 'open_averages'
    (for a summary file from the observing session) or 'open_spectra' (for data
    previously extracted from a summary file).
    
    @param workstation : name of host which received the data
    @type  workstation : str
    
    @param filename : file to open
    @type  filename : str
    """
    self.logger = logging.getLogger(logger.name +
                                                ".{}".format(self.__class__.__name__))
    self.logger.debug("__init__: logger is %s", self.logger.name)
    
    self.workstation = workstation
    if filename:
      self.filename = filename
      self.session_path = dirname(self.filename)+'/'
      if basename(filename)[:3] == "mon":
        self.open_averages()
      elif basename(filename)[:4] == "spec":
        self.open_spectra()
      else:
        self.logger.error("__init__: file name %s not recognized", filename)
      # find the alternate session path
      self.dss = {}
      for pol in self.signal.keys():
        self.dss[pol] = self.signal[pol][1:3]
    else:
      self.logger.warning("__init__: no file opened; use 'open_averages(name)'")

  def open_spectra(self):
    """
    open an averaged data file
    """
    self.logger.debug("open_spectra: for %s", self.filename)
    dfile = open(self.filename, "rb+")
    corename = splitext(basename(filename))[0]
    self.timestr = corename.split('-')[-1]
    data_dict = cPickle.load(dfile)
    self.signal = data_dict['signals']
    self.times = data_dict['times']
    self.data_array = data_dict['data']
    dfile.close()
    
  def open_averages(self):
    """
    open a file of time- averaged spectra from the monitoring session
    """
    self.logger.debug("open_averages: for %s", self.filename)
    dfile = open(self.filename, "rb+")
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
    num_spectra = self.data_array.shape[2]
    self.get_session_times()
    data_dict = OrderedDict()
    data_dict['signals'] = self.signal
    data_dict['times'] = self.times[:num_spectra]
    data_dict['data'] = self.data_array
    array_file_name = open(self.session_path+"spectra-"+self.timestr+".pkl",
                           "wb")
    cPickle.dump(data_dict, array_file_name)
    array_file_name.close()

  def get_session_times(self):
    """
    set start and end UNIXtime attributes and return matplotlibdates
    
    The time information is in the corresponding HDF5 raw data file.
    """
    # convert .pkl file name to .hdf5 file name
    h5name = basename(self.filename).replace('mon','kurt').replace('pkl','hdf5')
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
    #   truncate so the len of times matches the array time axis
    seconds = rawdata.keys()[:self.data_array.shape[2]] 
    seconds.sort()
    self.logger.debug("get_session_times: found %s", seconds)
    self.times = []
    for second in seconds:
      self.times.append(rawdata[second]['time'].value)
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
    extracts dynamic spectrum type 'power' or 'kurtosis' for input 'I' or 'Q'
    """
    self.dyn_spec = {}
    for moment in ["power","kurtosis"]:
      self.dyn_spec[moment] = {}
      for inpt in ["I","Q"]:
        index = KurtspecSummarizer.array_index[moment][inpt]
        if moment == "power":
          self.dyn_spec[moment][inpt] = numpy.log10(self.data_array[:,index,:])
        else:
          self.dyn_spec[moment][inpt] = self.data_array[:,index,:]
    return self.dyn_spec
    
  def average_power(self, inpt, log=False, first=0, end=1800):
    """
    The data array has shape (1024, 4, num_seconds)
    """
    index = KurtspecSummarizer.array_index["power"][inpt]
    # take the designated array of shape (1024, num_seconds) and average 
    # on the time axis between the specified times
    mean = self.data_array[:,index,first:end].mean(axis=1)
    if log:
      avg_pwr = numpy.log10(mean)
    else:
      avg_pwr = mean
    return avg_pwr

  def kurtosis_timeseries(self, inpt):
    """
    Calculate the frequency averaged kurtosis
    
    @param inpt : input I or Q
    @type  inpt : str
    """
    index = KurtspecSummarizer.array_index["kurtosis"][inpt]
    mean = self.data_array[50:-50,index,:].mean(axis=0)
    return mean

class KurtspecMonitor(KurtspecSummarizer):
  """
  Displays session summaries
  """
  def __init__(self, workstation=None, filename=None, thumbscale=0.1):
    """
    """
    KurtspecSummarizer.__init__(self, workstation=workstation, filename=filename)
    self.logger = logging.getLogger(logger.name +
                                                ".{}".format(self.__class__.__name__))
    self.logger.debug("__init__: logger is %s", self.logger.name)
    self.thumbpath = self.session_path + 'thumbnails/'
    self.thumbscale = thumbscale
    if exists(self.thumbpath):
      pass
    else:
      os.makedirs(self.thumbpath)
  
  def format_file_name(self, name, sigtype, UNIXtime):
    """
    format a filename with signal and time
    
    @param name : like 'passband'
    @type  name : str
    
    @param sigtype : BDDP for band, DSS, pol
    @type  sigtype : str
    
    @param time : UNIX time
    @type  time : float
    """
    return name+"-"+sigtype+"-"+time.strftime("%Y-%j-%H%M",
                                              time.gmtime(UNIXtime))
    
  def plot_spectra(self, RFin, save=True, log=True):
    """
    plot average passband
    
    @param RFin - I or Q
    @type  RFin - str
    
    @param save - save the plot to a file
    @type  save - bool
    
    @param log - logarithmic or linear
    @type  log - bool
    """
    times = self.get_MPL_times()
    for start_idx in range(0, len(times), 1800):
      self.logger.debug("plot_spectra: %s index from %d to %d",
                        self.signal[RFin], start_idx, start_idx+1800)
      last_idx = min(start_idx+1800, len(times)-1)
      avg_pwr = self.average_power(RFin, log=log, first=start_idx, end=last_idx)
      self.logger.debug("plot_spectra: average power: %s", avg_pwr)
      fig = self.plot_segment_spectra(avg_pwr, RFin, first=start_idx, end=last_idx,
                                      save=save, log=log)
      filename = self.format_file_name("passband",
                                       self.signal[RFin],
                                       self.times[start_idx])+".png"
      fig.savefig(self.session_path + filename)
      self.logger.debug("plot_spectra: saved to %s", 
                          self.session_path + filename)
      thumbfile = self.thumbpath + filename
      matplotlib.image.thumbnail(self.session_path+filename, thumbfile,
                                 scale=self.thumbscale)
      plt.close(fig)
    
  def plot_segment_spectra(self, avg_pwr, RFin, first=0, end=1800,
                           save=True, log=True):
    """
    Note that the default for 'log' differs from the one in 'average_power'
    """
    self.logger.debug("plot_segment_spectra: entered")
    fig = plt.figure(figsize=(1.25*12.80, 1.25*10.24), dpi=100)
    ax = fig.add_subplot(111)
    self.logger.debug("plot_segment_spectra: plotting...")
    ax.plot(numpy.arange(0, 650, 650./1024), avg_pwr)
    ax.set_xlim(0, 650)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Frequency (MHz)')
    ax.grid()
    ax.set_title("Average Passband " + self.signal[RFin] 
                 + " from " + time.strftime("%H%M", time.gmtime(first))
                 + " to "   + time.strftime("%H%M", time.gmtime(end))   )
    self.logger.debug("plot_segment_spectra: finished")
    return fig
  
  def plot_segment_kurtosis(self, kurtosis, sigtype,first=0, end=1800,
                            save=True, log=True):
    """
    plot the frequency averaged kurtosis for the designated times
    """
    width = end-first
    axes_width_inches = 18.0*width/1800
    axes_height_inches = 10.24
    figwidth_inches = 1.25*axes_width_inches
    figheight_inches = 1.25*axes_height_inches
    fig = plt.figure(figsize=(figwidth_inches, figheight_inches), dpi=100)
    ax = fig.add_subplot(111)
    times = self.get_MPL_times()
    #ax.plot(times[first:end], kurtosis['I'][first:end], label="I")
    ax.plot(times[first:end], kurtosis[first:end])
    ax.grid(True)
    avkurt = kurtosis.mean()
    krtrng = kurtosis.max()-kurtosis.min()
    ax.set_ylim(avkurt-5*krtrng, avkurt+5*krtrng)
    ax.set_title(   self.signal[sigtype]
                 +" from "+time.strftime("%H%M", time.gmtime(first))
                 +" to "  +time.strftime("%H%M", time.gmtime(end))      )
    ax.xaxis_date()
    fig.autofmt_xdate()
    return fig
    
  def plot_kurtosis(self, RFin, save=True, log=True):
    """
    """
    times = self.get_MPL_times()
    for start_idx in range(0, len(times), 1800):
      kurtosis = self.kurtosis_timeseries(RFin)
      self.logger.debug("plot_kurtosis: %s index from %d to %d",
                        self.signal[RFin], start_idx, start_idx+1800)
      last_idx = min(start_idx+1800, len(times)-1)
      fig = self.plot_segment_kurtosis(kurtosis, RFin,
                                       first=start_idx, end=last_idx,
                                       save=save, log=log)
      filename = self.format_file_name("kurtosis",
                                       self.signal[RFin],
                                       self.times[start_idx])+".png"
      fig.savefig(self.session_path + filename)
      self.logger.debug("plot_kurtosis: saved to %s", 
                          self.session_path + filename)
      thumbfile = self.thumbpath + filename
      matplotlib.image.thumbnail(self.session_path+filename, thumbfile,
                                 scale=self.thumbscale)
      plt.close(fig)
    
  def show_spectrogram(self, sigtype, RFin, first=0, end=1800, save=True):
    """
    show one dynamic spectrum
    
    This assumes that the times are separated by one second
    
    @param sigtype : 'power' or 'kurtosis'
    @type  sigtype : str
    
    @param first : index of first time to include
    @type  first : int
    
    @param last : index of last time (-1) to include
    @type  last : int
    """
    width = end-first
    axes_width_inches = 18.0*width/1800
    axes_height_inches = 10.24
    figwidth_inches = 1.25*axes_width_inches
    figheight_inches = 1.25*axes_height_inches
    times = self.get_MPL_times()
    self.logger.debug("show_spectrogram: %s index from %d to %d",
                      sigtype, first, end)
    self.logger.debug("show_spectrogram: %s time from %f to %f",
                      sigtype, times[first], times[end])
    fig = {}; img = {}
    fig = plt.figure(figsize=(figwidth_inches, figheight_inches), dpi=100)
    ax = plt.subplot(111)
    ax.set_position((0.1,0.2,0.8,0.8))
    img = ax.imshow(self.dyn_spec[sigtype][RFin][:,first:end+1],
                    extent=[times[first], times[end],
                            0,            650        ],
                    aspect='auto', origin='lower')
    ax.grid(True)
    ax.get_xaxis().set_visible(True)
    ax.xaxis_date()
    ax.set_ylabel("Frequency (MHz)")
    ax.set_title(self.signal[RFin])
    cbaxes = fig.add_axes([0.95, 0.2, 0.03, 0.8]) 
    cb = plt.colorbar(img, cax = cbaxes)
    fig.autofmt_xdate()
    return fig, img
    
  def show_spectrograms(self, sigtype, RFin):
    """
    for the y-axis, about a pixel for every two seconds
    
    @param sigtype - 'power' or 'kurtosus'
    @type  sigtype - str
    
    @param RFin - 'I' or 'Q'
    @type  RFin - str
    """
    times = self.get_MPL_times()
    for start_idx in range(0, len(times), 1800):
      self.logger.debug("show_spectrograms: %s index from %d to %d",
                        sigtype, start_idx, start_idx+1800)
      last_idx = min(start_idx+1800, len(times)-1)
      fig, img = self.show_spectrogram(sigtype,
                                       RFin,
                                       first=start_idx,
                                       end=last_idx)
      filename = self.format_file_name("specgram-"+sigtype,
                                      self.signal[RFin],
                                      self.times[start_idx])+".png"
      fig.savefig(self.session_path + filename)
      self.logger.debug("show_spectrograms: saved to %s", 
                        self.session_path + filename)
      thumbfile = self.thumbpath + filename
      matplotlib.image.thumbnail(self.session_path+filename,
                                 thumbfile, scale=self.thumbscale)
      plt.close(fig)
    self.logger.debug("show_spectrograms: %s completed", sigtype)
    
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
                 help = 'DSN station number for signal input I')
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
    mylogger.debug("processing %s", filename)
    corename = splitext(basename(filename))[0]
    mylogger.debug("core file name is %s", corename)
    timestr = corename.split('-')[-1]
    process_file = True
    if corename[:3] == "mon":
      # does the 'spectra' file exist?
      spec_exists = False
      for name in files:
        bname = basename(name)
        mylogger.debug("checking %s for 'spec' with %s", name, timestr)
        if bname[:4] == "spec" and bname.index(timestr) > 0:
          # don't reprocess this file; use the spectra file
          mylogger.info("ignoring %s", bname)
          spec_exists = True
          break
        else:
          mylogger.debug("%s has no corresponding spectra file", bname)
      if spec_exists:
        process_file = False
    if process_file:
      mylogger.debug("using %s", filename)
      mon = KurtspecMonitor(workstation=args.workstation, filename=filename)
      for RFin in "I","Q":
        mon.plot_spectra(RFin)
        mon.plot_kurtosis(RFin)
        mon.extract_dynamic_spectra()
        mon.show_spectrograms("power", RFin)
        mon.show_spectrograms("kurtosis", RFin)
    else:
      mylogger.debug("%s not processed", filename)
  mylogger.debug("done")
