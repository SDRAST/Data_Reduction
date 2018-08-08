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

from os.path import basename, dirname, exists, splitext
from time import gmtime, strftime

from Data_Reduction import get_obs_dirs
from Data_Reduction.DSN import path_to_remote
from DatesTimes import UnixTime_to_MPL
from support.logs import get_loglevel, initiate_option_parser, init_logging

logger = logging.getLogger(__name__)

# these were determined by inspection of the pasbands
#   at present the K data (usually in the Q channels) are the same as what is
# in the I channel
valid_freqs = {"S": (130,160), "X": (30,580), "K": (30,580)}    # MHz
  
class KurtspecSummarizer(object):
  """
  class to examine spectra of 1-sec kurtspec packet averages
  
  Attributes::
    data        - contents of HDF5 monitor file
    doy         - from file name
    dyn_spec    - freq (Y) and time (X) 2D image indexed on type and signal
    filename    - name of picklefile
    logger      - logging.Logger instance
    signal      - signal name as BDDP (band, dss, pol)
    times       - UNIX time for each average
    timestr     - time as HHMMSS from file name
    workstation - local host
    year        - file name
  """
  def __init__(self, workstation=None, filename=None, session_path=None):
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
      self.session_path = session_path
      self.open_averages()
    else:
      self.logger.warning("__init__: no file opened; use 'open_averages(name)'")
    
  def open_averages(self):
    """
    open a file of time- averaged spectra from the monitoring session
    """
    self.logger.debug("open_averages: for %s", self.filename)
    data = h5py.File(self.filename)
    self.signal = data.attrs['signal']
    self.band = self.signal[0]
    self.times = data['time'].value
    self.power = data['power'].value
    self.kurtosis = data['kurtosis'].value
  
  def get_MPL_times(self):
    """
    """
    return UnixTime_to_MPL(self.times)
        
  def average_power(self, log=False, first=0, end=1800):
    """
    Average power over the time axis to get a passband
    
    The power array has shape (num_seconds, 1024)
    """
    if log:
      avg_pwr = numpy.log10(self.power[first:end,:].mean(axis=0))
    else:
      avg_pwr = self.power[first:end,:].mean(axis=0)
    return avg_pwr
    
  def extract_dynamic_spectra(self):
    """
    extracts dynamic spectrum type 'power' or 'kurtosis'
    
    It appears that the last kurtosis spectrum is zeros
    """
    self.dyn_spec = {}
    for moment in ["power","kurtosis"]:
      self.dyn_spec[moment] = {}
      if moment == "power":
          self.dyn_spec[moment] = numpy.log10(self.power).T
      else:
          self.dyn_spec[moment] = self.kurtosis.T
    return self.dyn_spec


class KurtspecMonitor(KurtspecSummarizer):
  """
  Displays session summaries
  
  Attributes::
    logger     - logging.Logger object
    thumbpath  - path to thumbnails
    thunmscale - thumb image scale relative to original
  
  Attributes inherited from KurtspecSummarizer::
    data        - contents of HDF5 monitor file
    doy         - from file name
    dyn_spec    - freq (Y) and time (X) 2D image indexed on type and signal
    filename    - name of picklefile
    logger      - logging.Logger instance
    signal      - signal name as BDDP (band, dss, pol)
    times       - UNIX time for each average
    timestr     - time as HHMMSS from file name
    workstation - local host
    year        - file name
  """
  def __init__(self, workstation=None, filename=None, session_path=None, 
               thumbscale=0.1):
    """
    """
    KurtspecSummarizer.__init__(self, workstation=workstation, filename=filename,
                                session_path=session_path)
    self.logger = logging.getLogger(logger.name +
                                                ".{}".format(self.__class__.__name__))
    self.logger.debug("__init__: logger is %s", self.logger.name)
    self.thumbpath = self.session_path + 'thumbnails/'
    self.thumbscale = thumbscale
    self.logger.debug("__init__: thumbs in %s", self.thumbpath)
    if exists(self.thumbpath):
      pass
    else:
      os.makedirs(self.thumbpath, mode=0775)
  
  def format_file_name(self, name, UNIXtime):
    """
    format a filename with signal and time
    
    @param name : like 'passband'
    @type  name : str
   
    @param time : UNIX time
    @type  time : float
    """
    return name+"-"+self.signal+"-"+time.strftime("%Y-%j-%H%M",
                                              time.gmtime(UNIXtime))
    
  def plot_spectra(self, save=True, log=True):
    """
    plot average passband
    
    @param save - save the plot to a file
    @type  save - bool
    
    @param log - logarithmic or linear
    @type  log - bool
    """
    times = self.get_MPL_times()
    for start_idx in range(0, len(times), 1800):
      self.logger.debug("plot_spectra: %s index from %d to %d",
                        self.signal, start_idx, start_idx+1800)
      last_idx = min(start_idx+1800, len(self.times)-1)
      avg_pwr = self.average_power(log=log, first=start_idx, end=last_idx)
      self.logger.debug("plot_spectra: average power: %s", avg_pwr)
      fig = self.plot_segment_spectra(avg_pwr,
                                      first=start_idx, end=last_idx,
                                      save=save, log=log)
      filename = self.format_file_name("passband",
                                       self.times[start_idx])+".png"
      fig.savefig(self.session_path + filename)
      self.logger.debug("plot_spectra: saved to %s", 
                          self.session_path + filename)
      thumbfile = self.thumbpath + filename
      matplotlib.image.thumbnail(self.session_path+filename, thumbfile,
                                 scale=self.thumbscale)
      plt.close(fig)
    
  def plot_segment_spectra(self, avg_pwr, first=0, end=1800,
                           save=True, log=True):
    """
    Note that the default for 'log' differs from the one in 'average_power'
    """
    fig = plt.figure(figsize=(1.25*12.80, 1.25*10.24), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(numpy.arange(0, 650, 650./1024), avg_pwr)
    ax.set_xlim(0, 650)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Frequency (MHz)')
    ax.grid()
    ax.set_title("Average Passband " + self.signal 
                 + " from " + time.strftime("%H%M", time.gmtime(first))
                 + " to "   + time.strftime("%H%M", time.gmtime(end))   )
    return fig
  
  def plot_segment_kurtosis(self, first=0, end=1800, save=True):
    """
    plot the frequency averaged kurtosis for the designated times
    
    Frequency range for averaging is in global 'valid_freqs' determined by
    inspection of the passbands.
    """
    width = end-first
    axes_width_inches = 18.0*width/1800
    axes_height_inches = 10.24
    figwidth_inches = 1.25*axes_width_inches
    figheight_inches = 1.25*axes_height_inches
    fig = plt.figure(figsize=(figwidth_inches, figheight_inches), dpi=100)
    ax = fig.add_subplot(111)
    times = self.get_MPL_times()
    # take the mean over the valid frequency range
    minfreq = valid_freqs[self.band][0]
    maxfreq = valid_freqs[self.band][1]
    avg_kurtosis = self.kurtosis[first:end, minfreq:maxfreq].mean(axis=1)
    ax.plot(times[first:end], avg_kurtosis)
    ax.grid(True)
    ax.set_xlim(times[first], times[end])
    avkurt = avg_kurtosis.mean()
    krtrng = avg_kurtosis.max()-avg_kurtosis.min()
    ax.set_ylim(avkurt-2*krtrng, avkurt+2*krtrng)
    ax.set_title(   self.signal
                 +" from "+time.strftime("%H%M", time.gmtime(first))
                 +" to "  +time.strftime("%H%M", time.gmtime(end))      )
    ax.xaxis_date()
    fig.autofmt_xdate()
    return fig
    
  def plot_kurtosis(self, save=True):
    """
    display a spectrogram of the original data
    """
    times = self.get_MPL_times()
    for start_idx in range(0, len(times), 1800):
      self.logger.debug("plot_kurtosis: %s index from %d to %d",
                        self.signal, start_idx, start_idx+1800)
      last_idx = min(start_idx+1800, len(times)-1)
      fig = self.plot_segment_kurtosis(first=start_idx, end=last_idx,
                                       save=save)
      filename = self.format_file_name("kurtosis",
                                       self.times[start_idx])+".png"
      fig.savefig(self.session_path + filename)
      self.logger.debug("plot_kurtosis: saved to %s", 
                          self.session_path + filename)
      thumbfile = self.thumbpath + filename
      matplotlib.image.thumbnail(self.session_path+filename, thumbfile,
                                 scale=self.thumbscale)
      plt.close(fig)
    
  def show_spectrogram(self, sigtype, first=0, end=1800, save=True):
    """
    show one dynamic spectrum
    
    This assumes that the times are separated by one second
    
    For a typical image on an 8x6 in frame, the portion of the figure allocated
    to the image is 0.8 starting at 0.1. So that means the left side of the
    image is at 0.8 in and the right is at 7.2 in. The same proportion is used
    for high, so the image starts 0.6 in from the bottom and ends at 5.4 in.
    The calculations below should give roughly one pixel per second and per MHz
      I suspect that the figure dimensions are calculated so that the left and
    right margins are 0.1 of the width and height respectively
    
    Regarding aspect, we want one pixel per MHz and second.  The X-axis is in
    matplotlib datenums, which is days. So one pixel is 1.16e-5 datenums.
    'aspect' is the number of width units divided by the number of height
    units
    
    @param sigtype : 'power' or 'kurtosis'
    @type  sigtype : str
    
    @param first : index of first time to include
    @type  first : int
    
    @param end : index of last time (-1) to include
    @type  end : int
    """
    width = end-first
    axes_width_inches = 18.0*width/1800 # pixel/sec @ 100 dpi
    axes_height_inches = 10.24          # pixel/MHz @ 100 dpi
    self.logger.debug("show_spectrogram: image is %.1f x %.1f",
                      axes_width_inches, axes_height_inches)
    figwidth_inches =  0.8 + axes_width_inches/0.8
    figheight_inches = 0.6 + axes_height_inches/0.8
    self.logger.debug("show_spectrogram: figure is %.1f x %.1f",
                      figwidth_inches, figheight_inches)
    times = self.get_MPL_times()
    self.logger.debug("show_spectrogram: %s index from %d to %d",
                      sigtype, first, end)
    self.logger.debug("show_spectrogram: %s time from %f to %f",
                      sigtype, times[first], times[end])
    fig = {}; img = {}
    fig = plt.figure(figsize=(figwidth_inches, figheight_inches), dpi=100)
    ax = plt.subplot(111)
    ax.set_position((0.1,0.2,0.8,0.8))
    data = self.dyn_spec[sigtype][:,first:end]
    self.logger.debug("show_spectrogram: data shape is %s", data.shape)
    if sigtype == "power":
      img = ax.imshow(data,
                      extent=[times[first], times[end],
                              0,            650        ],
                      aspect='auto', origin='lower')
    elif sigtype == "kurtosis":
      minfreq = valid_freqs[self.band][0]
      maxfreq = valid_freqs[self.band][1]
      min_kurt = self.dyn_spec['kurtosis'][minfreq:maxfreq,:-1].min()
      max_kurt = self.dyn_spec['kurtosis'][minfreq:maxfreq,:-1].max()
      meankurt = self.dyn_spec['kurtosis'][minfreq:maxfreq,:-1].mean()
      img = ax.imshow(data,
                      extent=[times[first], times[end],
                              0,            650        ],
                      vmin=meankurt-1, vmax=meankurt+1,
                      aspect='auto', origin='lower')
    else:
      raise RuntimeError("invalid signal type:", sigtype)
    ax.grid(True)
    ax.get_xaxis().set_visible(True)
    ax.xaxis_date()
    ax.set_ylabel("Frequency (MHz)")
    ax.set_title(self.signal)
    cbaxes = fig.add_axes([0.95, 0.2, 0.03, 0.8]) 
    cb = plt.colorbar(img, cax = cbaxes)
    fig.autofmt_xdate()
    return fig, img
    
  def show_spectrograms(self, sigtype):
    """
    for the y-axis, about a pixel for every two seconds
    
    @param sigtype - 'power' or 'kurtosus'
    @type  sigtype - str
    """
    times = self.get_MPL_times()
    for start_idx in range(0, len(times), 1800):
      self.logger.debug("show_spectrograms: %s index from %d to %d",
                        sigtype, start_idx, start_idx+1800)
      last_idx = min(start_idx+1800, len(times)-1)
      fig, img = self.show_spectrogram(sigtype,
                                       first=start_idx,
                                       end=last_idx)
      filename = self.format_file_name("specgram-"+sigtype,
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
    
  projdatapath, projworkpath, datapath = \
                    get_obs_dirs("PESD", args.dss, year, doy)
  mylogger.debug("project data path: %s", projdatapath)
  mylogger.debug("project work path: %s", projworkpath)
  mylogger.debug("raw data path: %s", datapath)
  
  if args.workstation:
    datapath = path_to_remote(args.workstation, projdatapath)
    sessionpath = path_to_remote(args.workstation, projworkpath)
  else:
    datapath = projdatapath
    sessionpath = projworkpath
  mylogger.debug("path for data: %s", datapath)
  files = glob.glob(datapath+"*.hdf5")
  mylogger.debug("found: %s", files)
  
  for filename in files:
    mylogger.debug("processing %s", filename)
    corename = splitext(basename(filename))[0]
    mylogger.debug("core file name is %s", corename)
    timestr = corename.split('-')[-1]
    mylogger.debug("using %s", filename)
    mon = KurtspecMonitor(workstation=args.workstation, filename=filename, 
                            session_path=sessionpath)
    mon.plot_spectra()
    mon.plot_kurtosis()
    mon.extract_dynamic_spectra()
    mon.show_spectrograms("power")
    mon.show_spectrograms("kurtosis")
  mylogger.debug("done")
