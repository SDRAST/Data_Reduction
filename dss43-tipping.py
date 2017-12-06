"""
class for analyzing SAO/TAMS tipping curve data

This can be improved by using the support.tunneling module to check for a
tunnel and creating one if need.
"""
import calendar
import glob
import logging
import math
import numpy
import os
import time

from collections import OrderedDict
from pylab import *
from scipy.optimize import curve_fit

from Data_Reduction import get_obs_session
from Data_Reduction.tipping import airmass
from DatesTimes import calendar_date
from local_dirs import projects_dir
from support import nearest_index
from support.logs import initiate_option_parser, init_logging, get_loglevel
#from support.process import search_response

logger = logging.getLogger(__name__)

saologdir  = "/home/ops/roach_data/sao_test_data/log_dir/" # on crux

IFcolors = ['b', 'g', 'r', 'c']

class TidTipAnalyzer(object):
  """
  Attributes::
    am
    directory
    files
    logger
    tsys
  """
  def __init__(self, remote=True, date=None):
    """
    """
    self.logger = logging.getLogger(logger.name+".TidTipAnalyzer")
    if remote:
      self.directory = os.environ['HOME'] + "/mnt/crux" + saologdir
    else:
      self.directory = saologdir
    self.files = self.sort_by_time(subdir="2014-2016")
    self.files.update(self.sort_by_time())
    if date:
      self.get_session_files(date)

  def sort_by_time(self, name_template="tipping_data", subdir=None):
    """
    Returns a time-ordered dict of tipping files with UNIX times as keys
  
    @param directory : where the files are
    @type  directory : str
  
    @param name_template : part of name common to all files
    @type  name_template : str
  
    @return: dict keyed on filenames
    """
    if subdir:
      if subdir[-1] != "/":
        subdir += "/"
      path = self.directory + subdir
    else:
      path = self.directory
    files = glob.glob(path + name_template +"*.npy")
    filedict = {}
    for name in files:
      unixtime = float(os.path.basename(name)[12:-4])
      filedict[unixtime] = name
    keys = filedict.keys()
    keys.sort()
    ordered = OrderedDict()
    for key in keys:
      ordered[key] = filedict[key]
    return ordered

  def get_session_files(self, date):
    """
    Get the files for the specified data
    """
    self.date = date
    yearstr, doystr = date.split('/')
    self.year = int(yearstr)
    self.DOY = int(doystr)
    year,month,day = calendar_date(self.year, self.DOY)
    UNIXstart = calendar.timegm((self.year,month,day,  0, 0, 0, 0,0))
    UNIXend   = calendar.timegm((self.year,month,day, 23,59,59, 0,0))
    filekeys = self.files.keys()
    first_idx = nearest_index(filekeys, UNIXstart)
    last_idx  = nearest_index(filekeys, UNIXend)
    self.logger.debug("get_session_files: index from %d to %d", first_idx, last_idx)
    if filekeys[first_idx] < UNIXstart:
      first_idx += 1
    if filekeys[last_idx] < UNIXend:
      last_idx += 1
    files = []
    for index in range(first_idx,last_idx):
      files.append(self.files[filekeys[index]])
    self.datefiles = files
    return files

  def get_tipping_data(self, filename):
    """
    Gets tipping data from a file made with observatoryCtrl
    """
    tipfile = open(filename, 'r')
    data = numpy.load(tipfile)
    tipfile.close()
    elevations = data[0][:].astype(float)
    tsys = []
    for item in data[1].tolist():
      tsys.append(item.ravel())
    Tsys = numpy.array(tsys)
    return elevations, Tsys

  def plot_elevations(self, index, tipfiles=None):
    """
    Plots tipping curve data from a specified file
  
    @param tipfiles : dict keyed on filename with file times
    @type  tipfiles : dict
  
    @param index : file in tipfiles to be plotted
    @type  index : int
    """
    if tipfiles:
      filename = tipfiles[index]
    else:
      filename =  self.datefiles[index]
    e,t = self.get_tipping_data(filename)
    figure()
    for chan in range(4):
      plot(e,t[:,chan], IFcolors[chan]+'-')
      plot(e,t[:,chan], IFcolors[chan]+'.', label="IF"+str(chan+1))
    xlabel("Elevation")
    ylabel("System Temperature")
    grid()
    legend(numpoints=1)
    unixtime = float(os.path.basename(filename)[12:-4])
    title(time.ctime(unixtime))
    show()

  def fit_tipping_data(self, index, tipfiles=None, Tatm=250, project=None):
    """
    Fits tipping curve data from a specified file
  
    @param tipfiles : dict keyed on filename with file times
    @type  tipfiles : dict
  
    @param index : file in tipfiles to be plotted
    @type  index : int
  
    @param Tatm : average temperature of the atmosphere
    @type  Tatm : float
  
    @return: float (receiver temperature), float (atmospheric opacity)
    """
    def plot_data(am, tsys):
      """
      """
      fig = figure()
      for IF in range(4):
        plot(am, tsys[:,IF], IFcolors[IF]+'-')
        plot(am, tsys[:,IF], IFcolors[IF]+'.', label="IF%d" % IF)
      legend(loc="upper left", numpoints=1)
      grid()
      title(self.date)
      xlabel("airmass")
      show()
      return fig, gca()
              
    def opacity_fitting(e, Tsysvac, tau):
      """
      radiative transfer equation for atmosphere
      
      Tatm is the average atmospheric temperature along the line of sight
      
      @param e : elevations (deg)
      @type  e : numpy.array of float
      
      @param Tsysvac : system tenperature in a vacuum
      @type  Tsysvac : float
      
      @param tau : opacity
      @type  tau : float
      """
      am = airmass(e)
      return Tsysvac + Tatm*(1 - numpy.exp(-tau*am))
    
    # get the data
    if tipfiles:
      filename = tipfiles[index]
    else:
      filename = self.datefiles[index]
    e,t = self.get_tipping_data(filename)
    self.logger.debug("fit_tipping_data: %d elevations", len(e))
    self.logger.debug("fit_tipping_data: %d Tsys groups", len(t))
    am = airmass(e)
    fig, ax = plot_data(am, t)
    handles, labels = ax.get_legend_handles_labels()
    
    # get the airmass range
    min_am_idx = am.argmin()
    max_am_idx = am.argmax()
    self.logger.debug("fit_tipping_data: min am index: %d", min_am_idx)
    self.logger.debug("fit_tipping_data: max am index: %d", max_am_idx)
    
    # fit the data
    Trx = {}; sigTrx = {}
    tau = {}; sigtau = {}
    for IF in range(4):
      # estimate the parameters
      est_slope = (t[max_am_idx,IF] - t[min_am_idx,IF]) \
                 /   (am[max_am_idx] - am[min_am_idx])
      est_t0 = t[min_am_idx,IF] - est_slope*am[min_am_idx]
      est_tau = est_slope/Tatm
      self.logger.debug("fit_tipping_data: IF%d est. slope = %f", IF, est_slope)
      self.logger.debug("fit_tipping_data: IF%d est. intcp = %f", IF, est_t0)
      # fit data to curve
      popt, pcov = curve_fit(opacity_fitting, e, t[:,IF],
                             p0 = [est_t0, est_tau])
      self.logger.debug("fit_tipping_data: covar is %s", pcov)
      Trx[IF] = popt[0]
      tau[IF] = popt[1]
      sigTrx[IF], sigtau[IF] = numpy.sqrt(numpy.diag(pcov))
      # add fit to legend
      label = "IF%d  %5.1f$\pm$%3.1f  %5.3f$\pm$%5.3f" % (
                 IF, Trx[IF],sigTrx[IF], tau[IF],sigtau[IF])
      handles[IF].set_label(label)
    legend(loc="upper left", numpoints=1)
    if project:
      sessiondir = projects_dir+project+"/Observations/dss43/%4d/%03d/" % (
                                                           self.year, self.DOY)
      fig.savefig(sessiondir+"tipdatafit-"+str(index+1)+".png")
    return Trx, sigTrx, tau, sigtau


examples = """
Examples
========
fill in
"""
if __name__ == "__main__":
  p = initiate_option_parser(__doc__, examples)
  p.usage='dss43-tipping.py [options]'

  p.add_argument('--date',
               dest = 'date',
               type = str,
               default = None,
               help = 'Date of observation as YEAR/DOY string')
  args = p.parse_args()

  mylogger = logging.getLogger()
  # change default console logging
  if get_loglevel(args.console_loglevel) > get_loglevel("info"):
    loglevel = "info"
  else:
    loglevel = args.console_loglevel
  init_logging(mylogger,
                 loglevel = get_loglevel(args.file_loglevel),
                 consolevel = get_loglevel(loglevel),
                 logname = args.logpath+"TidTipAnalyzer.log")
  mylogger.debug("dss43-tipping args: %s", args)

  mylogger = logging.getLogger()
  mylogger.setLevel(logging.DEBUG)
  
  analyzer = TidTipAnalyzer(date=args.date)
  

