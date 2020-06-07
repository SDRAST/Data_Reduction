"""
Analyze tipping curve data taken with observatoryCtrl

A tipping array consists of two 1-D arrays, the first containing the elevation
in degrees and the second the system temperature for each of the four channels.

Older files (like 2014) have only one reading and it may be dBm instead of Tsys

"""
import logging
import numpy
import sys

from collections import OrderedDict
from os.path import basename
from pylab import *
from scipy.optimize import curve_fit

from local_dirs import log_dir
from Radio_Astronomy import dBm_to_watts
from support.logs import init_logging, get_loglevel
from support.options import initiate_option_parser
from support.process import search_response

logdir  = "/home/ops/roach_data/sao_test_data/log_dir/"
K2logs = "K2Observatory"
tiplogs = "tipping_data"

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def sort_by_time(directory, name_template):
  """
  Returns a time-ordered dict with times of files
  
  @param directory : where the files are
  @type  directory : str
  
  @param name_template : part of name common to all files
  @type  name_template : str
  
  @return: dict keyed on filenames
  """
  command = "ls -lrt "+logdir+name_template+'*'
  logger.debug("command: %s", command)
  lines = search_response(['ls','-l', '--full-time',directory],
                          ['grep', name_template])
  files = OrderedDict()
  for line in lines:
    parts = line.strip().split()
    files[parts[-1]] = datestr2num(' '.join(parts[-4:-1]))
  return files

def get_tipping_data(filename):
  """
  Gets tipping data from a file made with observatoryCtrl
  """
  tipfile = open(logdir+filename, 'r')
  data = numpy.load(tipfile)
  tipfile.close()
  elevations = data[0][:].astype(float)
  tsys = []
  for item in data[1].tolist():
    temps = item.ravel()
    if temps[0] > 0:
      # W or K
      tsys.append(temps)
    elif temps[0] < 0:
      # dBm
      power = array(temps)
      tsys.append(dBm_to_watts(power))
  return numpy.array(elevations), numpy.array(tsys)

def plot_tipcurves(e, t, plot_title="", fig_title=""):
  """
  Plots tipping curve data from a specified file
  
  @param tipfiles : dict keyed on filename with file times
  @type  tipfiles : dict
  
  @param index : file in tipfiles to be plotted
  @type  index : int
  """
  #filename = tipfiles.keys()[index]
  #e,t = get_tipping_data(filename)
  for chan in range(len(t[0])):
    plot(e,t[:,chan],label=str(chan+1))
  xlabel("Elevation")
  ylabel("System Temperature")
  grid()
  legend()
  title(plot_title)
  suptitle(fig_title)
  show()

def fit_tipping_data(e, t, Tamb=290):
  """
  Fits tipping curve data from a specified file
  
  @param tipfiles : dict keyed on filename with file times
  @type  tipfiles : dict
  
  @param index : file in tipfiles to be plotted
  @type  index : int
  
  @param Tamb : temperature of the ambient load
  @type  Tamb : float
  
  @return: float (receiver temperature), float (atmospheric opacity)
  """
  def opacity_fitting(x, a, tau):
    x_rad = numpy.deg2rad(x)
    x_sec = 1/(numpy.cos((numpy.pi/2) - x_rad))
    return a + Tamb*(1 - numpy.exp(-tau*x_sec))
  
  def nearest_index(value):
     abse = abs(e-value)
     nearest = abse.argmin()
     return nearest
     
  #filename = tipfiles.keys()[index]
  #e,t = get_tipping_data(filename)
  #emin = e.min()
  #emax = e.max()
  #indices = [nearest_index(emin), nearest_index(emax)]
  #indices.sort()
  #imin,imax = indices
  #imax +=1 # for python style indexing
  imin, imax = extract_simple_tip(e)
  Trx = []
  tau = []
  for chan in range(len(t[0])):
    popt, pcov = curve_fit(opacity_fitting,
                           e[imin:imax], t[imin:imax,chan], p0 = [0, 0])
    Trx.append(popt[0])
    tau.append(popt[1])
  return Trx, tau, imin, imax

def extract_simple_tip(e):
  """
  """
  emin = e.min()
  emax = e.max()
  indices = [nearest_index(emin), nearest_index(emax)]
  indices.sort()
  imin,imax = indices
  imax +=1 # for python style indexing
  return imin, imax

def airmass(elev):
  """
  @param elev : elevation above the horizon in degrees
  @type  elev : float

  @return: float
  """
  if type(elev) == list:
    elev = array(elev)
  return 1/sin(elev*pi/180)

def fit_airmasses(e,t):
  """
  """
  am = airmass(e)
  min_am_index, max_am_index = extract_simple(am)
  slope = []
  intercept = []
  for chan in range(len(t[0])):
    pw = t[chan]
    coefs = polyfit(am[min_am_index:max_am_index],
                    pw[min_am_index:max_am_index], 1)
    slope.append(coefs[0])
    intercept.append(coefs[1])
  return slope, intercept
  
  
if __name__ == "__main__":
  examples = """
  python tipping.py -h|--help
    show help
  """
  # for local testing
  # these data have only one reading, and in dBm not Tsys
  #logdir = "/home/kuiper/Projects/AGN/TAMS/Observations/2014-123/tipping_doy123/"
  # parse the arguments
  from optparse import OptionParser

  class MyParser(OptionParser):
    """
    Subclass of OptionParser which does not mess up examples formatting
    """
    def format_epilog(self, formatter):
      """
      Don't use textwrap; just return the string.
      """
      return self.epilog

  logging.basicConfig(level=logging.WARNING)
  mylogger = logging.getLogger()

  p = initiate_option_parser(__doc__, examples)
  #p = MyParser(epilog=examples)
  #p.set_usage('tipping.py [options]')
  #p.set_description(__doc__)

  opts, args = p.parse_args(sys.argv[1:])
  mylogger = init_logging(mylogger,
                          loglevel = logging.INFO,
                          consolevel = get_loglevel(opts.loglevel),
                          logname = log_dir+"tipping.log")
  
  sys.exit()
  tipfiles = sort_by_time(logdir, tiplogs)
  numfiles = len(list(tipfiles.keys()))
  #index = randint(numfiles)
  for index in range(-20,-2):
    filename = list(tipfiles.keys())[index]
    e,t = get_tipping_data(filename)
    datestr = num2date(tipfiles[filename]).ctime()
    print(("Selected file %d: %s from %s" %
              (index, filename, datestr)))
    print(("min/max elevation = %4.1f/%4.1f, airmass =  %4.1f/%4.1f" % 
          (e.min(), e.max(), airmass(e.min()), airmass(e.max()) ) ))
    figure(index)
    plot_tipcurves(e, t, plot_title=datestr, fig_title=filename)
    savefig(basename(filename)+".png")
    try:
      Trx,tau,imin,imax = fit_tipping_data(e, t)
    except Exception as details:
      print("Fit failed because %s" % details)
    print(("Trx = %s" % Trx))
    print(("tau = %s" % tau))
    print(("from points %d to %d" % (imin,imax)))
    slope, intercept = fit_airmasses(e,t)
    print("airmass slope = %s" % slope)
    print("airmass intercept = %s" % intercept)
    
