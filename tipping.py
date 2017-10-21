"""
Analyze tipping curve data taken with observatoryCtrl

A tipping array consists of two 1-D arrays, the first containing the elevation
in degrees and the second the system temperature for each of the four channels.

Older files (like 2014) have only one reading and it may be dBm instead of Tsys

"""
import logging
import numpy

from collections import OrderedDict
from pylab import *
from scipy.optimize import curve_fit

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
    tsys.append(item.ravel())
  return elevations, numpy.array(tsys)

def plot_tipcurves(tipfiles, index):
  """
  Plots tipping curve data from a specified file
  
  @param tipfiles : dict keyed on filename with file times
  @type  tipfiles : dict
  
  @param index : file in tipfiles to be plotted
  @type  index : int
  """
  filename = tipfiles.keys()[index]
  e,t = get_tipping_data(filename)
  for chan in range(4):
    plot(e,t[:,chan],label=str(chan+1))
  xlabel("Elevation")
  ylabel("System Temperature")
  grid()
  legend()
  title(num2date(tipfiles[filename]).ctime())
  show()

def fit_tipping_data(tipfiles, index, Tamb=290):
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
  def opacity_fitting(x, a, Tatm, tau):
    x_rad = numpy.deg2rad(x)
    x_sec = 1/(numpy.cos((numpy.pi/2) - x_rad))
    return a + Tatm*(1 - numpy.exp(-tau*x_sec))
    
  filename = tipfiles.keys()[index]
  e,t = get_tipping_data(filename)
  popt, pcov = curve_fit(opacity_fitting, e, t[0], p0 = [0, 0])
  Trx = popt[0]
  tau = popt[1]
  return Trx, tau

def airmass(elev):
  """
  @param elev : elevation above the horizon in degrees
  @type  elev : float

  @return: float
  """
  if type(elev) == list:
    elev = array(elev)
  return 1/sin(elev*pi/180)
  
if __name__ == "__main__":
  # for local testing
  # these data have only one reading, and in dBm not Tsys
  logdir = "/home/kuiper/Projects/AGN/TAMS/Observations/2014-123/tipping_doy123/"
  
  tipfiles = sort_by_time(logdir, tiplogs)
  numfiles = len(tipfiles.keys())
  index = randint(numfiles)
  filename = tipfiles.keys()[index]
  logger.info("Selected file %d: %s from %s",
              index, filename, num2date(tipfiles[filename]).ctime())
  plot_tipcurves(tipfiles, index)
  print fit_tipping_data(tipfiles, index)
