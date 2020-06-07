"""
Finds Tidbinbilla tipping curve files
 
Notes
=====

Instead of 'TidTipFinder' having an initialization argument 'remote', one could
just test with 'socket.gethostname' and do a remote mount if the response is
not 'crux'.
"""
import calendar
import glob
import logging
import numpy
import os
import time

from collections import OrderedDict
from pylab import *

from DatesTimes import calendar_date
from support import nearest_index
#from support.tunneling import RemoteDir

logger = logging.getLogger(__name__)

saologdir  = "/home/ops/roach_data/sao_test_data/log_dir/" # on crux

IFcolors = ['b', 'g', 'r', 'c']

class TidTipFinder(object):
  """
  Attributes::
    am
    directory
    files
    logger
    tsys
  """
  def __init__(self, dss=43, remote=True, date=None):
    """
    """
    self.logger = logging.getLogger(logger.name+".TidTipFinder")
    if remote:
      #mp = RemoteDir("crux")
      self.directory = os.environ['HOME'] + "/mnt/crux" + saologdir
    else:
      self.directory = saologdir
    self.files = self.sort_by_time(subdir="2014-2016") # from subdirec'y first
    self.files.update(self.sort_by_time())
    if self.files:
      pass
    else:
      self.logger.debug("__init__: found %d files; crux tunnel open?", len(self.files))
    if date:
      pass
    else:
      date = input("Session date (YYYY/DDD)? ")
    self.date = date
    yearstr, doystr = date.split('/')
    self.year = int(yearstr)
    self.DOY = int(doystr)
    if self.files:
      self.get_session_files(self.date)
      

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
    keys = list(filedict.keys())
    keys.sort()
    ordered = OrderedDict()
    for key in keys:
      ordered[key] = filedict[key]
    return ordered

  def get_session_files(self, date):
    """
    Get the files for the specified data
    """
    year,month,day = calendar_date(self.year, self.DOY)
    UNIXstart = calendar.timegm((self.year,month,day,  0, 0, 0, 0,0))
    UNIXend   = calendar.timegm((self.year,month,day, 23,59,59, 0,0))
    filekeys = list(self.files.keys())
    if filekeys:
      pass
    else:
      raise RuntimeError("no tipping files found")
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
    files.sort()
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

