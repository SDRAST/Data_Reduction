"""
Analyze tipping curve data taken with observatoryCtrl

This was superceded by class TidTipAnalyzer in Data_Reduction.dss43-tipping

A tipping array consists of two 1-D arrays, the first containing the elevation
in degrees and the second the system temperature for each of the four channels.

Older files (like 2014) have only one reading and it may be dBm instead of Tsys

"""
import logging
import numpy

from pylab import *

K2logs = "K2Observatory"
tiplogs = "tipping_data"

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
  
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
