"""
Functions for handling Stokes data
"""
import logging
import math

from numpy import array, ndarray, zeros

logger = logging.getLogger(__name__)

def get_source_averages(scandata, sources, first_scan, last_scan, datatype):
  """
  """
  try:
    dindex = scandata['header']['dtype'][3].index(datatype)
  except ValueError:
    logger.error('get_source_averages: data type %s not in data set')
    raise RuntimeError('wrong data set type')
  chankeys = list(scandata.keys())
  chankeys.remove('header')
  chankeys.sort()
  x = {}
  y = {}
  ref = {}
  for chan in chankeys:
    x[chan] = {}
    y[chan] = {}
    ref[chan] = {}
    for source in sources:
      f = first_scan[source]
      l = last_scan[source]
      x[chan][source] = scandata[chan]['freq'] # Hz
      y[chan][source] = zeros((131072,))
      ref[chan][source] = zeros((131072,))
      count = 0
      for scan in range(f,l,2):
        try:
          y[chan][source] += scandata[chan][scan]['data'][:,0,0,dindex]
          ref[chan][source] += scandata[chan][scan+1]['data'][:,0,0,dindex]
          count +=1
        except (KeyError, ValueError):
          # no or invalid data
          logger.debug("get_source_averages: no valid data for chan %s scan %d", chan, scan)
      if count:
        y[chan][source] /= count
        ref[chan][source] /= count
      else:
        y[chan][source] = None
        ref[chan][source] = None
  return x, y, ref

