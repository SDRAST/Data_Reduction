# -*- coding: utf-8 -*-
"""
Modules to support data reduction in Python.
"""
from numpy import array, loadtxt
import logging
import datetime
import scipy.fftpack

module_logger = logging.getLogger(__name__)

def nearest_index(np_array,value):
  """
  Return the index of the element of the array with the nearest value

  Note
  ====
  This has a problem with arrays of datetime.datetime() objects.

  @param np_array : an array of value
  @type  np_array : numpy array

  @param value : a single value
  @type  value : the same type as values in np_array

  @return: index of item in np_array closest to value
  """
  # Convert to to numpy array if necessary
  if type(np_array) == list:
    np_array = array(np_array)
  # find the index of the array value nearest the test value
  if type(value) == datetime.datetime:
    data_array = date2num(np_array)
    ref_value = date2num(value)
    index = abs(data_array-ref_value).argmin()
  else:
    index = abs(np_array-value).argmin()
  # discard points beyond the ends of the array
  if value < np_array[0] or value > np_array[-1]:
    return -1
  else:
    return index

def load_csv_with_header(filename,delimiter=""):
  """
  Takes a text table with headers and converts it into an ASCII array

  There must be a header for each column and the separator must not be
  ambiguous, for example, a space if the column headers also include spaces.

  The data are ASCII and the program or programmer will have to convert
  them to an appropriate form, e.g. int or float or datetime, if necessary.

  A column can be extracted using data[label].

  @param filename : name of the data file
  @type  filename : str

  @param delimiter : column separator
  @type  delimiter : str

  @return: column headers (str), data array(str)
  """
  datafile = open(filename,"r")
  labels = datafile.readline().strip().split(',')
  num_cols = len(labels)
  fmts = ('S20',)*num_cols
  datafile.close()
  data = loadtxt(filename,skiprows=1,delimiter=',',
                 dtype = {'names': tuple(labels),
                          'formats':fmts})
  return labels,data

def unpack_to_complex(rawdata):
  """
  Converts a sequence of alternating real/imag samples to complex

  @param rawdata : alternating real and imaginary bytes
  @type  rawdata : numpy array of signed int8

  @return: numpy array of complex
  """
  datalen = len(rawdata)
  real = rawdata[0:datalen:2]
  imag = rawdata[1:datalen:2]
  data = real + 1j*imag
  return data


def sideband_separate(data):
  """
  Converts a complex array time series and returns two reals with USB and LSB

  This applies a Hilbert transform to the complex data.
  """
  usb = (data.real + scipy.fftpack.hilbert(data).imag)
  lsb = (scipy.fftpack.hilbert(data).real + data.imag)
  return lsb,usb
