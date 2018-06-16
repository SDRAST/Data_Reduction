"""
Convert raw kurtosis spectrometer data into 1 sec averages
"""
import cPickle
import h5py
import numpy

data = h5py.File("../145/kurt-2018-145-192439.hdf5")

def merge_1_sec(data, group):
  """
  Average power and kurtosis data for one second
  """
  one_second = data[group]
  for key in ['pwr-I', 'krt-I', 'pwr-Q', 'krt-Q']:
    array1d = numpy.array(one_second[key]).mean(axis=0)
    array2d = array1d.reshape(array1d.shape[0],1)
    if key == 'pwr-I':
      merged = array2d
    else:
      merged = numpy.append(merged, array2d, axis=1)
  return merged

def merge_monitor_data(data, groupkeys):
  """
  collect 1-sec averages
  """
  groupkeys = data.keys()
  groupkeys.sort()
  for group in groupkeys:
    merged = merge_1_sec(data, group)
    if group == groupkeys[0]:
      bigmerge = merged.reshape(merged.shape+(1,))
    else:
      bigmerge = numpy.append(bigmerge, merged.reshape(merged.shape+(1,)), axis=2)
  return bigmerge
