import cPickle
import h5py
import logging
import numpy

from pylab import *

from DatesTimes import UnixTime_to_MPL

logger = logging.getLogger(__name__)

class KurtspecMonitor(object):
  """
  class to examine spectra of 1-sec kurtspec packet averages
  """
  def __init__(self, filename=None):
    """
    initialize a KurtspecMonitor object and optionally open a file
    
    If a filename is not given, use 'open_datafile()'
    
    @param filename : name of an HDF5 monitor file
    """
    if filename:
      data = self.open_datafile(filename)
    pass
  
  def open_datafile(self, filename):
    """
    open an HDF5 datafile
    """
    data = cf5py.Fil(filename)
    self.data_array = data.reshape(data.shape + (1,))
    self.signal = self.data_array.attrs['signal']
    self.time = self.data_array['time']
    self.power = self.data_array['power']
    self.kurtosis = self.data_array['kurtosis']
    return self.data_array
    
  def average_power(self, log=False):
    """
    make dynamic spectra averaged over one second
    """
    mean_pwr = self.power.mean()
      if log:
        return numpy.log10(mean_pwr)
      else:
        return mean

  def plot_spectra(self, log=True):
    """
    Note that this default differs from the one in 'average_power'
    """
    avg_pwr = self.average_power(log=log)
    mpldates = UnixTime_to_MPL(self.time)
    fig = figure()
    plot_date(mpldates, avg_pwr)
    legend()
    grid()
    title("Average Passband")
    fig.autofmt_xdate()
    show()

  def show_spectrogram(self, ptype='power'):
    """
    display a spectrogram of the original data
    """
    fig = figure(figsize=(30,12))
    if ptype == 'power':
      imshow(self.power)
    else:
      imshow(self.kurtosis)
    colorbar()
    show()
    
if __name__ == "__main__":
    mon = KurtspecMonitor()
    
