import cPickle
import glob
import logging
import numpy
import socket

from pylab import *

from local_dirs import wvsr_dir

logger = logging.getLogger(__name__)

datapath = "/data/HDF5/dss14/2018/"
filename = datapath + "169/mon-2018-169-004039.bin"

def get_data_path(workstation, dss, year, doy):
  """
  """
  thishost = socket.gethostname()
  if thishost == workstation:
    datapath = wvsr_dir
  else:
    workstationpath = "/home/kuiper/mnt/"+workstation
    if thishost == 'dto':
      datapath = workstationpath+wvsr_dir
    else:
      dtopath = "/home/kuiper/mnt/dto"
      if thishost == 'kuiper':
        datapath = dtopath +workstationpath+ wvsr_dir
      else:
        return None
  datapath += "HDF5/dss"+str(dss)+"/"+str(year)+("/%03d/" % doy)
  return datapath

class KurtspecMonitor(object):
  """
  class to examine spectra of 1-sec kurtspec packet averages
  """
  array_index = {"power":    {"I": 0,
                               "Q": 2},
                 "kurtosis": {"I": 1,
                               "Q": 3}}

  def __init__(self, filename=None):
    """
    """
    if filename:
      data_array = self.open_datafile(filename)
    pass
  
  def open_datafile(self, filename):
    """
    """
    dfile = open(filename, "rb+")
    data = cPickle.load(dfile)
    self.data_array = data.reshape(data.shape + (1,))
    count = 0
    reading = True
    while reading:
      try:
        data = cPickle.load(dfile)
        self.data_array = numpy.append(self.data_array, 
                                       data.reshape(data.shape + (1,)),
                                       axis=2)
        count += 1
        print "\r"+str(count),
      except EOFError:
        reading = False
    dfile.close()
    return self.data_array

  def extract_dynamic_spectra(self):
    """
    """
    self.dyn_spec = {}
    for moment in ["power","kurtosis"]:
      self.dyn_spec[moment] = {}
      for inpt in ["I","Q"]:
        index = KurtspecMonitor.array_index[moment][inpt]
        if moment == "power":
          self.dyn_spec[moment][inpt] = numpy.log10(self.data_array[:,index,:])
        else:
          self.dyn_spec[moment][inpt] = self.data_array[:,index,:]
    return self.dyn_spec
    
  def average_power(self, log=False):
    """
    """
    avg_pwr = {}
    for inpt in ["I","Q"]:
      index = KurtspecMonitor.array_index["power"][inpt]
      mean = self.data_array[:,index,:].mean(axis=1)
      if log:
        avg_pwr[inpt] = numpy.log10(mean)
      else:
        avg_pwr[inpt] = mean
    return avg_pwr

  def plot_spectra(self, log=True):
    """
    Note that this default differs from the one in 'average_power'
    """
    avg_pwr = self.average_power(log=log)
    figure()
    plot(avg_pwr['I'], label="I")
    plot(avg_pwr['Q'], label="Q")
    legend()
    grid()
    title("Average Passband")
    show()

  def show_spectrogram(self, type='power'):
    """
    """
    self.extract_dynamic_spectra()
    fig, ax = subplots(nrows=2, ncols=1, figsize=(30,12))
    ax[0].imshow(self.dyn_spec[type]['I'])
    title(type+" I")
    colorbar()
    ax[1].imshow(self.dyn_spec[type]['Q'])
    title(type+" Q")
    colorbar()
    show()
    
if __name__ == "__main__":
  workstation = "gpu2"
  dss=14
  year=2018
  doy=179
  datapath = get_data_path(workstation, dss, year, doy)
  files = glob.glob(datapath+"*.pkl")
  mon = KurtspecMonitor(files[0])
    
