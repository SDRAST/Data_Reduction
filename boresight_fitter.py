import logging
from math import log10, pi
from matplotlib.dates import date2num
import numpy as NP
from scipy import polyfit, polyval

from Data_Reduction.maps import center_data
from DatesTimes import UnixTime_to_MPL
from Math.least_squares import fit_gaussian, gaussian_error_function, st_dev
from Radio_Astronomy import HPBW
from support import nearest_index
from support.Ephem import calibrator, DSS, EphemException

logger = logging.getLogger(__name__)
dss28 = DSS(28)
beam_limit = 2.5 # beam widths

def DSS28_beamtaper(freq):
  """
  ad hoc fit to beamwidth vs frequency plot
  """
  if freq < 7:
    taper=0
  else:
    taper = 50*(log10(freq)-log10(7))
  return taper
  
def DSS28_beamwidth(freq):
  """
  beamwidth in deg. with edge taper
  """
  return HPBW(DSS28_beamtaper(freq), 0.3/float(freq), 34)*180/pi

class ScanFitter(object):
  """
  Create an object to fit a scan in one direction to a baseline and a Gaussian
  """
  def __old__(self, date_nums, ras, decs, data, source, direction, freq):
    """
    Initiate a scan fitter
    """
    self.logger = logging.getLogger(logger.name+".ScanFitter")
    self.calibrator = calibrator(source)
    self.direction = direction
    self.freq = freq
    self.data = data
    self.dxdecs,self.ddecs = center_data(date_nums,
                                         ras, decs,
                                         self.calibrator,
                                         dss28)
    #self.logger.debug("__init__: delta-X-dec: %s", self.dxdecs)
    #self.logger.debug("__init__: delta-dec: %s", self.ddecs)

  def __init__(self, scan):
    """
    """
    self.logger = logging.getLogger(logger.name+".ScanFitter")
    self.calibrator = calibrator(scan.source)
    self.direction = scan.axis
    self.freq = scan.freq
    self.data = scan.tsys
    self.atten = scan.conv_cfg['atten']
    self.dxdecs,self.ddecs = center_data(scan.date_nums,
                                         scan.ras, scan.decs,
                                         self.calibrator,
                                         dss28)
    
  def fit_gaussian(self):
    # Extract the appropriate data
    #    For raster scans, 'xdec' means that 'xdec' stays fixed while the
    #    antenna moves up and down; 'dec' means that 'dec' stays fixed while the
    #    left and right
    self.logger.debug("fit_gaussian: direction is %s", self.direction)
    if self.direction.lower() == 'xdec':
      x = NP.array(self.ddecs)
    else:
      x = NP.array(self.dxdecs)
    self.logger.debug("fit_gaussian: selected x: %s", x)
    # define the domain of the Gaussian fit:
    beam_index = NP.array(self.data).argmax()
    self.logger.debug("fit_gaussian: peak at index %d", beam_index)
    beam_center = x[beam_index]
    self.logger.debug("fit_gaussian: peak at x = %f", beam_center)
    beamwidth = DSS28_beamwidth(self.freq/1000)
    self.logger.debug("fit_gaussian: beamwidth = %f deg", beamwidth)
    lower_limit  = beam_center - beam_limit*beamwidth
    upper_limit = beam_center + beam_limit*beamwidth
    self.logger.debug("fit_gaussian: center left: %f", lower_limit)
    self.logger.debug("fit_gaussian: center right: %f", upper_limit)
    # Define baseline ranges for the lower end and the upper end of the spectrum
    #  * 'lower_baseline' and 'upper_baseline' are 2-item lists
    #  * assume that there are at least 5 data points for each baseline section
    if x[0] < x[-1]: # increasing X-coordinate
      # scans go from low sample to high sample
      if lower_limit < x[5]: # lower baseline segment
        lower_baseline = [0,5]
      else:
        lower_baseline = [0, nearest_index(x, lower_limit)]
      if upper_limit > x[-5]: # upper baseline segment
        upper_baseline = [-6,-1]
      else:
        upper_baseline = [nearest_index(x, upper_limit), -1]
    else:
      # scans go from high sample to low sample
      if upper_limit > x[5]:
        upper_baseline = [0,5]
      else:
        upper_baseline = [0, nearest_index(x,upper_limit)]
      if upper_limit < x[-5]:
        upper_baseline = [-6,-1]
      else:
        upper_baseline = [nearest_index(x,lower_limit), 0]
    self.logger.debug("fit_gaussian: lower baseline: %s", lower_baseline)
    self.logger.debug("fit_gaussian: upper baseline: %s", upper_baseline)
    # define the baseline data
    xdata = NP.append(x[lower_baseline[0]:lower_baseline[1]],
                      x[upper_baseline[0]:upper_baseline[1]]).astype(float)
    ydata = NP.append(self.data[lower_baseline[0]:lower_baseline[1]],
                      self.data[upper_baseline[0]:upper_baseline[1]]).astype(float)
    #   Fit baseline
    self.baseline_pars = polyfit(xdata,ydata,1)
    self.logger.debug("fit_gaussian: baseline parameters: %s", self.baseline_pars)
    #   Fit the beam
    zdata = NP.array(self.data).astype(float)
    self.logger.debug("fit_gaussian: zdata: %s", zdata)
    height = zdata[beam_index] - polyval(self.baseline_pars, x[beam_index])
    self.logger.debug("fit_gaussian: height: %s", height)
    sigma = st_dev(beamwidth)
    initial_guess = [height, beam_center, sigma]
    # in this case we only fit out to one beamwidth
    if x[0] < x[-1]:
      xfit =  x[nearest_index(x,beam_center-beamwidth):nearest_index(x,beam_center+beamwidth)]
      y = zdata[nearest_index(x,beam_center-beamwidth):nearest_index(x,beam_center+beamwidth)]
    else:
      xfit =  x[nearest_index(x,beam_center+beamwidth):nearest_index(x,beam_center-beamwidth)]
      y = zdata[nearest_index(x,beam_center+beamwidth):nearest_index(x,beam_center-beamwidth)]
    self.pars, err = fit_gaussian(gaussian_error_function,
                            initial_guess,
                            xfit,
                            y-polyval(self.baseline_pars,xfit))
    return self.baseline_pars, self.pars, err
    
