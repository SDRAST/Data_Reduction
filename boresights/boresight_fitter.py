"""
This is supposed to be a general purpose boresight fitter but it has too many
DSS-28 dependencies.
"""

import logging
import numpy as NP
import scipy

import Astronomy.Ephem as Aeph
import Data_Reduction.maps as DRm
import Math.least_squares as Mlsq
import support

logger = logging.getLogger(__name__)
dss28 = Aeph.DSS(28)

class ScanFitter(object):
  """
  Create an object to fit a scan in one direction to a baseline and a Gaussian
  
  Public attributes::
  
    atten         - (float) receiver channel attenuation for this can
    baseline_pars - (nparray) polynomial parameters for baseline
    calibrator    - (Astronomy.Ephem.calibrator) calibrator source
    data          - (nparray) VFC count
    ddecs         - (nparray) declination offsets
    direction     - (str) scan axis
    dxdecs        - (nparray) cross-declination offsets
    logger        - (logging.Logger)
    pars          - (nparray) Gaussian parameters
    
  """

  def __init__(self, scan):
    """
    Initiate a scan fitter
    
    This takes a 'scan' object with these attributes::
    
      axis
      datenums
      conv_cfg
      decs
      freq
      ras
      source
      tsys
      
    """
    self.logger = logging.getLogger(logger.name+".ScanFitter")
    # the following returns an ephem planet or quasar
    self.calibrator = Aeph.calibrator(scan.source)
    self.axis = scan.axis
    self.freq = scan.freq
    self.tsys = scan.tsys
    self.dxdecs,self.ddecs = DRm.center_data(scan.date_nums,
                                         scan.ras, scan.decs,
                                         self.calibrator,
                                         dss28)
    
  def fit_gaussian(self, beam_limit=2.5):
    """
    Extract the appropriate data::
        For raster scans, 'xdec' means that 'xdec' stays fixed while the
        antenna moves up and down; 'dec' means that 'dec' stays fixed while the
        left and right.
     The Gaussian is assumed to fit the inner five beamwidths of the data,
     though that limit can be adjusted.  The baseline is the rest of the data,
     although the lower baseline includes at least data[:5] and the upper
     baseline includes data[-5:]
    """
    self.logger.debug("fit_gaussian: direction is %s", self.axis)
    if self.axis.lower() == 'xdec':
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
        lower_baseline = [0, support.nearest_index(x, lower_limit)]
      if upper_limit > x[-5]: # upper baseline segment
        upper_baseline = [-6,-1]
      else:
        upper_baseline = [support.nearest_index(x, upper_limit), -1]
    else:
      # scans go from high sample to low sample
      if upper_limit > x[5]:
        upper_baseline = [0,5]
      else:
        upper_baseline = [0, support.nearest_index(x,upper_limit)]
      if upper_limit < x[-5]:
        upper_baseline = [-6,-1]
      else:
        upper_baseline = [support.nearest_index(x,lower_limit), 0]
    self.logger.debug("fit_gaussian: lower baseline: %s", lower_baseline)
    self.logger.debug("fit_gaussian: upper baseline: %s", upper_baseline)
    # define the baseline data
    xdata = NP.append(x[lower_baseline[0]:lower_baseline[1]],
                      x[upper_baseline[0]:upper_baseline[1]]).astype(float)
    ydata = NP.append(self.tsys[lower_baseline[0]:lower_baseline[1]],
                      self.tsys[upper_baseline[0]:upper_baseline[1]]).astype(float)
    #   Fit baseline
    self.baseline_pars = scipy.polyfit(xdata,ydata,1)
    self.logger.debug("fit_gaussian: baseline parameters: %s", self.baseline_pars)
    #   Fit the beam
    zdata = NP.array(self.tsys).astype(float)
    self.logger.debug("fit_gaussian: zdata: %s", zdata)
    height = zdata[beam_index] - scipy.polyval(self.baseline_pars, x[beam_index])
    self.logger.debug("fit_gaussian: height: %s", height)
    sigma = Mlsq.st_dev(beamwidth)
    initial_guess = [height, beam_center, sigma]
    # in this case we only fit out to one beamwidth
    if x[0] < x[-1]:
      xfit =  x[support.nearest_index(x,beam_center-beamwidth):support.nearest_index(x,beam_center+beamwidth)]
      y = zdata[support.nearest_index(x,beam_center-beamwidth):support.nearest_index(x,beam_center+beamwidth)]
    else:
      xfit =  x[support.nearest_index(x,beam_center+beamwidth):support.nearest_index(x,beam_center-beamwidth)]
      y = zdata[support.nearest_index(x,beam_center+beamwidth):support.nearest_index(x,beam_center-beamwidth)]
    self.pars, err = Mlsq.fit_gaussian(Mlsq.gaussian_error_function,
                            initial_guess,
                            xfit,
                            y-scipy.polyval(self.baseline_pars,xfit))
    return self.baseline_pars, self.pars, err
    
