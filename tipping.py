"""
Analyze tipping curve data taken with observatoryCtrl

This was superceded by class TidTipAnalyzer in Data_Reduction.dss43-tipping

A tipping array consists of two 1-D arrays, the first containing the elevation
in degrees and the second the system temperature for each of the four channels.

Older files (like 2014) have only one reading and it may be dBm instead of Tsys

"""
import logging
import math
import numpy

from scipy.optimize import curve_fit

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def airmass(elev):
  """
  @param elev : elevation above the horizon in degrees
  @type  elev : float

  @return: float
  """
  if type(elev) == list:
    elev = numpy.array(elev)
  return 1/numpy.sin(elev*math.pi/180)

def fit_tipcurve_data(elev, Tsys, Tatm=250, linear=True):
  """
  fits system temperatures to airmass
  
  The length of 'elev' and the first axis of 'Tsys' must be the same
  
  If the Tsys array has two or more dimensions then a fit is done on the data
  over the first index for all the values of the other indices and the results
  will have the same shape as the indices 1..N
  
  If Tatm is given, the radiative transfer equation will be used and K/am will
  be Tatm*tau(zenith).  If it is not given, the straight-line low-tau equation
  is used and Tatm*tau(zenith) will be the slope of the line.
  
  @param elev : 1D array of elevations
  @type  elev : list or 1D numpy.array of float
  
  @param Tsys : ND array of system temperatures
  @type  Tsys : list or numpy.array of float
  
  @param Tatm : mean air temperature along path
  @type  Tatm : float
  
  @param linear : use the linear (low tau) approximation
  @type  linear : True
  """
  def est_pars(am, t, Tatm=None):
    """
    estimate the fit parameters
    
    @param am : air mass
    @type  am : numpy 1D array of float
    
    @param t : system temperature or power
    @type  t : numpy 1D array of float
    
    @param Tatm : air temperature along path for zenith tau
    @type  Tatm : int or float
    """
    # estimate the parameters
    est_slope = (t[max_am_idx] - t[min_am_idx]) \
              /   (am[max_am_idx] - am[min_am_idx])
    est_t0 = t[min_am_idx] - est_slope*am[min_am_idx]
    if Tatm:
      est_tau = est_slope/Tatm
      return est_t0, est_tau
    else:
      return est_t0, est_slope
              
  def rad_transf_eqn(am, Tsysvac, tau):
    """
    radiative transfer equation for atmosphere
      
    Tatm is the average atmospheric temperature along the line of sight
      
    @param am : airmass
    @type  am : numpy.array of float
      
    @param Tsysvac : system tenperature in a vacuum
    @type  Tsysvac : float
      
    @param tau : opacity
    @type  tau : float
    """
    return Tsysvac + Tatm*(1 - numpy.exp(-tau*am))
  
  def lin_transf_eqn(am, Pvac, Ppam):
    """
    linear (small tau) radiative transfer equation for atmosphere
    
    @param am : airmass
    @type  am : numpy.array of float
      
    @param Pvac : system power in a vacuum
    @type  Pvac : float
      
    @param Ppam : added power per airmass
    @type  Ppam : float
    """
    return Ppam + Ppam*am
    
  # find the indices of the smallest and largest airmass
  if type(elev) == list:
    elev = numpy.array(elev)
  am = airmass(elev)
  min_am_idx = am.argmin()
  max_am_idx = am.argmax()
  logger.debug("fit_tipcurve_data: min am at %d, max am at %d",
               min_am_idx, max_am_idx)
  if type(Tsys) == list:
    Tsys = numpy.array(Tsys)
  # fit data to small tau function
  if linear:
    # estimate the fit parameters
    est_t0, est_slope = est_pars(am, Tsys)
    logger.debug("fit_tipcurve_data: estimated Tsys(0)=%f, slope=%f",
               est_t0, est_slope)
    popt, pcov = numpy.polyfit(am, Tsys, 1, cov=True)
  else:
    # estimate the fit parameters
    est_t0, est_tau = est_pars(am, Tsys, Tatm=Tatm)
    logger.debug("fit_tipcurve_data: estimated Tsys(0)=%f, tau(1)=%f",
               est_t0, est_tau)
    # fit the data
    popt, pcov = curve_fit(lin_transf_eqn, am, Tsys, p0 = [est_t0, est_tau])
  logger.debug("fit_tipcurve_data: popt = %s", popt)
  logger.debug("fit_tipcurve_data: pcov = %s", pcov)
  if linear:
    Prx = popt[1]    # noise power from receiver in above atmosphere
    Ppam = popt[0]   # atmosphere power contribution per airmass
    sigPpam, sigPrx = numpy.sqrt(numpy.diag(pcov))
    return Prx, sigPrx, Ppam, sigPpam
  else:  
    Trx = popt[0]
    tau = popt[1]
    sigTrx, sigtau = numpy.sqrt(numpy.diag(pcov))
    return Trx, sigTrx, tau, sigtau
  
  
