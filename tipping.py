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

class AtmosAnalyzer(object):
  """
  Analyzer for environmental data contained in a weather data cube
  
  Weather data cubes are created with the DSNFITSexaminer method 
  get_wx_datacubes(), or the SessionAnalyzer method get_good_weather_data().
  A weather data cube is a dict with keys 'TAMBIENT', 'WINDDIRE', 'UNIXtime',
  'TSYS', 'HUMIDITY', 'PRESSURE', 'ELEVATIO', 'WINDSPEE'.  The data asociated
  with each key is a dict with numpy array for (SIG state) True and for False.
  The 'TSYS' array has four axes representing::
        time index   - 0-based sequence in order of matplotlib datenum 
        subchannel   - CYCLE value
        beam         - 1-based number sequence
        IF           - 1-based number sequence, usually representing pol
  The other keys have only a time axis.
  """
  def __init__(self, datacube):
    """
    """
    self.logger = logging.getLogger(logger.name+".AtmosAnalyzer")
    self.data = datacube

  def fit_Tsys_to_airmass(self, Tatm=250, linear=True):
    """
    Fit tipping curve data implicit in sessions elev and Tsys data
    
    Returns numpy arrays with indices for sig/ref state, subchannel, beam, IF.
    
    When linear=True then a straight line fit is performed. This would be
    appropriate when the TSYS units are not kelvin but count or something else.
    The parameters are zero airmass power intercept and its standard deviation,
    and power per airmass and its standard deviation.
    
    If the Tsys units are in K, then a linear=False fit is appropriate, which
    fits the data to a radiative transfer model. The returned parameters are
    then the system temperature above the atmosphere, its standard deviation,
    and the optical depth per airmass and its standard deviation. The average
    of the physical temperature of the atmosphere defaults to 250 K.
        
    @param Tatm : air temperature along line of sight
    @type  Tatm : float
  
    @param linear : use the linear (low tau) approximation
    @type  linear : True
    """
    # get Tsys data structure  
    if 'TSYS' in self.data:
      nrows, num_cy, num_bm, num_pl = self.data['TSYS'][True].shape
    else:
      self.logger.error("fit_Tsys_to_airmass: no system temperature data")
      return None
    states = list(self.data['TSYS'].keys())
    num_st = len(states)
    
    if 'ELEVATIO' in self.data:
      pass
    else:
      self.logger.error("fit_Tsys_to_airmass: no elevation data")
      return None
    
    # fit the data
    param_shape = (num_st, num_cy, num_bm, num_pl)
    # what type of fit?
    if linear:
      keys = ["pwr0", "pwr0stdv", "slope", "slopestdv"]
    else:
      keys = ["Tsys0", "Tsys0stdev", "tau0", "tau0stdev"]
    # intialize
    parameters = {}
    for key in keys:
      parameters[key] = numpy.zeros(param_shape)
                    
    for sig in states:              # sig/ref state first
      for subch in range(num_cy):   # subchannel second
        for beam in range(num_bm):  # beam third
          for pol in range(num_pl): # pol fourth
            el = self.data['ELEVATIO'][sig]
            Tsys = self.data['TSYS'][sig][:,subch,beam,pol]
            pars = fit_tipcurve_data(el, Tsys, Tatm=Tatm, linear=linear)
            for key in keys:
              index = keys.index(key)
              parameters[key][sig, subch, beam, pol] = pars[index]
    return parameters
    
def airmass(elev):
  """
  @param elev : elevation above the horizon in degrees
  @type  elev : float

  @return: float
  """
  if type(elev) == list:
    elev = numpy.array(elev)
  return 1/numpy.sin(elev*math.pi/180)

def fit_tipcurve_data(elev, Tsys, Tatm=250, method="linear"):
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
    
  def transfer_appr(e, Tsysvac, tau):
    """
    Taylor expansion approximation of transfer equation
    """
    am = airmass(e)
    return Tsysvac + Tatm*tau*am - 0.5*(Tatm*tau*am)**2
  
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
    
  if type(elev) == list:
    elev = numpy.array(elev)
  if type(Tsys) == list:
    Tsys = numpy.array(Tsys)
  # find the indices of the smallest and largest airmass
  am = airmass(elev)
  min_am_idx = am.argmin()
  max_am_idx = am.argmax()
  logger.debug("fit_tipcurve_data: min am at %d, max am at %d",
               min_am_idx, max_am_idx)
    
  # estimate the parameters
  if method == 'exact' or method == 'quadratic':
    est_t0, est_tau = est_pars(am, Tsys, Tatm=Tatm)
    logger.debug("fit_tipcurve_data: estimated Tsys(0)=%f, tau(1)=%f",
               est_t0, est_tau)
  else:
    est_t0, est_slope = est_pars(am, Tsys)
  # perform the fit
  if method == 'linear':
    popt, pcov = numpy.polyfit(am, Tsys, 1, cov=True)
    logger.debug("fit_tipcurve_data: popt = %s", popt)
    logger.debug("fit_tipcurve_data: pcov = %s", pcov)
    Prx = popt[1]    # noise power from receiver in above atmosphere
    Ppam = popt[0]   # atmosphere power contribution per airmass
    sigPpam, sigPrx = numpy.sqrt(numpy.diag(pcov))
    return Prx, sigPrx, Ppam, sigPpam
  elif method == 'quadratic':
    popt, pcov = numpy.polyfit(am, Tsys, 2, cov=True)
    logger.debug("fit_tipcurve_data: popt = %s", popt)
    logger.debug("fit_tipcurve_data: pcov = %s", pcov)
    Trx = popt[2]
    tau = 2*popt[0]/popt[1]
    Tatm = popt[1]/tau
    sigP0, sigP1, sigP2 = numpy.sqrt(numpy.diag(pcov))
    # the derivation for following is in the log for 2017 Dec 19
    sigTrx = sigP0
    sigtau = math.sqrt(tau**2*(sigP1**2 + sigP2**2/popt[2]**2))
    sigTatm = math.sqrt(((1+popt[1]**2)*sigP1**2
                        + popt[1]**2*sigP2**2/sigP1**2)/tau**2)
    return Trx, sigTrx, tau, sigtau, Tatm, sigTatm
  else:
    popt, pcov = curve_fit(lin_transf_eqn, am, Tsys, p0 = [est_t0, est_tau])
    logger.debug("fit_tipcurve_data: popt = %s", popt)
    logger.debug("fit_tipcurve_data: pcov = %s", pcov)
    Trx = popt[0]
    tau = popt[1]
    sigTrx, sigtau = numpy.sqrt(numpy.diag(pcov))
    return Trx, sigTrx, tau, sigtau
  
  
