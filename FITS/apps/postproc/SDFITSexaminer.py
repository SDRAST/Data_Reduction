"""
Provides class to examine DSN SDFITS data files
"""
import datetime
import logging

from numpy import array, ndarray

from Astronomy.redshift import doppler_radio, doppler_optical, doppler_relat
from MonitorControl.BackEnds import get_freq_array
from support import mkdir_if_needed, nearest_index


logger = logging.getLogger(__name__)

class Data(object):
  """
  A subset of data extracted from an average difference spectrum
  
  Attributes::
    channels -
    frame    -
    logger   -
    x        -
    y        -
  """
  def __init__(self, x, y):
    """
    @param x : X-values
    @type  x : list or nparray
    """
    self.logger = logging.getLogger(logger.name+".Data")
    if type(x) == list:
      self.x = array(x)
    elif type(x) == ndarray:
      self.x = x
    else:
      self.logger.error("__init__: type %s is not valid for X", type(x))
    #if type(y) == dict:
    #  if y.keys() == [0, 1]:
    #    self.y = y
    #  else:
    #    self.logger.error("__init__: pols keys must be [0,1]")
    #elif type(y) == list:
    #  if type(y[0]) == list:
    #    self.y = array(y)
    #  else:
    #    self.logger.error("__init__: %s is not a dict of valid Y-values")
    self.y = y
    self.frame = None
    self.channels = None

class DSNFITSexaminer(object):
  """
  Class for examining SAOdataset objects
  """
  def __init__(self, FITSextension=None, dss=None):
    """
    Create an SAOexaminer object
    
    Attributes::
      avg_diff     - average of all scans for each pol
      datafiles    - list of pkl files in datapath
      datapath     - location of project data for this observation session
      dataset      - collection of scana and metadata taken from HDF5 file
      logger       - logging.Logger object
      projworkpath - obs. session sub-dir. in project directory
      rawdatapath  - local of HDF5 files for this obs. session
      tau          - total integration times for averaged spectra
    
    Methods::
      align_spectrum
      beam_switched_average - average spectra for all the scan
      load_ephemeris        - get the ephemeris for the designated source
    
    @param FITSextension : SDFITS file BINTABLE extension
    @type  FITSextension : pyfits HDU
    """
    self.logger = logging.getLogger(logger.name+'.DSNFITSexaminer')
    if FITSextension and dss:
      self.avg_diff = {}
      self.dataset = FITSextension.data
      self.nscans, self.nbeams, self.nrecs, self.npols, nra, ndec, self.nchans \
          = self.dataset['SPECTRUM'].shape
      self.header = FITSextension.header
      self.tau = {}      # cumulative integration time
      self.frame = None
      self.dss = dss
      # the following is a HACK until all the SDFITS files are fixed
      try:
        self.SB = self.dataset['SIDEBAND']
      except KeyError:
        if self.dataset['CDELT1'][0] > 0:
          self.SB = 'U'
        else:
          self.SB = 'L'
    else:
      self.logger.error("__init__: extension and dss are required")
      raise RuntimeError("insufficient arguments")
  
  def extract_window(self, row, xlimits, frame="RADI-LSR"):
    """
    Extracts a subset of a avg_diff spectrum
    """
    if frame:
      if frame == "RADI-LSR":
        x = self.compute_X_axis("RADI-LSR")
      elif frame == "RADI-OBJ":
        if re.match('67P', self.dataset['OBJECT'][0]):
          vspline = self.load_ephemeris('67P')
          x = self.compute_X_axis("RADI-OBJ", vspline=vspline)
        else:
          self.logger.error("extract_window: no ephemeris for %s", self.object)
      else:
        x = self.compute_X_axis(frame)
    else:
      frame = self.frame
    self.logger.debug("extract_window: x=%s", x)
    self.frame = frame
    # extract the designated range
    x1,x2 = xlimits
    if x[0] < x[-1]:
      ch1 = nearest_index(x, x1)
      ch2 = nearest_index(x, x2)
    else:
      ch1 = nearest_index(x, x2)
      ch2 = nearest_index(x, x1)
    #for pol in self.avg_diff.keys():
    #  extract[pol] = self.avg_diff[pol][ch1:ch2]
    d = Data(x[ch1:ch2], self.dataset['SPECTRUM'][row][:,:,:,0,0,ch1:ch2])
    d.frame = frame
    d.channels = (ch1,ch2)
    return d
  
  def compute_X_axis(self, frame, ref_freq=None, vspline=None, 
                     obstime=None):
    """
    Computes the appropriate X-axis for the averaged difference spectrum.
     
    The header parameter 'velocity' is set to the value of 'vobj' which is 
    velocity of the object in the rest velocity of frame specified.
    
    Acceptable frames are defined in the SAOhdf5.rel_freq_units() docstring.
    In addition we allow here RADI-OBJ which defines the rest frame of the
    object.
    
    @param frame : the rest frame and the X-axis type (freq or vel)
    @type  frame : str
    
    @param ref_freq : frequency in MHz for computing relative frequencies
    @type  ref_freq : float
    
    @param vspline : radial velocity of a moving body as a function of time
    @type  vspline : function
    
    @param obstime : the time at which the spline is to be evaluated
    @type  obstime : UNIX seconds):
    """
    if not obstime:
      index = (0,0,0,0,0,0,0)
      obstime = datetime.datetime.fromtimestamp(
                                       self.dataset['UNIXtime'][index])
    self.logger.debug("compute_X_axis: using time %s", obstime)
    if ref_freq:
      f_ref = ref_freq
    else:
      f_ref = self.dataset['RESTFREQ'][0]/1e6 # MHz
    v_ref = self.dataset['VELOCITY'][0]
    if frame:
      self.frame = frame
    self.logger.debug("compute_X_axis: requested frame is %s", frame)
    self.logger.debug("compute_X_axis: reference frequency is %10.3f", f_ref)
    if frame == "CHAN-OBS" or frame == "FREQ-OBS" or frame == "RELA-OBS":
      # all these are in the observer's rest frame
      x = self.rel_freq_units(frame=frame, ref_freq=f_ref)
      vobj = None
    elif frame == "FREQ-LSR":
      vobj = self.V_LSR(obstime, self.dss) # km/s
      C = c/1000 # km/s
      delf = -self.V_LSR(obstime, self.dss)*f_ref/C
      self.logger.debug(" compute_X_axis: freq offset = %f", delf)
      x = self.rel_freq_units(frame="FREQ-OBS", ref_freq=f_ref)-delf
    elif frame == "RADI-OBS":
      # radial velocity referred to the observer
      vobj = self.V_LSR(obstime, self.dss)
      x = self.rel_freq_units(frame=frame, ref_freq=f_ref)
    elif frame == "RADI-LSR":
      x = self.rel_freq_units(frame=frame, ref_freq=f_ref, 
                              v_frame=self.V_LSR(obstime, self.dss))
      vobj = v_ref
      self.logger.debug("compute_X_axis: vobj = %.2f", vobj)
    elif frame == "RADI-OBJ":
      # This is the object's rest frame
      self.logger.debug("compute_X_axis: vspline = %s", vspline)
      self.logger.debug("compute_X_axis: time = %s", obstime)
      if vspline and obstime:
        vobj = vspline(time.mktime(obstime.timetuple()))
        x = -(c/1000)*self.rel_freq_units(frame="DELF-OBS")/f_ref - vobj
      else:
        vobj = self.header['VELOCITY']
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref,
                                v_frame=self.V_LSR(obstime, self.dss) + vobj)
    else:
      self.logger.warning(" frame %s is not valid", frame)
      return
    self.header['velocity'] = vobj
    return x
                         
  def rel_freq_units(self, frame="FREQ-OBS", ref_freq=None, v_frame=None, 
                     row=0):
    """
    Returns the X-axis units for various reference frames
    
    The recognized frames are::
      CHAN-OBS - channel number
      FREQ-OBS - frequency (MHz) in the observatory frame
      DELF-OBS - freq. in the observatory frame relative to the reference freq.
      RADI-OBS - * radio Doppler shift (km/s) in the observatory frame
      OPTI-OBS - * optical Doppler shift (km/s) in the observatory frame
      RELA-OBS - * relativistic Doppler shift (km/s) in the observatory frame
      RADI-LSR - * radio Doppler shift (km/s) in the frame of the LSR
      OPTI-LSR - * optical Doppler shift (km/s) in the frame of the LSR
      RELA-LSR - * relativistic Doppler shift (km/s) in the frame of the LSR
    All the ones marked with * are valid SDFITS values for the VELDEF keyword.
    In addition, the prefix HELI- is allowed but not yet implemented here. The
    first three values are not for velocity scales and so are not valid VELDEF
    values.
   
    @param frame : reference frame for the reference frequency
    @type  frame : str
    
    @param v_frame : the velocity of the observer relative to the frame.
    @type  v_frame : float or None
    
    @param ref_freq : frequency used to compute frequency offsets
    @type  ref_freq : float
    
    @param scan : the scan to use to get observation data
    @type  scan : int
    """
    if ref_freq:
      f_ref = ref_freq
    else:
      f_ref = self.dataset['RESTFREQ'][row]/1e6 # reference frequency in MHz
    self.logger.info("rel_freq_units: reference frequency is %.3f", f_ref)
    rel_freqs = self.freqs(row)-f_ref          # channel frequencies relative
                                                # to the reference frequency
    v_ref = self.dataset['VELOCITY'][row]       # velocity of the object in LSR
    self.logger.info("rel_freq_units: reference velocity is %.3f", v_ref)
    self.frame = frame
    if frame == "CHAN-OBS":
      return range(len(rel_freqs))
    elif frame == "FREQ-OBS":
      return self.freqs()
    elif frame == "DELF-OBS":
      return rel_freqs
    elif frame == "RADI-OBS":  
      return doppler_radio(rel_freqs, f_ref)
    elif frame == "OPTI-OBS":
      return doppler_optical(rel_freqs, f_ref)
    elif frame == "RELA-OBS":
      return doppler_relat(rel_freqs, f_ref)
    elif frame == "RADI-LSR":
      return doppler_radio(rel_freqs, f_ref) - v_frame
    elif frame == "OPTI-LSR":
      return doppler_optical(rel_freqs, f_ref) - v_frame
    elif frame == "RELA-LSR":
      return doppler_relat(rel_freqs, f_ref) - v_frame
  
  def freqs(self, row=None):
    """
    Computes frequencies for spectra.
    
    Assumes scan 1 is typical of all scans in the dataset
    """
    if row == None:
      row = 0
    basefreq = self.dataset['OBSFREQ'][row]
    freqs = get_freq_array(self.dataset['BANDWIDT'][row], self.nchans)
    self.frame = "CHAN-OBS"
    if self.SB == "U":
      return (basefreq + freqs)/1e6
    elif self.SB == "L":
      return (basefreq - freqs)/1e6
    else:
      raise RuntimeError("freqs: sideband %s not not valid", self.SB)
   
  
  def get_spectra(self, row):
    """
    returns each record for each IF as a spectrum
    
    Basically, it just removes the degenerate coordinate axes
    
    Example::
      In [9]: row = 0; dd = examiner[0].dataset['SPECTRUM'][row,:,:,:,0,0,:]
      In [10]: dd.shape
      Out[10]: (2, 7, 2, 32768)
    """
    return self.dataset['SPECTRUM'][row,:,:,:,0,0,:]
    
  def diff_spectra(self, row):
    """
    to subtract beam 2 from beam 1 for each record and each pol
    """
    spectra = self.get_spectra(row)
    return spectra[0,:,:,:]-spectra[1,:,:,:]
  
  def ratio_spectra(self, row):
    """
    returns ratio of beam 1 divided by beam 2
    """
    spectra = self.get_spectra(row)
    return spectra[0,:,:,:]/spectra[1,:,:,:]
  
  def mean_spectra(self, row, type="ratio"):
    """
    returns the means of the ratio or difference spectra for each pol
    """
    if type == "ratio":
      D = self.ratio_spectra(row)
    else:
      D = self.diff_spectra(row)
    return D.mean(axis=0)
    
  def switched_spectra(self, on_row, off_row):
    """
    returns on-beam/off-beam - 1 averaged over two scans
    
    one SIG=True average divided by one SIG=False average
    """
    Don = self.mean_spectra(on_row)
    Doff = self.mean_spectra(off_row)
    return Don-Doff
  
  def load_ephemeris(self, sourcename):
    """
    Loads the ephemeris and returns a spline for the radial velocity with date
    """
    # the source for which the ephemeris was calculated
    if sourcename == "67P":
      # only works in /home/kuiper/Projects/Comets/67P for now
      import imp
      ephemeris_67P = imp.load_source('ephemeris_67P',
                                    '/usr/local/projects/67P/ephemeris_67P.py')
      from ephemeris_67P import get_ephemeris
      return get_ephemeris()
    else:
      self.logger.error("load_ephemeris: only have ephemeris for Comet 67P")
      raise RuntimeError("No ephemeris for %s" % opts.sourcename)
    
  def beam_switched_average(self, mode='LINE-PBSW'):
    """
    Returns pol average spectrum and total integration of given dataset(s)
    
    Also computes the average spectrum and total integration time for each pol
    and accumulates them in the SAOexaminer attributes avg_diff and tau.
    
    @param index : dataset index
    @type  index : int
    
    @param mode : observing mode (which aught to be in the header)
    @type  mode : str
    """
    first_scan = self.dataset.keys()[0]
    num_chan   = self.dataset[first_scan].shape[0]
    if mode == 'LINE-PBSW':
      spec, tau = self.dataset.norm_spectra()
      if self.avg_diff == {}:
        # cumulative average
        # xan't do it in __init__ because num_chan is not known
        self.avg_diff = {0: zeros(num_chan), 1: zeros(num_chan)}
        self.tau = {0: 0.0, 1: 0.0}      # cumulative integration time
      for pol in [0,1]:
        self.avg_diff[pol] += spec[pol]
        self.tau[pol] += tau[pol]
    else:
        self.logger.error("beam_switched_average: mode %s not yet implemented", mode)
    for pol in [0,1]:
      self.avg_diff[pol] /= len(self.indices)

