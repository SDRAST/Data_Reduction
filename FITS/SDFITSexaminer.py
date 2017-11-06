"""
Provides class to examine DSN SDFITS data files

At present there are two types of DSN SDFITS formats.

For files obtained from WVSRs, there is one feed, no time axis (i.e. no
records with a scan) but typically two subchannels with the band::
    subch 1 is in the even numbered rows
    subch 2 is in the odd numbered rows    
    the first two rows in a set of four are on source 
    the second two rows in a set of four are off source
    subch1 on  is [0::4] (subch 0, pos 0)
    subch2 on  is [1::4] (subch 1, pos 0)
    subch1 off is [2::4] (subch 0, pos 1)
    subch2 off is [3::4] (subch 1, pos 1)
A SPECTRUM cell data cube has four axes::
   (freq, RA, dec, Stokes)

For files obtained from the SAO spectrometers, there is only one CYCLE but
there are two feeds. A SPECTRUM cell has six axes::
   (freq, RA, dec, Stokes, time, beam)
"""
import astropy.io.fits as pyfits
import astropy.units as u
import datetime
import logging
import numpy
import os
import re
import sys
import time

from astropy.coordinates import FK4, FK5, SkyCoord
from copy import copy
from matplotlib.dates import date2num

from scipy import interpolate
from scipy.optimize import curve_fit

# Third party packages
from novas import compat as novas
from novas.compat import eph_manager

import support.lists

from Astronomy import c, MJD, v_sun
from Astronomy.redshift import doppler_radio, doppler_optical, doppler_relat
from Data_Reduction import clobber
from Data_Reduction.FITS.DSNFITS import get_indices, get_table_stats 
from DatesTimes import UnixTime_to_datetime
from MonitorControl.BackEnds import get_freq_array
from MonitorControl.Configurations.CDSCC.FO_patching import DistributionAssembly
from MonitorControl.Configurations.coordinates import DSS
from Radio_Astronomy import rms_noise
from support import mkdir_if_needed, nearest_index


jd_start, jd_end, number = eph_manager.ephem_open()
logger = logging.getLogger(__name__)

class DSNFITSexaminer(object):
  """
  Class for examining SAOdataset objects
  """
  def __init__(self, parent=None, FITSfile=None, hdulist=None):
    """
    Create an SAOexaminer object
    
    Attributes::
      file    - FITS file path and name
      header  - file primary HDU header
      hdulist - list of HDUs in file
      logger  - logging.Logger object
      tables  - private class Table (SINGLE DISH) objects
    """
    self.parent = parent
    self.logger = logging.getLogger(logger.name+'.DSNFITSexaminer')
    if FITSfile:
      self.hdulist = pyfits.open(FITSfile)
      # hdulist[0] has the file header
      self.tables = self.get_SDFITS_tables()
      self.file = FITSfile
      self.logger.info("__init__: %d tables found in %s",
                     len(self.tables), os.path.basename(FITSfile))
    elif hdulist:
      self.hdulist = hdulist
      self.tables = self.get_SDFITS_tables()
    else:
      self.logger.error("__init__: requires filename or HDU list")
    self.header = self.hdulist[0].header
       
  def get_SDFITS_tables(self):
    """
    Finds the SINGLE DISH extensions in a FITS file
    """
    tables = {}
    index = 0
    self.logger.debug("get_SDFITS_tables: checking %s", self.hdulist[1:])
    for extension in self.hdulist[1:]:
      if extension.header['extname'] == 'SINGLE DISH':
        self.logger.debug("get_SDFITS_tables: found %s", extension)
        tables[index] = DSNFITSexaminer.Table(self, extension)
        self.logger.debug("get_SDFITS_tables: created table %s", tables[index])
        # the following is a HACK until all the SDFITS files are fixed
        try:
          tables[index].SB = numpy.unique(tables[index].data['SIDEBAND'])
        except KeyError:
          if tables[index].data['CDELT1'][0] > 0:
            tables[index].SB = 'U'
          else:
            tables[index].SB = 'L'
        index += 1
    return tables
  
  def get_sources(self):
    """
    """
    sources = []
    for key in self.tables.keys():
      for source in self.tables[key].sources:
        if source: # is not empty
          sources.append(source)
    return support.lists.unique(sources)
  
  def save(self, filename):
    """
    save examiner in session directory
    
    If the examiner has a parent, the file is saved in the session working
    directory. Else, if the filename includes a path, that path is used.
    Otherwise, it is saved in the user's home directory.
    """
    if self.parent:
      savepath = self.parent.projworkpath
    elif os.path.dirname(filename):
      savepath = os.path.dirname(filename)
    else:
      savepath = os.environ['HOME']
    if filename[-5:] == ".fits":
      savename = savepath + filename
    else:
      savename = savepath + filename+".fits"
    self.logger.info('consolidate: writing %s', savename)
     
  class Table(pyfits.BinTableHDU):
    """
    Public Attributes::
      BE          - Backend subclass associated with this table
      dss         - DSN antenna
      cycles      - list of unique cycle numbers
      cycle_keys  - 
      data        - contents of FITS binary table
      frame       - reference frame for X-axis
      header      - FITS table header
      logger      - logging.Logger object
      num_indices - number of indices in the SPECTRUM cube
      obsmodes    - 
      obs_freqs   - non-redundant list of observing frequencies in table
      props       - table properties
      rows        - rows with valid data
      row_keys    - ordered list of row numbers
      scan_keys   - ordered list of scan numbers
      sources     - non-redundant list of sources in table
    
    Methods::
      align_spectrum
      freqs
      load_ephemeris        - get the ephemeris for the designated source
    These are in addition to the attributes of a pyfits.hdu.table.BinTableHDU
    """
    def __init__(self, parent, extension):
      """
      Create a FITS table object
    
      @param extension : SDFITS file BINTABLE extension
      @type  extension : pyfits HDU

      """
      pyfits.BinTableHDU.__init__(self, data=extension.data,
                                        header=extension.header)
      self.logger = logging.getLogger(parent.logger.name+".Table")
      self.logger.debug("__init__: created %s", self)
      # add some attributes for convenience
      d = UnixTime_to_datetime(self.get_first_value('UNIXtime',0)).timetuple()
      self.year = d.tm_year
      self.DOY = d.tm_yday
      self.datestr = "%4d/%03d" % (self.year,self.DOY)
      self.timestr = "%02d:%02d" % (d.tm_hour, d.tm_min)
      self.dss = int(self.header['TELESCOP'].split('-')[1])
      self.logger.debug("__init__: station is DSS-%d", self.dss)
      self.scan_keys, self.cycle_keys, self.row_keys, self.obsmodes \
                                                        = get_table_stats(self)
      self.rows = self.data['CYCLE'].nonzero()[0] # valid data
      self.props, self.num_indices = self.properties()
      datashape = extension.data['SPECTRUM'].shape
      self.BE = self.header['BACKEND']
      if self.BE == 'SAO spectrometer':
        # hack until SAO2SDFITS.py is fixed
        self.cycles = [1]
      else:
        self.cycles = list(numpy.unique(self.data['CYCLE'][self.rows]))
      self.sources = numpy.unique(self.data['OBJECT'][self.rows])
      self.sources.sort()

    def get_rows(self, keyword, value):
      """
      """
      return numpy.where(self.data[keyword] == value)[0]
    
    def report_table(self):
      """
      """
      text = "%3d rows with valid data\n" % len(self.rows)
      text += "cycles:  %s\n" % self.cycles
      text += "frequencies: %s\n" % self.get_obs_freqs()
      text += "sources: %s\n" % self.sources
      text += "table properties:  %s" % self.props
      self.logger.info("report_table:\n%s", text)
      return text
        
    def get_obs_freqs(self):
      """
      gets the values for the frequency or velocity channels
      
      The right way to do this is to look at CRTYP1 which can have the 
      following values (though not from WBDC2 which always has FREQ-OBS)::
        axis types
          'FREQ' - frequency (Hz)
          'VELO' - velocity (m/s) (radio convention, unless overridden by use 
                   of the VELDEF SHARED keyword)
          'FELO' - regularly gridded in frequency but expressed as velocity in
                   the optical convention (m/s)
        reference frames
          '-LSR' : Local Standard of Rest
          '-HEL' : heliocentric (barycentric)
          '-OBS' : the frame of rest of the observer/telescope (topocentric)
          'LSRK' : LSR as a kinematical definition
          '-GEO' : Geocentric
          'REST' : rest frequency
          '-GAL' : Galactocentric
      
      Then compute the frequency for each channel as::
        f(pixel) = CRVAL1 + CDELT1 * (pixel - CRPIX1)
        
      Reference::
        https://casa.nrao.edu/aips2_docs/notes/236/node14.html
      """
      self.obs_freqs = {}
      self.logger.debug("get_obs_freqs: checking cycles %s", self.cycles)
      self.logger.debug("get_obs_freqs: OBSFREQ shape is %s",
                        self.data['OBSFREQ'].shape)
      self.logger.debug("get_obs_freqs: 'rows' shape is %s", self.rows.shape)
      for cycle in self.cycles:
        self.logger.debug("get_obs_freqs: getting freqs for cycle %d", cycle)
        self.obs_freqs[cycle] = self.data['OBSFREQ'][self.rows]\
                                             [numpy.equal(self.rows, cycle)][0]
      return self.obs_freqs

    def get_index_keys(self):
      """
      get the row and construct a tuple for a data array index
  
      The maximum is (beam,record,IF,0,0) which will return a spectrum or whatever
      is in the first column of the multi-dimensional array.
      """
      print "These are the indices to provide to retrieve a spectrum:",
      if self.num_indices < 3:
        raise RuntimeWarning("minimum data array index is 'row, RA, dec'")
        return None
      elif self.num_indices == 3:
        keys = "row", "dec", "RA"
      elif self.num_indices == 4:
        keys = "row", "pol", "dec", "RA"
      elif self.num_indices == 5:
        keys = "row", "record", "pol", "dec", "RA"
      elif self.num_indices == 6:
        keys = "row", "beam", "record", "pol", "dec", "RA"
      else:
        raise RuntimeWarning("cannot have more that 6 axes in SPECTRUM array")
        return None
      return keys
    
    def get_indices(self, scan=1, cycle=1, pol=1, beam=1, record=1,
                    trimmed=False):
      """
      returns indices for getting one spectrum from SPECTRUM column
      
      @param scan : SCAN number
      @type  scan : int
  
      @param cycle : CYLE number
      @type  cycle : int
  
      @param beam : BEAM axis value
      @type  beam : int
  
      @param IF : IF number or STOKES axis value
      @type  IF : int
  
      @param record : 1-based TIME axis index (FITS/FORTRAN convention) 
      @type  record : int
  
      @param trimmed : return tuple with 'RA' and 'dec' indices removed (always 0)
      @type  trimmed : bool
      """
      scan_idx = self.scan_keys.index(scan)
      cycle_idx = self.cycle_keys.index(cycle)
      beam_idx = beam-1
      IF_idx = pol-1
      record -= 1
      return get_indices(self.num_indices, self.props, scan_idx=scan_idx, 
                         cycle_idx=cycle_idx, beam_idx=beam_idx, IF_idx=IF_idx,
                         record=record, trimmed=trimmed)
  
    def freqs(self, row=0):
      """
      Computes frequencies for spectra in MHz
    
      Assumes scan 1 is typical of all scans in the dataset
      """
      freqs = self.data['CRVAL1'][row] + self.data['CDELT1'][row] * \
               (numpy.arange(self.props['num chans'])-self.data['CRPIX1'][row])
      return freqs/1e6
                         
    def rel_freq_units(self,frame="FREQ-OBS", ref_freq=None, v_frame=0, row=0):
      """
      Returns the X-axis units for various reference frames
    
      The recognized frames are::
        CHAN-OBS - channel number
        FREQ-OBS - frequency (MHz) in the observatory frame
        DELF-OBS - freq. in the observatory frame relative to the ref. freq.
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
        f_ref = self.data['RESTFREQ'][row]/1e6 # ref freq in MHz
      rel_freqs = self.freqs(row)-f_ref          # channel frequencies relative
                                                 # to the reference frequency
      v_ref = self.data['VELOCITY'][row]         # vel. of object in LSR
      self.logger.debug("rel_freq_units: reference frequency is %.3f", f_ref)
      self.logger.debug("rel_freq_units: reference velocity is %.3f", v_ref)
      self.logger.debug("rel_freq_units: frame velocity is %.3f", v_frame)
      self.frame = frame
      if frame == "CHAN-OBS": # channel number
        return range(len(rel_freqs))
      elif frame == "FREQ-OBS": # receiver frequency in the observer's frame
        return self.freqs(row=row)
      elif frame == "DELF-OBS": # frequency relative to the reference pixel
        return rel_freqs
      elif frame == "RADI-OBS": # velocity relative to the reference pixel
        return doppler_radio(rel_freqs, f_ref)
      elif frame == "OPTI-OBS": # same but with optical (wavelength) equation
        return doppler_optical(rel_freqs, f_ref)
      elif frame == "RELA-OBS": # same but with full relativity 
        return doppler_relat(rel_freqs, f_ref)
      elif frame == "RADI-LSR": # velocity w.r.t. LSR
        return doppler_radio(rel_freqs, f_ref) - v_frame
      elif frame == "OPTI-LSR": # same but with optical equation
        return doppler_optical(rel_freqs, f_ref) - v_frame
      elif frame == "RELA-LSR": # same but with full 
        return doppler_relat(rel_freqs, f_ref) - v_frame
  
    def compute_X_axis(self, row, frame='RADI-OBS', ref_freq=None,
                       vspline=None, obstime=None):
      """
      Computes the appropriate X-axis for a spectrum.
     
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
      @type  obstime : UNIX seconds
    
      @return: numpy.ndarray
      """
      if not obstime:
        # get time from the first row
        obstime = self.get_first_value('UNIXtime', row)
      if ref_freq:
        f_ref = ref_freq
      else:
        f_ref = self.data['RESTFREQ'][0]/1e6 # MHz
      v_ref = self.data['VELOCITY'][0]
      if frame:
        self.frame = frame
      self.logger.debug("compute_X_axis: requested frame is %s", frame)
      self.logger.debug("compute_X_axis: reference frequency is %10.3f", f_ref)
      if frame == "CHAN-OBS" or frame == "FREQ-OBS" or frame == "RELA-OBS":
        # all these are in the observer's rest frame
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref)
        vobj = None
      elif frame == "FREQ-LSR":
        vobj = self.V_LSR(row) # km/s
        C = c/1000 # km/s
        delf = -self.V_LSR(row)*f_ref/C
        self.logger.debug(" compute_X_axis: freq offset = %f", delf)
        x = self.rel_freq_units(frame="FREQ-OBS", ref_freq=f_ref)-delf
      elif frame == "RADI-OBS":
        # radial velocity referred to the observer
        vobj = self.V_LSR(row)
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref)
      elif frame == "RADI-LSR":
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref, 
                                v_frame=self.V_LSR(row))
        vobj = v_ref
        self.logger.debug("compute_X_axis: vobj = %.2f", vobj)
      elif frame == "RADI-OBJ":
        # This is the object's rest frame
        self.logger.debug("compute_X_axis: vspline = %s", vspline)
        self.logger.debug("compute_X_axis: time = %s", obstime)
        if vspline and obstime:
          # convert obstime to UNIX time
          vobj = vspline(obstime)
          x = -(c/1000)*self.rel_freq_units(frame="DELF-OBS")/f_ref - vobj
        else:
          vobj = self.data['VELOCITY'][row]
          x = self.rel_freq_units(frame=frame, ref_freq=f_ref,
                                  v_frame=self.V_LSR(row) + vobj)
      else:
        self.logger.warning(" frame %s is not valid", frame)
        return
      return x
  
    def V_LSR(self, row=0):
      """
      Computes the velocity of the local standard of rest w.r.t. the observer
    
      @param dt : date/time of the observation
      @type  dt : datetime object
    
      @param dss : DSN station
      @type  dss : int
    
      @return: float (km/s)
      """
      try:
        if self.data['EQUINOX'][row] == 1950:
          position = SkyCoord(self.data['CRVAL2'][row],
                              self.data['CRVAL3'][row],
                              frame=FK4, unit=(u.hourangle, u.deg))
      except KeyError:
        # assume J2000
        pass
      else:
        position = SkyCoord(self.data['CRVAL2'][row],
                            self.data['CRVAL3'][row],
                            frame=FK5, unit=(u.hourangle, u.deg))
      ra = position.ra.hour
      dec = position.dec.deg
      self.logger.debug("V_LSR: ra = %f, dec = %f", ra, dec)
      cat_entry = novas.make_cat_entry(self.data["OBJECT"][row],
                                       "", 0, ra, dec, 0, 0, 0, 0)
      source = novas.make_object(2, 0, self.data["OBJECT"][0], cat_entry)
      longdeg = self.header['SITELONG']
      self.logger.debug("V_LSR: longitude in degrees W = %f", longdeg)
      if longdeg > 180:
        longdeg -= 360
      self.logger.debug("V_LSR: longitude in degrees E = %f", longdeg)
    
      DSS43 = novas.make_observer_on_surface(self.header['SITELAT'], longdeg,
                                             self.header['SITEELEV'], 0, 0)
      dt = UnixTime_to_datetime(self.get_first_value('UNIXtime', row))
      self.logger.debug("V_LSR: computing for %s", dt.ctime())
      jd = novas.julian_date(dt.year,dt.month,dt.day,dt.hour+dt.minute/60.)
      mjd = MJD(dt.year,dt.month,dt.day)
      earth = novas.make_object(0, 3, 'Earth', None)
      urthpos,urthvel = novas.ephemeris((jd,0), earth, origin=0)
      (obspos,obsvel) = novas.geo_posvel(jd,0,DSS43,0)
      self.logger.debug("V_LSR: Earth velocity = %s", urthvel)
      self.logger.debug("V_LSR: observer velocity = %s", obsvel)
      totvel = tuple(numpy.array(urthvel)+numpy.array(obsvel))
      self.logger.debug("V_LSR: total velocity = %s", totvel)
      (srcpos,srcvel) = novas.starvectors(cat_entry)
      self.logger.debug("V_LSR: source velocity = %s", srcvel)
      V = novas.rad_vel(source, srcpos, srcvel, totvel,0,0,0)
      self.logger.debug("V_LSR: velocity of LSR = %.2f km/s", V)
      return V+v_sun(mjd, position.ra.hour, position.dec.deg)
        
    def make_directory(self, dest=sys.stdout):
      """
      """
      labels = "Row Scan ch      Source       Sig Freq      intg"
      flines = "--- ---- -- ---------------- ---- --------- ----"
      lbform = "%3d  %3d %2d %16s %5s %9.3f %4d"
      print >> dest, labels
      print >> dest, flines
      for row in self.rows:
        print >> dest, lbform % (row, self.data['SCAN'][row],
                                      self.data['CYCLE'][row],
                                      self.data['OBJECT'][row],
                                      self.data['SIG'][row],
                                      self.data['OBSFREQ'][row]/1e6,
                                      self.data['EPOSURE'][row])        

    def get_good_rows(self):
      """
      eliminates rows which have bad data, like 0 for CYCLE, nan, etc
      
      Note that this provides a method for flagging bad data by setting the
      CYCLE value of a row to 0.
      
      The good data flags are True if there are any good data and False if
      there are no good data for that column.
      
      Good data values are returned in the dict 'good_data'
      """
      scan_keys, cycle_keys, row_keys, obsmodes = get_table_stats(self)
      if len(obsmodes) > 1:
        raise RuntimeErrror(
                   "get_good_rows: multiple observing modes not yet supported")
      num_scans = len(scan_keys)
      num_cycles = len(cycle_keys)
      num_rows = len(row_keys)
      props, num_indices = self.properties() # table properties
      # check the data:
      #    elevations
      elevation = self.data['ELEVATIO']
      mask = ~(numpy.isnan(elevation) | numpy.equal(elevation, 0))
      indices = numpy.where(mask)
      good_el_data = mask.any()
      if not good_el_data:
        self.logger.warning("get_good_rows: elevation data is bad")
      else:
        self.logger.debug("get_good_rows: elevations shape: %s",
                          elevation.shape)
      #    system temperature
      if numpy.equal(self.data['TSYS'], 0.).all():
        good_tsys_data = False
        self.logger.warning("get_good_rows: Tsys data is bad")
      else:
        good_tsys_data = True
        self.logger.debug("get_good_rows: system temperatures shape: %s",
                       self.data['TSYS'].shape)
      #    ambient temperature
      tamb = self.data['TAMBIENT']
      mask = ~(numpy.isnan(tamb) | numpy.equal(tamb, 0))
      good_tamb_data = mask.any()
      if not good_tamb_data:
        self.logger.warning("get_good_rows: ambient temperature data is bad")
      #    pressure
      pres = self.data['PRESSURE']
      mask = ~(numpy.isnan(pres) | numpy.equal(pres, 0))
      good_pres_data = mask.any()
      if not good_pres_data:
        self.logger.warning("get_good_rows: pressure data is bad")
      #    humidity
      humi = self.data['HUMIDITY']
      mask = ~(numpy.isnan(humi) | numpy.equal(humi, 0))
      good_humi_data = mask.any()
      if not good_humi_data:
        self.logger.warning("get_good_rows: humidity data is bad")
      #    windspeed
      wspe = self.data['WINDSPEE']
      mask = ~(numpy.isnan(wspe) | numpy.equal(wspe, 0))
      good_wspe_data = mask.any()
      if not good_wspe_data:
        self.logger.warning("get_good_rows: windspeed data is bad")
      #    winddir
      wdir = self.data['WINDDIRE']
      mask = ~(numpy.isnan(wdir) | numpy.equal(wdir, 0))
      good_wdir_data = mask.any()
      if not good_wdir_data:
        self.logger.warning("get_good_rows: wind direction data is bad")
        
      # Initialize for extracting data from simple (single value) columns
      good_data = {}
      good_data['mpltime'] = []
      if good_el_data:
        good_data['elev'] = []
      if good_tamb_data:
        good_data['Tambient'] = []
      if good_pres_data:
        good_data['pressure'] = []
      if good_humi_data:
        good_data['humidity'] = []
      if good_wspe_data:
        good_data['windspeed'] = []
      if good_wdir_data:
        good_data['winddirec'] = []
      if good_tsys_data:
        good_data['TSYS'] = {}
      # now make the multi-level dict
      for cycle in self.cycle_keys: # this loops over subchannels
        cycle_idx = cycle_keys.index(cycle)
        if good_data['TSYS'].has_key(cycle_idx):
          pass
        else:
          good_data['TSYS'][cycle_idx] = {}
        for beam_idx in range(self.props["num beams"]):
          if good_data['TSYS'][cycle_idx].has_key(beam_idx):
            pass
          else:
            good_data['TSYS'][cycle_idx][beam_idx] = {}
          for pol_idx in range(self.props["num IFs"]):
            if good_data['TSYS'][cycle_idx][beam_idx].has_key(pol_idx):
              pass
            else:
              good_data['TSYS'][cycle_idx][beam_idx][pol_idx] = [] 
      # process all scans
      for row in row_keys:
        # these are simple columns with multiple dimensions
        midnight_unixtime = time.mktime(time.strptime(
                                       self.data['DATE-OBS'][row], "%Y/%m/%d"))
        scan = self.data['SCAN'][row]
        cycle = self.data['CYCLE'][row]
        cycle_idx = cycle - 1
        if self.props['time axis'] == True:
          nrecs = self.props['num records'][cycle]
          for rec in range(nrecs):
            first_time = self.data['CRVAL5'][row]
            rectime = first_time + rec*self.data['CDELT5'][row] # numpy.array
            unixtime = midnight_unixtime + rectime
            datime = datetime.datetime.fromtimestamp(unixtime) # datetime
            good_data['mpltime'].append(date2num(datime))
            good_data['elev'].append(self.data['ELEVATIO'][row,0,rec,0,0,0,0])
            for beam_idx in range(self.props["num beams"]):
              beam = beam_idx+1
              for pol_idx in range(self.props["num IFs"]):
                pol = pol_idx+1
                indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                           beam=beam, record=rec)
                good_data['TSYS'][cycle_idx][beam_idx][pol_idx].append(
                                                    self.data['TSYS'][indices])
            good_data['Tambient'].append(self.data['TAMBIENT'][indices])
            good_data['pressure'].append(self.data['PRESSURE'][indices])
            good_data['humidity'].append(self.data['HUMIDITY'][indices])
            if good_wspe_data:
              good_data['windspeed'].append(self.data['WINDSPEE'][indices])
            if good_wdir_data:
              good_data['winddirec'].append(self.data['WINDDIRE'][indices])
        else:
          unixtime = self.data['UNIXtime'][row]
          datime = datetime.datetime.fromtimestamp(unixtime) # datetime object
          good_data['mpltime'].append(date2num(datime))
          good_data['elev'].append(self.data['ELEVATIO'][row])
          for beam_idx in range(self.props["num beams"]):
            beam = beam_idx+1
            for pol_idx in range(self.props["num IFs"]):
              pol = pol_idx+1
              self.logger.debug(
                           "get_good_rows: scan=%d, cycle=%d, beam=%d, pol=%d",
                           scan, cycle, beam, pol)
              indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                         beam=beam)
              self.logger.debug("get_good_rows: indices are %s", indices)
              good_data['TSYS'][cycle_idx][beam_idx][pol_idx].append(
                                                    self.data['TSYS'][indices])
          good_data['Tambient'].append(self.data['TAMBIENT'][row])
          good_data['pressure'].append(self.data['PRESSURE'][row])
          good_data['humidity'].append(self.data['HUMIDITY'][row])
          if good_wspe_data:
            good_data['windspeed'].append(self.data['WINDSPEE'][row])
          if good_wdir_data:
            good_data['winddirec'].append(self.data['WINDDIRE'][row])
      return good_data

    def get_first_value(self, column, row):
      """
      Get the first value in a vector cell.
  
      When there is a time axis, some column cells are vectors along the time axis
      of the data cube.  Often, just the first value is needed corresponding to the
      start of the scan.  This relieves the user of having to known how to access
      that datum.
      """
      try:
        cellshape = self.data[column][row].shape
      except AttributeError:
        return self.data[column][row]
      else:
        if cellshape == ():
          return self.data[column][row]
        else:
          idx = len(cellshape)*[[0]]
          return self.data[column][row][idx][0]
      
    def prepare_summary_arrays(self, num_chans):
      """
      Initiate dict of empty spectra indexed by sub-channel, beam and polariz'n
      
      Creates an empty spectrum for every scan, subchannel, beam and pol
      combination. This can be the start for averaging over all records for
      each combination. The averaged spectra can be used to produce a summary 
      plot of all the scans for each combination.
      """        
      spectra = {}
      for scan in self.scan_keys:
        scan_idx = self.scan_keys.index(scan) # top level dict for scans
        spectra[scan_idx] = {}
        for subch in self.cycles:
          subch_idx = subch-1
          spectra[scan_idx][subch_idx] = {} # 2e level dict is for subchannel
          for beam_idx in range(self.props["num beams"]):
            beam = beam_idx+1
            spectra[scan_idx][subch_idx][beam_idx] = {} # 4e levl dict for beam
            for IF_idx in range(self.props["num IFs"]): 
              pol = IF_idx+1
              spectra[scan_idx][subch_idx][beam_idx][IF_idx] = \
                                                      numpy.zeros((num_chans,))
              self.logger.debug(
               "prepare_summary_arrays: for scan %d subch %d, beam %d, pol %d",
               scan, subch, beam, pol)
              self.logger.debug("prepare_summary_arrays: spectrum shape is %s",
                          spectra[scan_idx][subch_idx][beam_idx][IF_idx].shape)
      return spectra
      
    def prepare_summary_images(self, num_chans):
      """
      Initiate dict of image arrays for SPECTRUM data
      
      The dicts are indexed by sub-channel, beam and polarization. Each one has
      a single row of zeros.  Each record of each scan for each combination can
      then be appended to form a 2D spectrum-time plot (dynamic spectrum).
      """
      images = {}
      for subch in self.cycles:
        subch_idx = subch-1
        images[subch_idx] = {}
        for beam_idx in range(self.props["num beams"]):
          beam = beam_idx+1
          images[subch_idx][beam_idx] = {}
          for IF_idx in range(self.props["num IFs"]): 
            pol = IF_idx+1
            images[subch_idx][beam_idx][IF_idx] = \
                               numpy.zeros((num_chans, 1))
            self.logger.debug(
                       "prepare_summary_images: for subch %d, beam %d, pol %d",
                       subch, beam, pol)
            self.logger.debug("prepare_summary_images: image shape is %s",
                              images[subch_idx][beam_idx][IF_idx].shape)
      return images

    def get_data_stats(self):
      """
      get the minimum, maximum, mean and std of values in the spectrum array
      """
      good_data = numpy.nonzero(self.data['CYCLE'])
      all_samples = self.data['SPECTRUM'][good_data]
      return all_samples.min(),  all_samples.max(), \
             all_samples.mean(), all_samples.std()

    def average_records(self, datasource, scan, subch, beam, pol):
      """
      average records in a scan
      
      To fill up the images arrays, the very first spectrum 
      
      @param images : spectra from all records in a 2D array
      @param datasource : SPECTRUM, or ifpower if Stokes
      @param scan : SCAN
      @param subch : CYCLE
      @param beam : 1-based number
      @param pol : 1-based number, i.e., 1 or 2
      """
      self.logger.debug(
                  "average_records: for scan %d subch%d beam%1d pol%d from %s",
                  scan, subch, beam, pol, datasource)
       
      IF_idx = pol-1
      cycle_idx = subch-1
      num_records = self.props["num records"][subch]
      first_spectrum = True # first record's spectrum
      for record in range(num_records):
        indices = self.get_indices(scan=scan, cycle=subch, beam=beam, IF=pol,
                                   record=record)
        spectrum = self.data[datasource][indices].reshape(
                                                     self.props['num chans'],1)
        # if image array is empty, initialize with first spectrum
        if first_spectrum:
          sum_spectrum = copy(spectrum)
          # we want time to be the row dimension (first index) and
          # frequency the column dimension (second index)
          image = spectrum.transpose()
          first_spectrum = False
        else:
          sum_spectrum += spectrum
          image = numpy.append(image, spectrum.transpose(), axis=0)
      # normalize
      sum_spectrum /= num_records
      return image, sum_spectrum[:,0]
   
    def get_spectra(self, rows):
      """
      returns each record for each IF as a spectrum
    
      Basically, it just removes the degenerate coordinate axes.  Use
      'get_index_keys' to get the order or the indices of the original array.
      
      Example::
        In [67]: tb0.data['SPECTRUM'].shape
        Out[67]: (12, 2, 22, 2, 1, 1, 32768)
        In [68]: specs = tb0.get_spectra(tb0.row_keys)
        In [69]: specs.shape
        Out[69]: (12, 2, 22, 2, 32768)
        
      @param rows : a list of row numbers
      @type  rows : list of int
      """
      nrows = len(rows)
      if (self.data['OBSMODE'][rows] == numpy.array(nrows*['LINEPBSW'])).all():
        # this requires a beam axis so...
        return self.data['SPECTRUM'][rows,:,:,:,0,0,:]
      elif (self.data['OBSMODE'][rows] == numpy.array(nrows*['LINEPSSW'])).all():
        if self.num_indices == 6:
          return self.data['SPECTRUM'][rows,:,:,:,0,0,:]
        elif self.num_indices == 5:
          # has a time axis
          return self.data['SPECTRUM'][rows,:,:,0,0,:]
        elif self.num_indices == 4:
          # has pol axis
          return self.data['SPECTRUM'][rows,:,0,0,:]
        else:
          # only one spectrum
          return self.data['SPECTRUM'][rows,0,0,:]
    
    def normalized_beam_diff(self, rows, remove_warning=False):
      """
      (on-source - off-source)/off_source for each record and each pol
      
      Example::
        In [70]: tb0.data['SPECTRUM'].shape
        Out[70]: (12, 2, 22, 2, 1, 1, 32768)
        In [71]: norms = tb0.normalized_beam_diff(tb0.row_keys)
        In [72]: norms.shape
        Out[72]: (12, 22, 2, 32768)
      This removes the beam axis.
      
      Note on EXPOSURE
      ================
      The SAO spectrometer nominally take 5 sec records, which in fact take
      longer. So fewer records are read out during a scan but we assume that
      the total integration time for the scan is still approximately the time
      that was intended. Then the effective integration time per record is
      5*number_of_intended_records/number of actual_records.
        
      @param rows : a list of row numbers
      @type  rows : list of int
      """
      if self.num_indices > 4: # 5 or more axes
        spectra = self.get_spectra(rows)
        # old TAMS spectra may have records with zeros at the end
        # any record may have a zero in the first 100 channels
        while (spectra[:,:,:,:,100:] == 0).any() == True:
          # remove last records to get rid of empty records
          nrecs = spectra.shape[2]
          spectra = spectra[:,:,:-1,:,:]
          if remove_warning:
            self.logger.warning("normalized_beam_diff: removed last records")
          # compensate integration time
          n_TAMS_recs = self.data['EXPOSURE']/5
          TAMS_intgr = 5*nrecs/n_TAMS_recs
          self.data['EXPOSURE'] -= TAMS_intgr
        # spectra indices: row, beam, record, pol, freq
        diff = spectra[:,0,:,:,:]-spectra[:,1,:,:,:] # beam 1 - beam 2
        # diff indices: row, record, pol, freq
        if self.data[rows]['SIG'][::2].all(): # source in beam 1
          diff[::2,:,:,:] /= spectra[::2,1,:,:,:]
        else:
          # source not in beam (i.e. in off beam)
          self.logger.error("normalized_beam_diff: %s has a non-SIG scan",
                            rows[::2])
          return None
        if self.data[rows]['SIG'][1::2].all() == False: # source in beam 2
          diff[1::2,:,:,:] /= -spectra[1::2,0,:,:,:]
        else:
          self.logger.error("normalized_beam_diff: %s has a non-REF scan",
                            rows[1::2])
          return None
        return diff
      else:
        raise RuntimeWarning("normalized_beam_diff: data cube has no beam axis")
        return None
    
    def BPSW_spectra(self, rows, Tsys=None, weighted=False):
      """
      returns record-by-record on-beam minus off-beam
      
      The optional Tsys array has, for each row pair, an array of off-source
      system temperatures to be used to normalize the beam 1 and beam 2
      normalized diffrences respectively.
      
      Example::
        In [70]: tb0.data['SPECTRUM'].shape
        Out[70]: (12, 2, 22, 2, 1, 1, 32768)
        In [74]: bpsw0 = tb0.BPSW_spectra(tb0.row_keys)
        In [75]: bpsw0.shape
        Out[75]: (6, 22, 2, 32768)
      The rows were reduced by pairs of one on-source and one off-source.
      
      For the first years of operation, the SAO configuration did not have
      beam 1 system temperatures for the 22H sub-band. Two power meters were
      assigned to the two pols of beam 2.  To construct a Tsys array which
      uses the beam 2 data for both beam difference spectra, do this::
        In [76]: tb0.data['TSYS'][tb0.row_keys].shape
        Out[76]: (12, 2, 22, 2, 1, 1, 1)
        In [77]: tsys0 = (tb0.data['TSYS'][::2,1,:,:,0,0] \
                        + tb0.data['TSYS'][1::2,1,:,:,0,0])/2
        In [78]: tsys0.shape
        Out[78]: (6, 22, 2, 1)
        In [79]: new = tsys0.reshape(6,1,22,2,1)
        In [80]: newnew = numpy.append(new,new,axis=1)
        In [81]: newnew.shape
        Out[81]: (6, 2, 22, 2, 1)
        
      @param rows : row with source alternating in beam 1 then beam 2
      @type  rows : list of int
      
      @param Tsys : an array of system temperatures with shape (6, 2, 22, 2, 1)
      @type  Tsys : nparray of float
      """
      on_rows = rows[::2]
      Don = self.normalized_beam_diff(rows)[::2]
      # SPECTRUM array may have had last records removed
      num_recs = Don.shape[1]
      if type(Tsys) == numpy.ndarray:
        Don *= Tsys[:,1,:num_recs,:,:]
      #self.logger.debug("BPSW_spectra: Don is now %s", Don)
      off_rows = numpy.array(on_rows)+1
      Doff = self.normalized_beam_diff(rows)[1::2]
      if type(Tsys) == numpy.ndarray:
        Doff *= Tsys[:,0,:num_recs,:,:]
      if weighted:
        # use off-source Tsys
        weight_on = 1./Tsys[:,1,:num_recs,:,:]**2
        weight_off = 1./Tsys[:,0,:num_recs,:,:]**2
        sum_weights = weight_on + weight_off
        weight_on /= sum_weights
        weight_off /= sum_weights
        return weight_on*Don + weight_off*Doff
      else:
        return (Don + Doff)/2
  
    def BPSW_average(self, rows, weighted=False, TAMS_hack=True):
      """
      Produce scaled, averaged spectra for each scan pair
      
      TAMS_hack uses beam 2 for both the SIG and REF Tsys because the beam 1
      data are not available.
      
      Returns a tuple with::
        * an array with spectrum for each pol for all scans (n_pairs,2,n_chans)
        * an array of averaged system temperatures   shape: (n_pairs,2,1)
        * an array of integration times              shape: (n_pairs,2,1)
      
      @param rows : rows from a table, starting with a SIG row
      @type  rows : numpy.ndarray
      
      @param weighted : spectra will be weighted with 1/Tsys^2 if True
      @type  weighted : bool
      
      @param TAMS_hack : only use beam 2 system temperatures
      @type  TAMS_hack : bool
      
      @return: (spectra, Tsys)
      """
      SIG_ref_beam = 1
      if TAMS_hack:
        REF_ref_beam = 1
      else:
        REF_ref_beam = 0
      # get beam 2 system temperatures
      #   rows must be alternating sig and ref
      sig_rows = rows[::2]
      ref_rows = rows[1::2]
      # indices: rows, beams, records, pols, dec, RA
      tsys_B2 = (self.data['TSYS'][sig_rows,SIG_ref_beam,:,:,0,0] + 
                 self.data['TSYS'][ref_rows,REF_ref_beam,:,:,0,0])/2
      # insert index for sig/ref beams
      npairs = tsys_B2.shape[0]
      nrecs = tsys_B2.shape[1]
      t = tsys_B2.reshape(npairs,1,nrecs,2,1)
      Tsys = numpy.append(t, t, axis=1)
      # one always weights by integration time
      integr = self.data['EXPOSURE'].reshape((npairs,2,1))
      # get the scaled difference spectra
      bpsw = self.BPSW_spectra(rows, Tsys=Tsys, weighted=weighted)
      return bpsw.mean(axis=1), Tsys.mean(axis=2).mean(axis=1), integr
    
    def scans_average(self, rows, weighted=True, TAMS_hack=True):
      """
      averages the scans weighted by int. time and 1/Tsys^2
      
      Returns a tuple with::
        * an array with spectrum for each pol             shape: (2, num_chans)
        * an array of averaged system temperatures        shape: (2,)
        * an array of averaged system temperatures        shape: (2,)
      """
      bpsw, Tsys, intgr = self.BPSW_average(rows, 
                                     weighted=weighted, TAMS_hack=TAMS_hack)
      weight = intgr/Tsys**2
      weight /= weight.sum(axis=0)
      wbpsw = weight*bpsw
      avewbpsw = wbpsw.sum(axis=0)
      # the integration time for each pol
      sumintgr = intgr.sum(axis=0).reshape((2,))
      # Tsys for each beam and pol
      aveTsys = (intgr*Tsys).sum(axis=0)/sumintgr
      # Tsys averaged over beams
      meanTsys = aveTsys.mean(axis=1)
      return avewbpsw, meanTsys, sumintgr
          
    def get_HP_PM_ID(self, row, pol, beam):
      """
      See if there is a power meter attached to a particular IF.
      
      This should probably be a DistributionAssembly method
      """
      da = DistributionAssembly()
      d = UnixTime_to_datetime(self.get_first_value('UNIXtime',0)).timetuple()
      da.get_sheet_by_date("%4d/%03d" % (d.tm_year, d.tm_yday))
      pms = da.get_signals('Power Meter')
      pm_keys = pms.keys()
      pm_keys.sort()
      for pm_key in pm_keys:
        self.logger.debug("get_HP_PM_ID: processing %s", pm_key)
        pm_index = pm_keys.index(pm_key)
        self.logger.debug("get_HP_PM_ID: processing index %d", pm_index)
        if   pms[pm_key]['Band'] == int(self.data['OBSFREQ'][row]/1e9) and \
           ((pms[pm_key]['IF'] == 'U' and self.data['SIDEBAND'][row] == +1) or \
            (pms[pm_key]['IF'] == 'L' and self.data['SIDEBAND'][row] == -1)) and \
           ((pms[pm_key]['Pol'] == 'E' and pol == 1) or \
            (pms[pm_key]['Pol'] == 'H' and pol == 2)) and \
             pms[pm_key]['Receiver'] == beam:
          return pm_index, pm_key
      return None
      
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
  
    def extract_window(self, rows, data=None, Tsys=None, intgr=None,
                       xlimits=(0,32767), frame="CHAN-OBS", source="67P"):
      """
      Extracts a subset of a SPECTRUM column.
      
      The datacube can be specified by a list of rows, in which case the
      spectra will be extracted from the SPECTRUM column, or else a data cube
      returned from 'normalized_beam_diff', 'BPSW_spectra', 'BPSW_average', or
      'scans_average'.  In the latter two cases, the averaged system 
      temperatures and integration times returned by these methods can also be
      provided.
    
      The X limits are in units appropriate to the frame. The corresponding
      channel numbers will be computed.
      
      If 'data' are provided then only the first row number is used to get the
      time of the observations.
      
      @param rows : required rows to select
      @type  rows : list or 1D nparray
      
      @param data : optional data cube
      @type  data : nparray
      
      @param xlimits : range of X-values to extract (min, max)
      @type  xlimits : tuple 
      """
      row = rows[0]
      # This computes values for the entire 32768 point spectrum
      if frame == "RADI-OBJ":
        # requires a spline computed from an ephemeris
        if re.match(source, self.data['OBJECT'][row]):
          vspline = self.load_ephemeris(source)
          x = self.compute_X_axis(row, "RADI-OBJ", vspline=vspline)
        else:
          self.logger.error("extract_window: %s does not match %s in row %d", 
                            source, self.data['OBJECT'][row], row)
          return None
      else:
        x = self.compute_X_axis(row, frame)
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
      if type(data) == numpy.ndarray:
        spectra = data
      else:
        # extract the rows
        spectra = self.data['SPECTRUM'][rows]
      # simplify the spectra
      #   in all cases the first ':' is for the sacn numbers
      if len(spectra.shape) == 7:
        # SAO-type original spectrum
        # indices: scan,beam,record,pol,dec,RA,freq
        s = spectra[:,:,:,0,0,ch1:ch2]
      elif len(spectra.shape) == 5:
        if spectra.shape[-3] == 1 and spectra.shape[-2] == 1:
          # dual pol WVSR-type spectrum
          # indices: scan,pol,dec,RA,freq
          s = spectra[:,:,0,0,ch1:ch2]
        else:
          self.logger.warning("extract_window: shape %s not recognized",
                              spectra.shape)
          s = data
      elif len(spectra.shape) == 4:
        # single pol WVSR-type original spectrum
        # indices: scan,pol,dec,RA,freq
        s = spectra[:,0,0,ch1:ch2]
      elif len(spectra.shape) == 3:
        # reduced dual pol WVSR spectrum
        # indices: scan,pol,freq
        s = spectra[:,:,ch1:ch2]
      elif len(spectra) == 2:
        # reduced single pol WVSR spectrum
        s = spectra
      else:
        raise RuntimeError("invalid SPECTRUM array shape: %s", spectra.shape)
      d = DSNFITSexaminer.Table.Data(self, x[ch1:ch2], s[:,ch1:ch2])
      d.frame = frame
      d.channels = (ch1,ch2)
      d.rows = list(rows)
      if type(Tsys) == numpy.ndarray:
        d.tsys = Tsys
      if type(intgr) == numpy.ndarray:
        d.intgr = intgr
      return d
      
    def reduce_line(self, rows=[], window=[-100,100], frame="RADI-OBJ",
                    source='67P'):
      """
      """
      if rows == []:
        rows = self.row_keys
      scanave, Tsys, intgr = self.scans_average(rows)
      subset = self.extract_window(rows, data=scanave, Tsys=Tsys,intgr=intgr,
                                   xlimits=window, frame=frame, source=source)
      if subset:
        x, y, rms = subset.remove_baseline()
        return x, y, rms, Tsys, intgr
      else:
        self.logger.error("reduce_line: no data window found")
        return None
                                      
    def fit_mean_power_to_airmass(self, Tvac_func, replace=False):
      """
      fits the mean power data vs airmass to the radiative transfer equation
    
      This assumes that every IF has a way of measuring power. The measured power
      is a single value along the first axis of the data array (or last index in
      a C/Python array).  If there are multiple records then they will be
      averaged.
    
      @param Tvac_func - a function for system temperature with no atmosph or CBR
      @type  Tvac_func - function(beam,pol)
      """
      def opacity_fitting(x, a, tau):
        """
        @param x : elevation in degrees
        @type  x : numpy.array of float
        
        @param a : intercept
        @type  a : float
        
        @param tau : slope
        @type  tau : slope
        """
        x_rad = numpy.deg2rad(x)   # radians
        x_sec = 1/numpy.sin(x_rad) # secant
        return a + tau*x_sec
    
      for beam_idx in range(self.props['num beams']):
        for IFidx in range(self.props['num IFs']):
          self.logger.debug('fit_mean_power_to_airmass: processing IF%d', IFidx+1)
          Tvac = Tvac_func(beam=beam_idx, pol=IFidx)
          pol = ['L','R'][IFidx]
          msg = "estimated %sCP zenith Tsys in vacuum= %6.2f" % (pol, Tvac)
          self.header.add_history(msg)
          # average records here if needed.
          subch_IDs = range(self.props['num cycles'])
          self.logger.debug("fit_mean_power_to_airmass: subchannels: %s",
                          subch_IDs)
          for subch in subch_IDs:
            subchannel = subch+1
            # Get the data for this subchannel.
            #   the following can be expressed as
            #     mean_power = self.data['avgpower'][subch::2,IFidx,0,0,0]
            #   or
            #     mean_power = self.data['avgpower'][:,IFidx,0,0,0][subch::2]
            #   or
            #     mean_power = self.data[subch::2]['avgpower'][:,IFidx,0,0,0]
            # assuming the first row is the same as the subchannel.  We are
            # using the first form.  In this form
            # WVSR data has indices (row,IF,dec,RA,chan)
            # SAO  data has indices (row,beam,time,IF,dec,RA,chan)
            if self.props['num beams'] > 1:
              mean_power = self.data['TSYS'][subch::2,beam_idx,:,IFidx,0,0,0]
              elevation  = self.data['ELEVATIO'][subch::2,0,:,0,0,0,0]
            elif self.props['time axis'] == False:
              if 'avgpower' in self.data.columns.names:
                mean_power = self.data['avgpower'][subch::2,IFidx,0,0,0]
              else:
                mean_power = self.data['TSYS'][subch::2,IFidx,0,0,0]
              elevation  = self.data['ELEVATIO'][subch::2]
            else:
              self.logger.error("fit_mean_power_to_airmass: unknown data shape")
              raise RuntimeError("unknown data shape")
            pw_shape = mean_power.shape
            el_shape = elevation.shape
            self.logger.debug("fit_mean_power_to_airmass: mean_power shape is %s",
                              pw_shape)
            self.logger.debug("fit_mean_power_to_airmass: elevation shape is %s",
                              el_shape)
            if pw_shape != el_shape:
              self.logger.error("fit_mean_power_to_airmass: shapes must be the same")
              raise RuntimeError("shapes must be the same")
            # get an elevation array with the 'nan' and zero values removed
            mask = ~(numpy.isnan(elevation) | numpy.equal(elevation, 0))
            self.logger.debug("fit_mean_power_to_airmass: mask shape is %s",
                              mask.shape)
            elv = elevation[mask]
            # remove the items with the same indices from mean_power
            pwr = mean_power[mask]
            #self.logger.debug('fit_mean_power_to_airmass: elevations: %s', elv)
            self.logger.debug('fit_mean_power_to_airmass: subch %d pwr shape: %s',
                          subchannel, pwr.shape)
            # fit the data
            popt, pcov = curve_fit(opacity_fitting, elv, pwr, p0=[0, 0])
            intercept, slope = popt[0], popt[1]
            self.logger.info(
             "fit_mean_power_to_airmass: B%dP%d sch%d intercept, slope: %f, %f",
                         beam_idx+1, IFidx+1, subchannel, intercept, slope)
            msg = \
            "IF%d, subch%d gain=%9.3e counts, gain_slope=%9.3e counts/airmass"\
                % (IFidx+1, subchannel, intercept, slope)
            self.header.add_history(msg)
            if replace:
              gain = Tvac/intercept
              K_per_am = gain*slope
              self.logger.debug(
           "fit_mean_power_to_airmass: convert power to Tsys for subch%s %sCP",
                          subch+1, pol)
              new_indices = numpy.where(self.data['CYCLE'] == subchannel)[0]
              self.logger.debug(
                          "fit_mean_power_to_airmass: table rows for Tsys: %s",
                          new_indices)
              self.logger.debug(
                          "fit_mean_power_to_airmass: destination shape is %s",
                          self.data['TSYS'][new_indices,IFidx,0,0,0].shape)
              self.data['TSYS'][new_indices,IFidx,0,0,0] = gain * pwr
            #else:
            #  self.logger.warning(
            #            "fit_mean_power_to_airmass: failed; Tsys not computed")

    def properties(self):
      """
      get table properties
  
      properties::
        num chans      - number of spectrometer channels in diagnostic spectra
        num IFs        - at most two, one per pol
        full Stokes    - four Stkes parameters instead of an IF for each pol
        time axis      - True if scan is divided into records
        num beams      - number of beams
        num records    - if 'time axis', number of records in each scan
      Notes
      =====
       * In a DSN SDFITS file there is only one backend
       * If there is a time axis, the number of records per scan may differ
       * A subchannel is a piece of band for a spectrometer with coarse and fine
       channels
       * cycles are used for dfferent subchannels and for position switching
      """
      props = {}
      # most parameters must be the same for all rows in a session
      spectrumshape = self.data['SPECTRUM'][0].shape # for common dims
      props["num cycles"] = \
               len(numpy.unique(self.data['CYCLE'][self.data['CYCLE'] != 0]))
      if len(spectrumshape) == 3: # (dec, RA, freq)
        # bare minimum SPECTRUM dimensions
        props["num chans"] = int(spectrumshape[-1])
        props["num IFs"] = 1        # no polarization axis
        props["num beams"] = 1      # no beam axis
      elif len(spectrumshape) >= 4: # (at least pol, dec, RA, freq)
        # one beam with polarization at least
        props["num chans"] = int(spectrumshape[-1])
        props["num IFs"] = 2
        if spectrumshape[-4] == 4:
          # has STOKES dimension
          props["num beams"] = 1     # no beam axis
          props["full Stokes"] = True
          if 'IFSPECTR' in self.data.columns.names:
            props["num IFs"] = 2
            IFspecshape = self.data['IFSPECTR'][0].shape
            props["num IFspec chans"] = int(IFspecshape[-1])
          else: # must be 3 or less
            props["num _IFs"] = 1
            self.logger.warning(
                    "properties: no IF data; will use Stokes I for monitoring")
            props["num IFspec chans"] = props["num chans"]
        elif spectrumshape[-4] == 2:
          # has straight pols (L,R or H,V)
          props["num IFs"] = 2
          props["full Stokes"] = False
          props["num beams"] = 1     # no beam axis
        if len(spectrumshape) >= 5: # (time, pol, dec, RA, freq)
          props["time axis"] = True
          if len(spectrumshape) == 5:
            props["num beams"] = 1     # no beam axis
        else:
          props["time axis"] = False
        if len(spectrumshape) == 6: # (beam, time, pol, dec, RA, freq)
          props["num beams"] = int(spectrumshape[0])
        else:
          props["num beams"] = 1
        # time axis length may vary due to corrupted records
        if len(spectrumshape) >= 5: # (time, pol, dec, RA, freq)
          props["num records"] = {}
          cycle_indices = range(props["num cycles"])
          for cycle_idx in cycle_indices: # check each cycle
            cycle = self.data['CYCLE'][cycle_idx]
            spectrumshape = self.data['SPECTRUM'][cycle_idx].shape
            # just do the first scan
            props["num records"][cycle] = int(spectrumshape[1])
        else:
          props["num records"] = {}
          for cycle in self.cycle_keys:
            props["num records"][cycle] = 1
        self.logger.debug("properties:\n %s", props)
        return props, len(spectrumshape)

    def remove_tones(self, rows=None):
      """
      """
      if rows:
        pass
      else:
        rows = self.row_keys
      # the nearest tone below the center frequency
      for row in rows:
        freq_offset = self.data['OBSFREQ'][row] - \
                                            round(self.data['OBSFREQ'][row],-6)
        self.logger.debug("remove_tones: frequency offset = %d", freq_offset)
        num_chans = self.props['num chans']
        # the center of the band is the lowest channel of the upper half
        cntr_chan = num_chans/2
        # the nearest tone below the center channel
        chan_offset = round(num_chans*freq_offset/self.data['BANDWIDT'][row])
        self.logger.debug("remove_tones: channel offset = %d", chan_offset)
        number_of_tones = int(self.data['BANDWIDT'][row]/1e6)
        tone_spacing = self.props['num chans']/number_of_tones
        
        for tone in range(-number_of_tones/2, number_of_tones/2):
          if chan_offset > 0:
            # there are more tones above the tone nearest the center than below
            tone += 1
          tone_channel = cntr_chan - chan_offset + tone*tone_spacing
          self.logger.debug("remove_tones: tone %d channel = %d",
                            tone, tone_channel)
          if self.props['full Stokes']:
            pols = range(4)
          else:
            pols = self.props['num IFs']
          for pol in pols:
            self.data['SPECTRUM'][row][pol,0,0] = clobber(
                             self.data['SPECTRUM'][row][pol,0,0], tone_channel)
      
    class Data(object):
      """
      A subset of data extracted from an average difference spectrum
  
      Attributes::
        channels - channel range selected from the original spectrum
        frame    - reference frame for the X-axis
        logger   - logging.Logger
        rows     - rows selected from parent
        x        - X-axis values
        y        - Multi-dimensional data array
      """
      def __init__(self, parent, x, y):
        """
        initiate a Data object
        
        @param x : X-values
        @type  x : list
        """
        self.logger = logging.getLogger(parent.logger.name+".Data")
        self.parent = parent
        if type(x) == list:
          self.x = numpy.array(x)
        elif type(x) == numpy.ndarray:
          self.x = x
        else:
          self.logger.error("__init__: type %s is not valid for X", type(x))
        self.y = y
        self.frame = None
        self.channels = None
        self.rows = None
      
      def __add__(self, other):
        """
        combine two data objects
        """
        if type(other) != DSNFITSexaminer.Table.Data:
          raise RuntimeError("__add__: other type wrong: %s", type(other))
        if other.frame != self.frame:
          raise RuntimeError("__add__: other must have frame type %s",
                             self.frame)
        new_y = numpy.append(self.y, other.y, axis=0)
        new_data = DSNFITSexaminer.Table.Data(self.parent, self.x, new_y)
        new_data.rows = self.rows + other.rows
        new_data.frame = self.frame
        new_data.channels = self.channels
        return new_data
      
      def remove_baseline(self, exclude=None):
        """
        """
        if exclude:
          minexcl, maxexcl = exclude
          x = numpy.append(self.x[:minexcl], self.x[maxexcl:])
          y = numpy.append(self.y[:minexcl], self.y[maxexcl:])
        else:
          x = self.x
          y = self.y
        if x[0] > x[-1]:
          reverse = True
          x = x[::-1]
          y = y[:,::-1]
        else:
          reverse = False
        new_y = numpy.empty_like(self.y)
        rms = []
        for pol in range(self.y.shape[0]):
          spline = interpolate.splrep(x, y[pol], k=5, s=10)
          residuals = y[pol] - interpolate.splev(x, spline, der=0)
          rms += [residuals.std()]
          # apply baseline correction
          if reverse:
            new_y[pol] = self.y[pol] - interpolate.splev(self.x[::-1], spline, der=0)[::-1]
          else:
            new_y[pol] = self.y[pol] - interpolate.splev(self.x, spline, der=0)
        return self.x, new_y, numpy.array(rms)
 
      
  # convert these to Table methods as needed.
    
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
  
def get_power_range(table, subch, beam, pol):
            """
            calculate limits for power levels
            
            the idea here is to avoid huge spikes compressing the spectra
            will need some tuning
            
            @param subch : sub-channel number
            @param beam : beam number (1-based)
            @param pol   : polarization number (not NRAO pol code)
            """
            subch_idx = subch-1
            beam_idx = beam-1
            IF_idx = pol-1
            # gets the indices for scan 1
            indices = table.get_indices(cycle=subch, beam=beam, IF=pol)
            logger.debug("get_power_range: indices: %s", indices)
            # gets spectra from all rows
            spectrum = table.data['SPECTRUM'][:,indices[1:]]
            logger.debug("get_power_range: spectrum shape is %s", spectrum.shape)
            mean_pwr = spectrum.mean()
            pwr_std = spectrum.std()
            max_pwr = spectrum.max()
            min_pwr = spectrum.min()
            if max_pwr > mean_pwr + 4*pwr_std:
              ymax = mean_pwr + pwr_std
            else:
              ymax = mean_pwr + 4*pwr_std
            if min_pwr < mean_pwr - 4*pwr_std:
              ymin = mean_pwr - pwr_std
            else:
              ymin = mean_pwr - 4*pwr_std
            return ymin, ymax

