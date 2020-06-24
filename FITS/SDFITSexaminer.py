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

Notes
=====
Other FITS Formats
------------------
These are not handled right now except to give some information::
In [1]: from Data_Reduction.FITS.SDFITSexaminer import DSNFITSexaminer
In [3]: mylogger.setLevel(logging.INFO)
In [4]: ex = DSNFITSexaminer(FITSfile="Carina_H92a_Tid_Gal.fits")
INFO:Data_Reduction.FITS.SDFITSexaminer.DSNFITSexaminer:__init__:
                                     0 tables found in Carina_H92a_Tid_Gal.fits
WARNING:Data_Reduction.FITS.SDFITSexaminer.DSNFITSexaminer:__init__:
                                                      this has no SDFITS tables
INFO:Data_Reduction.FITS.SDFITSexaminer.DSNFITSexaminer:__init__:
                                                     this is a simple FITS file
INFO:Data_Reduction.FITS.SDFITSexaminer.DSNFITSexaminer:__init__:
                                         data matrix shape is (1, 1800, 21, 31)
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
#from matplotlib.dates import date2num
 
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.stats import linregress

# Third party packages
from novas import compat as novas
from novas.compat import eph_manager

logger = logging.getLogger(__name__)

import support.lists

from Astronomy import c, v_sun
from Astronomy.DSN_coordinates import DSS
from Astronomy.redshift import doppler_radio, doppler_optical, doppler_relat
from Data_Reduction.FITS.DSNFITS import get_indices
from Data_Reduction.tipping import airmass, fit_tipcurve_data
from DatesTimes import ISOtime2datetime, UnixTime_to_datetime
from Data_Reduction import get_freq_array
from DatesTimes import MJD
# I don't think these classes are needed 2020/06/06
#try:
#  from MonitorControl.Configurations.CDSCC.FO_patching import DistributionAssembly
#  from MonitorControl.FrontEnds.K_band import K_4ch
#  K_4ch = None
#except ImportError:
#  logger.warning("SDFITSexaminer cannot handle DSS-43 K-band: no configuration")
#from MonitorControl.FrontEnds.DSN import DSN_fe
from Radio_Astronomy import rms_noise
from support import clobber, mkdir_if_needed, nearest_index
from support.dicts import make_key_if_needed


jd_start, jd_end, number = eph_manager.ephem_open()

sw_state = {True: "sig", False: "ref"}
SBcode = {1: 'U', 0:None, -1: 'L'}

class DSNFITSexaminer(object):
  """
  Class for examining SAOdataset objects
  """
  def __init__(self, parent=None, FITSfile=None, hdulist=None,
               dataname="DATA"):
    """
    Create an SAOexaminer object
    
    Attributes::
      file     - FITS file path and name
      header   - file primary HDU header
      hdulist  - list of HDUs in file
      logger   - logging.Logger object
      tables   - private class Table (SINGLE DISH) objects
      tctables - tipping curve tables
    """
    self.parent = parent
    self.logger = logging.getLogger(logger.name+'.DSNFITSexaminer')
    if FITSfile:
      self.logger.debug("__init__: opening %s", FITSfile)
      self.hdulist = pyfits.open(FITSfile)
      self.file = FITSfile
    elif hdulist:
      self.logger.debug("__init__: got %s", hdulist)
      self.hdulist = hdulist
    else:
      self.logger.error("__init__: requires filename or HDU list")
    self.header = self.hdulist[0].header
    self.tables = self.get_SDFITS_tables()
    self.logger.debug("__init__: tables found: %s", self.tables)
    if self.tables == {}:
      self.logger.warning("__init__: this has no SDFITS tables")
      if self.header['SIMPLE']: 
        self.logger.info("__init__: this is a simple FITS file")
        self.logger.info("__init__: data matrix shape is %s",
                         self.hdulist[0].data.shape)
      else:
        self.logger.warning("__init__: unknown FITS format")
    self.tctables = self.get_tipcurve_tables()
       
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
          tables[index].SB = SBcode[tables[index].get_first_good_value('CDELT1')]
        index += 1
    return tables
       
  def get_tipcurve_tables(self):
    """
    Finds the SINGLE DISH extensions in a FITS file
    """
    tctables = {}
    index = 0
    self.logger.debug("get_tipcurve_tables: checking %s", self.hdulist[1:])
    for extension in self.hdulist[1:]:
      if extension.header['extname'] == 'TIPPING CURVE':
        self.logger.debug("get_tipcurve_tables: found %s", extension)
        tctables[index] = TidTipAnalyzer(extension)
        self.logger.debug("get_tipcurve_tables: created table %s", tctables[index])
        index += 1
    return tctables
  
  def get_sources(self):
    """
    """
    sources = []
    for key in list(self.tables.keys()):
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
    # where is the save command?
     
  class Table(pyfits.BinTableHDU):
    """
    Public Attributes::
      BE          - Backend subclass associated with this table
      dss         - DSN antenna
      cycle_keys  - list of unique CYCLE values
      data        - contents of FITS binary table
      dataname    - whether it is DATA or SPECTRUM
      datestr     - like 2017/020
      DOY         - day of year
      dss         - DSN station number
      FE          - FrontEnd instance
      frame       - reference frame for X-axis
      header      - FITS table header
      logger      - logging.Logger object
      num_indices - number of indices in the DATA cube
      obsmodes    - observing modes used in this session
      obs_freqs   - non-redundant list of observing frequencies in table
      parent      - the DSNFITSexaminer object to which this belongs
      props       - table properties
      rows        - ordered list of rows with valid data
      acs_rows    - ordered list of rows for scans which have all their cycles
      scan_keys   - ordered list of scan numbers
      sources     - non-redundant list of sources in table
      timestr     - HH:MM
      year        - 4-digits
    
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
      self.parent = parent
      # add some attributes for convenience
      try:
        d = UnixTime_to_datetime(self.get_first_good_value('UNIXtime')).timetuple()
      except KeyError:
        # old deprecated format
        date_obs = self.parent.hdulist[1].get_first_good_value('DATE-OBS')
        time_obs = self.parent.hdulist[1].get_first_good_value('TIME')
        try:
          time_at_midnight = \
          time.strptime(date_obs, "%Y/%m/%d")
          d = unixtime_at_midnight + time_obs
        except ValueError:
          # acceptable format
          try:
            time_at_midnight = \
            time.strptime(date_obs, "%Y-%m-%d")
            d = unixtime_at_midnight + time_obs
          except ValueError:
            # ISOtime
            d = ISOtime2datetime(date_obs).timetuple()
      # what is the data column called?
      if 'SPECTRUM' in self.data.columns.names:
        self.dataname = 'SPECTRUM'
      elif 'DATA' in self.data.columns.names:
        self.dataname = 'DATA'
      else:
        self.logger.error("__init__: table has no DATA or SPECTRUM column")
        raise RuntimeError("no DATA or SPECTRUM column")
      self.year = d.tm_year
      self.DOY = d.tm_yday
      self.datestr = "%4d/%03d" % (self.year,self.DOY)
      self.timestr = "%02d:%02d" % (d.tm_hour, d.tm_min)
      try:
        self.dss = int(self.header['TELESCOP'].split('-')[1])
      except IndexError:
        # must be a space instead of a dash
        self.dss = int(self.header['TELESCOP'].split(' ')[1])
      self.logger.debug("__init__: station is DSS-%d", self.dss)
      self.logger.debug("__init__: front end is %s", self.header['FRONTEND'])
      if self.header['FRONTEND'] == "K":
        self.FE = K_4ch("K")
      else:
        self.FE = DSN_fe(self.header['FRONTEND'])
      try:
        #self.BE = self.header['BACKEND']
        #if self.BE == 'SAO spectrometer':
        # trying to see if this works without the actual BE class 2020/06/06
        if self.header['BACKEND'] == 'SAO spectrometer':
          # hack until SAO2SDFITS.py is fixed
          self.cycle_keys = [1]
        else:
          rows = list(range(len(self.data)))
          self.cycle_keys = list(numpy.unique(self.data['CYCLE'][rows]))
      except KeyError:
        # must be a column
        backends = self.data['BACKEND']
        if backends[0] == '':
          # deduce from number of channels
          numchans = self.data[self.dataname].shape[-1]
          # I don't think BE is used so trying without 2020/06/06
          #if numchans == 16384:
          #  self.BE = 'DAVOS'
          #elif numchans == 256:
          #  self.BE = 'SpectraData'
          #else:
          #  self.logger.error("__init__: no name for backend with %d channels",
          #                    numchans)
          #  raise RuntimeError("no name for backend")
          self.data['CYCLE'] = numpy.ones(len(self.data))
      # more convenience attributes
      self.get_table_stats()
      #self.rows = range(len(self.data))
      self.scans = numpy.unique(self.data['SCAN'])
      self.props, self.num_indices = self.properties()
      datashape = extension.data[self.dataname].shape
      self.sources = numpy.unique(self.data['OBJECT'][self.rows])
      self.sources.sort()

    def get_table_stats(self):
      """
      get number of scans, cycles, rows and observing mode
  
      Notes
      =====
      * each subchannel gets its own cycle
      * each on-source or off-source position gets its own scan 
      """
      # I suspect the following works because all are the same
      indices_for_nonzero_freqs = self.data['OBSFREQ'].nonzero()[0]
      indices_for_nonzero_cycles = self.data['CYCLE'].nonzero()[0]
      indices_for_OK_obsmodes = numpy.where(self.data['OBSMODE'] != "")[0]
      row_indices = self.data['SCAN'].nonzero()[0]
      row_indices = numpy.intersect1d(row_indices, indices_for_nonzero_freqs)
      row_indices = numpy.intersect1d(row_indices, indices_for_nonzero_cycles)
      row_indices = numpy.intersect1d(row_indices, indices_for_OK_obsmodes)
      self.rows = row_indices
      logger.debug("get_table_stats: rows with valid data: %s", self.rows)
  
      self.cycle_keys = list(numpy.unique(self.data['CYCLE'][self.rows]))
      try:
        self.cycle_keys.remove(0)
      except ValueError:
        # no CYCLE=0
        pass
      logger.debug("get_table_stats: nonzero cycles: %s", self.cycle_keys)
      num_cycles = len(self.cycle_keys)
      logger.debug("get_table_stats: %d cycles per scan", num_cycles)

      # frequencies with the same first 3 digits are the same to within 300km/s
      nonzero_subch = list(numpy.unique(self.data['OBSFREQ'][row_indices]))
      num_subch = len(nonzero_subch)
  
      if num_cycles != num_subch:
        self.logger.warning(
        "get_table_stats: number of cycles not equal to number of subchannels")
  
      # scans with valid data for at least one cycle
      scan_keys = numpy.unique(self.data['SCAN'][self.rows])
      num_scans = len(scan_keys)
  
      n_rows = num_scans * num_cycles
      if n_rows != len(self.data):
        # scans with two or more cycles, each with its own row
        diff = len(self.data) - n_rows
        logger.info("get_table_stats: there are %d scans without all its cycles",
                    diff)
        # find the scans that have all the subchannels
        complete_scans = []
        for scan in scan_keys:
          # select the rows with non-zero CYCLE for each of the scans
          scan_rows = numpy.where(self.data['SCAN'] == scan)[0]
          num_rows = len(self.data['CYCLE'][scan_rows].nonzero()[0])
          if num_rows == num_cycles:
            complete_scans.append(scan)
        logger.debug("get_table_stats: complete scans: %s", complete_scans)
        # find the rows that have complete scans
        complete_scan_rows = []
        for row in row_indices:
          this_scan = self.data['SCAN'][row]
          if this_scan in complete_scans:
            complete_scan_rows.append(row)
        self.scan_keys = list(complete_scans)
        self.acs_rows = complete_scan_rows
      else:
        # only one cycle per scan
        self.scan_keys = list(scan_keys)
        self.acs_rows = list(row_indices)
      self.logger.debug("get_table_stats: rows with complete scans: %s",
                        self.acs_rows)
      self.obsmodes = numpy.unique(self.data['OBSMODE'][self.acs_rows])
    
    def report_table(self):
      """
      """
      text = "%3d rows with valid data\n" % len(self.rows)
      text += "cycles:  %s\n" % self.cycle_keys
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
      self.logger.debug("get_obs_freqs: checking cycles %s", self.cycle_keys)
      self.logger.debug("get_obs_freqs: OBSFREQ shape is %s",
                        self.data['OBSFREQ'].shape)
      # get a frequency from each CYCLE
      self.logger.debug("get_obs_freqs: rows: %s", self.rows)
      for cycle in self.cycle_keys:
        self.logger.debug("get_obs_freqs: getting freqs for cycle %d", cycle)
        self.obs_freqs[cycle] = self.data['OBSFREQ'][self.rows]\
                         [numpy.equal(self.data['CYCLE'][self.rows], cycle)][0]
        #self.obs_freqs[cycle] = self.data['OBSFREQ'][self.rows]\
        #                                     [numpy.equal(self.rows, cycle)][0]
      return self.obs_freqs

    def get_index_keys(self):
      """
      get the row and construct a tuple for a data array index
  
      The maximum is (beam,record,IF,0,0) which will return a spectrum or whatever
      is in the first column of the multi-dimensional array.
      """
      print("These are the indices to provide to retrieve a spectrum:", end=' ')
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
  
      @param pol : IF number or STOKES axis value
      @type  pol : int
  
      @param beam : BEAM axis value
      @type  beam : int
  
      @param record : 1-based TIME axis index (FITS/FORTRAN convention) 
      @type  record : int
  
      @param trimmed : return tuple with 'RA' and 'dec' indices removed (always 0)
      @type  trimmed : bool
      """
      scan_idx = list(self.scans).index(scan)
      cycle_idx = self.cycle_keys.index(cycle)
      beam_idx = beam-1
      IF_idx = pol-1
      record -= 1
      return get_indices(self.num_indices, self.props, scan_idx=scan_idx, 
                         cycle_idx=cycle_idx, beam_idx=beam_idx, IF_idx=IF_idx,
                         record=record, trimmed=trimmed)
  
    def freqs(self, row=0, num_chans=None):
      """
      Computes frequencies for spectra in MHz
    
      Assumes scan 1 is typical of all scans in the dataset
      """
      if row == 0:
        bandwidth = self.get_first_good_value('BANDWIDT')
        ref_pix   = self.get_first_good_value('CRPIX1')
        freq      = self.get_first_good_value('CRVAL1')
        freq_step = self.get_first_good_value('CDELT1')
      else:
        bandwidth = self.data['BANDWIDT'][row]
        ref_pix   = self.data['CRVAL1'][row]
        freq      = self.data['CRVAL1'][row]
        freq_step = self.data['CDELT1'][row]
      if num_chans:
        num_chans = num_chans
        delf = bandwidth/num_chans
        refpix = ref_pix*num_chans/self.props['num chans']
      else:
        num_chans = self.props['num chans']
        delf = freq_step
        refpix = ref_pix
      freqs = freq + delf * \
               (numpy.arange(num_chans)-refpix)
      return freqs/1e6
                         
    def rel_freq_units(self,frame="FREQ-OBS", ref_freq=None, v_frame=0, row=0,
                       num_chans=None):
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
        if row == 0:
          f_ref = self.get_first_good_value('RESTFREQ')/1e6
        else:
          f_ref = self.data['RESTFREQ'][row]/1e6 # ref freq in MHz
      rel_freqs = self.freqs(row, num_chans=num_chans) - f_ref # channel
                                                 # frequencies relative
                                                 # to the reference frequency
      v_ref = self.data['VELOCITY'][row]         # vel. of object in LSR
      #self.logger.debug("rel_freq_units: reference frequency is %.3f", f_ref)
      #self.logger.debug("rel_freq_units: reference velocity is %.3f", v_ref)
      #self.logger.debug("rel_freq_units: frame velocity is %.3f", v_frame)
      self.frame = frame
      if frame == "CHAN-OBS": # channel number
        return list(range(len(rel_freqs)))
      elif frame == "FREQ-OBS": # receiver frequency in the observer's frame
        return self.freqs(row=row, num_chans=num_chans)
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
                       vspline=None, obstime=None, num_chans=None):
      """
      Computes the appropriate X-axis for a spectrum.
     
      The header parameter 'velocity' is set to the value of 'vobj' which is
      velocity of the object in the rest velocity of frame specified.
    
      Acceptable frames are defined in the self.rel_freq_units() docstring.
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
        obstime = self.get_first_good_value('UNIXtime')
      if ref_freq:
        f_ref = ref_freq
      else:
        f_ref = self.get_first_good_value('RESTFREQ')/1e6 # data['RESTFREQ'][0]/1e6 # MHz
      v_ref = self.get_first_good_value('VELOCITY') # data['VELOCITY'][0]
      if frame:
        self.frame = frame
      #self.logger.debug("compute_X_axis: requested frame is %s", frame)
      #self.logger.debug("compute_X_axis: reference frequency is %10.3f", f_ref)
      if frame == "CHAN-OBS" or frame == "FREQ-OBS" or frame == "RELA-OBS" or \
         frame == "DELF-OBS":
        # all these are in the observer's rest frame
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref,
                                num_chans=num_chans)
        vobj = None
      elif frame == "FREQ-LSR":
        vobj = self.V_LSR(row) # km/s
        C = c/1000 # km/s
        delf = -self.V_LSR(row)*f_ref/C
        #self.logger.debug(" compute_X_axis: freq offset = %f", delf)
        x = self.rel_freq_units(frame="FREQ-OBS", ref_freq=f_ref,
                                num_chans=num_chans)-delf
      elif frame == "RADI-OBS":
        # radial velocity referred to the observer
        vobj = self.V_LSR(row)
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref,
                                num_chans=num_chans)
      elif frame == "RADI-LSR":
        x = self.rel_freq_units(frame=frame, ref_freq=f_ref, 
                                v_frame=self.V_LSR(row),
                                num_chans=num_chans)
        vobj = v_ref
        #self.logger.debug("compute_X_axis: vobj = %.2f", vobj)
      elif frame == "RADI-OBJ":
        # This is the object's rest frame
        self.logger.debug("compute_X_axis: vspline = %s", vspline)
        self.logger.debug("compute_X_axis: time = %s", obstime)
        if vspline and obstime:
          # convert obstime to UNIX time
          vobj = vspline(obstime)
          x = -(c/1000)*self.rel_freq_units(frame="DELF-OBS",
                                            num_chans=num_chans)/f_ref - vobj
        else:
          vobj = self.data['VELOCITY'][row]
          x = self.rel_freq_units(frame=frame, ref_freq=f_ref,
                                  v_frame=self.V_LSR(row) + vobj,
                                  num_chans=num_chans)
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
      if row == 0:
        equinox =  self.get_first_good_value('EQUINOX')
        ra =       self.get_first_good_value('CRVAL2')
        dec =      self.get_first_good_value('CRVAL3')
        source =   self.get_first_good_value('OBJECT')
        unixtime = self.get_first_good_value('UNIXtime')
      else:
        equinox  = self.data['EQUINOX'][row]
        ra       = self.data['CRVAL2'][row]
        dec      = self.data['CRVAL3'][row]
        source   = self.data["OBJECT"][row]
        unixtime = self.data["UNIXtime"][row]
      try:
        if equinox == 1950:
          position = SkyCoord(ra, dec, frame=FK4, unit=(u.hourangle, u.deg))
      except KeyError:
        # assume J2000
        pass
      else:
        position = SkyCoord(ra, dec, frame=FK5, unit=(u.hourangle, u.deg))
      ra = position.ra.hour
      dec = position.dec.deg
      #self.logger.debug("V_LSR: ra = %f, dec = %f", ra, dec)
      cat_entry = novas.make_cat_entry(source, "", 0, ra, dec, 0, 0, 0, 0)
      source = novas.make_object(2, 0, source, cat_entry)
      longdeg = self.header['SITELONG']
      #self.logger.debug("V_LSR: longitude in degrees W = %f", longdeg)
      if longdeg > 180:
        longdeg -= 360
      #self.logger.debug("V_LSR: longitude in degrees E = %f", longdeg)
    
      DSS43 = novas.make_observer_on_surface(self.header['SITELAT'], longdeg,
                                             self.header['SITEELEV'], 0, 0)
      dt = UnixTime_to_datetime(unixtime)
      #self.logger.debug("V_LSR: computing for %s", dt.ctime())
      jd = novas.julian_date(dt.year,dt.month,dt.day,dt.hour+dt.minute/60.)
      mjd = MJD(dt.year,dt.month,dt.day)
      earth = novas.make_object(0, 3, 'Earth', None)
      urthpos,urthvel = novas.ephemeris((jd,0), earth, origin=0)
      (obspos,obsvel) = novas.geo_posvel(jd,0,DSS43,0)
      #self.logger.debug("V_LSR: Earth velocity = %s", urthvel)
      #self.logger.debug("V_LSR: observer velocity = %s", obsvel)
      totvel = tuple(numpy.array(urthvel)+numpy.array(obsvel))
      #self.logger.debug("V_LSR: total velocity = %s", totvel)
      (srcpos,srcvel) = novas.starvectors(cat_entry)
      #self.logger.debug("V_LSR: source velocity = %s", srcvel)
      V = novas.rad_vel(source, srcpos, srcvel, totvel,0,0,0)
      #self.logger.debug("V_LSR: velocity of LSR = %.2f km/s", V)
      return V+v_sun(mjd, position.ra.hour, position.dec.deg)
        
    def make_directory(self, dest=sys.stdout):
      """
      """
      labels = "Row Scan ch      Source       Sig  Freq      intg      date/time"
      flines = "--- ---- -- ---------------- ----- --------- ---- -------------------"
      lbform = "%3d  %3d %2d %16s %5s %6.0f %s %4d %20s"
      print(labels, file=dest)
      print(flines, file=dest)
      for row in self.rows:
        print(lbform % (row, self.data['SCAN'][row],
                                      self.data['CYCLE'][row],
                                      self.data['OBJECT'][row],
                                      self.data['SIG'][row],
                                      self.data['OBSFREQ'][row]/1e6,
                                      SBcode[self.data['SIDEBAND'][row]],
                                      self.data['EXPOSURE'][row],
                                      time.ctime(self.data['UNIXtime'][row])[4:]), file=dest)        
    def validate(self, colname, allow_zero=False):
      """
      """
      data = self.data[colname]
      if allow_zero:
        mask = ~numpy.isnan(data)
      else:
        mask = ~(numpy.isnan(data) | numpy.equal(data, 0))
      return mask.any(), numpy.where(mask)
        
    def validate_wx_data(self):
      """
      ensure that the weather data are usable
      
      check the data to make sure at least some are good and mask for those.
      The first item of the returned tuple reports whether any data point did
      not meet the test (of being numpy.nan or zero). The second returns the
      indices of the valid data.
      """
      valid_data = {}
      indices = {}
      #    time
      valid_data['UNIXtime'], indices['UNIXtime'] = self.validate('UNIXtime',
                                                             allow_zero=True)
      
      #    elevations
      valid_data['ELEVATIO'], indices['ELEVATIO'] = self.validate('ELEVATIO')
      if not valid_data['ELEVATIO']:
        self.logger.warning("validate_wx_data: elevation data is bad")

      #    system temperature
      valid_data['TSYS'], indices['TSYS'] = self.validate('TSYS', allow_zero=True)
                       
      #    ambient temperature
      valid_data['TAMBIENT'], indices['TAMBIENT'] = self.validate('TAMBIENT',
                                                             allow_zero=True)
      if not valid_data['TAMBIENT']:
       self.logger.warning("validate_wx_data: ambient temperature data is bad")
       
      #    pressure
      valid_data['PRESSURE'], indices['PRESSURE'] = self.validate('PRESSURE')
      if not valid_data['PRESSURE']:
        self.logger.warning("validate_wx_data: pressure data is bad")
        
      #    humidity
      valid_data['HUMIDITY'], indices['HUMIDITY'] = self.validate('HUMIDITY')
      if not valid_data['HUMIDITY']:
        self.logger.warning("validate_wx_data: humidity data is bad")
        
      #    windspeed
      valid_data['WINDSPEE'], indices['WINDSPEE'] = self.validate('WINDSPEE',
                                                             allow_zero=True)
      if not valid_data['WINDSPEE']:
        self.logger.warning("validate_wx_data: windspeed data is bad")
        
      #    winddir
      valid_data['WINDDIRE'], indices['WINDDIRE'] = self.validate('WINDDIRE',
                                                             allow_zero=True)
      if not valid_data['WINDDIRE']:
            self.logger.warning("validate_wx_data: wind direction data is bad")
      return valid_data, indices

    def get_wx_datacubes(self):
      """
      Create data for analyzing environmental conditions.
      
      Returns a dict with keys 'TAMBIENT', 'WINDDIRE', 'UNIXtime', 'TSYS',
      'HUMIDITY', 'PRESSURE', 'ELEVATIO', 'WINDSPEE'.  The data asociated with
      each key is a dict with numpy array for (SIG state) True and for False.
      The 'TSYS' array has four axes representing::
        time index   - 0-based sequence in order of matplotlib datenum 
        subchannel   - CYCLE value
        beam         - 1-based number sequence
        IF           - 1-based number sequence, usually representing pol
      The other keys have only a time axis.
      """
      # the following give rows for both cycles independent of spectrum quality
      #
      #on_rows = numpy.array(list(set(self.rows).intersection(
      #                         set(numpy.where(self.data['SIG'] == True)[0]))))
      goodtime = numpy.intersect1d(numpy.unique(self.data['UNIXtime'],
                                   return_index=True)[1],
                                   numpy.nonzero(self.data['UNIXtime'])[0])
      onweather = numpy.intersect1d(numpy.unique(self.data['UNIXtime'],
                                                  return_index=True)[1],
                                    numpy.where(self.data['SIG'] == True)[0])
      onweather = numpy.intersect1d(onweather, goodtime)
      self.logger.debug("get_wx_datacubes: rows for SIG=True: %s", onweather)
      #of_rows = numpy.array(list(set(self.rows).intersection(
      #                        set(numpy.where(self.data['SIG'] == False)[0]))))
      ofweather = numpy.intersect1d(numpy.unique(self.data['UNIXtime'],
                                                  return_index=True)[1],
                                     numpy.where(self.data['SIG'] == False)[0])
      ofweather = numpy.intersect1d(ofweather, goodtime)
      self.logger.debug("get_wx_datacubes: rows for SIG=False: %s", ofweather)
      
      valid_data, indices =  self.validate_wx_data()
      param_keys = list(valid_data.keys())
      param_keys.remove('TSYS') # to be handled separately
      
      # populate the dict with the simple parameters
      # these are the same for both cycles
      wx_data = {}
      goodonweather = numpy.intersect1d(onweather, numpy.nonzero(self.data['CYCLE'])[0])
      goodofweather = numpy.intersect1d(ofweather, numpy.nonzero(self.data['CYCLE'])[0])
      for key in param_keys:
        #wx_data[key] = {True: self.data[key][on_rows[::num_cycles]].flatten()}
        wx_data[key] = {True: self.data[key][goodonweather].flatten()}
        self.logger.debug("get_wx_datacubes: got %s for SIG=True", key)
        # it is possible that there was no position switching, e.g. pulsars
        if len(goodofweather):
          #wx_data[key][False] = self.data[key][of_rows[::num_cycles]].flatten()
          wx_data[key][False] = self.data[key][goodofweather].flatten()
          self.logger.debug("get_wx_datacubes: got %s for SIG=False", key)

      # now organize extra TSYS dimensions
      #   initialize
      wx_data['TSYS'] = {}
      # this assumes that all subchannels with a TIME axis have the same number
      # of records
      num_records = self.props['num records'][1]
      tsys_shape = {}
      goodons = numpy.intersect1d(onweather, numpy.nonzero(self.data['CYCLE'])[0])
      goodofs = numpy.intersect1d(ofweather, numpy.nonzero(self.data['CYCLE'])[0])
      #tsys_shape[True] = num_records*(len(on_rows)/self.props['num cycles'],
      tsys_shape[True] = num_records*(len(goodons),
                         self.props['num cycles'],
                         self.props['num beams'],
                         self.props['num IFs'])
      #tsys_shape[False] = num_records*(len(of_rows)/self.props['num cycles'],
      tsys_shape[False] = num_records*(len(goodofs),
                          self.props['num cycles'],
                          self.props['num beams'],
                          self.props['num IFs'])
      self.logger.debug("get_wx_datacubes: weather TSYS shape is %s", tsys_shape)
      wx_data['TSYS'][True] = numpy.zeros(tsys_shape[True])
      wx_data['TSYS'][False] = numpy.zeros(tsys_shape[False])
      # TSYS is handled differently because of a possible TIME axis
      if self.header['BACKEND'] == 'SAO spectrometer':          # has time axis
        # get the number of scans, records per scan, and sig/ref states
        if self.props['num beams'] > 1:
          time_axis = 2
        else:
          time_axis = 1
        #num_scans = self.data['TSYS'][on_rows].shape[0]
        #num_recs  = self.data['TSYS'][on_rows].shape[time_axis]
        num_scans = self.data['TSYS'][goodonweather].shape[0]
        num_recs  = self.data['TSYS'][goodonweather].shape[time_axis]
        num_on = num_scans*num_recs        # there is always SIG=True data
        #if len(self.data['SIG'][of_rows]):
        if len(self.data['SIG'][goodofweather]):
          # there are SIG=False rows
          num_states = 2
          #num_of = self.data['TSYS'][of_rows].shape[0] \
          #        *self.data['TSYS'][of_rows].shape[time_axis]
          num_of = self.data['TSYS'][goodofweather].shape[0] \
                  *self.data['TSYS'][goodofweather].shape[time_axis]
        # fill in the TSYS data
        self.logger.debug("get_wx_datacubes: %s props: %s", self, self.props)
        for subch_idx in range(self.props['num cycles']):
          for beam_idx in range(self.props['num beams']):
            for IF_idx in range(self.props['num IFs']):
              if self.props['num beams'] > 1:
                # this therefore has a time axis too
                data_on = \
                   self.data['TSYS'][goodonweather,beam_idx,:,IF_idx,0,0,0].flatten()
                   #self.data['TSYS'][on_rows,beam_idx,:,IF_idx,0,0,0].flatten()
              else:
                #data_on = self.data['TSYS'][on_rows,:,IF_idx,0,0,0].flatten()
                data_on = self.data['TSYS'][goodonweather,:,IF_idx,0,0,0].flatten()
              wx_data['TSYS'][True][:,subch_idx,beam_idx,IF_idx] = data_on
              if num_states == 2:
                if self.props['num beams'] > 1:
                  data_of = \
                   self.data['TSYS'][goodofweather,beam_idx,:,IF_idx,0,0,0].flatten()
                   #self.data['TSYS'][of_rows,beam_idx,:,IF_idx,0,0,0].flatten()
                else:
                  #data_of = self.data['TSYS'][of_rows,:,IF_idx,0,0,0].flatten()
                  data_of = self.data['TSYS'][goodofweather,:,IF_idx,0,0,0].flatten()
                wx_data['TSYS'][False][:,subch_idx,beam_idx,IF_idx] = data_of            
      elif self.header['BACKEND'][:4].upper() == 'WVSR':
        # this currently has no time axis but it has cycles, and only one beam
        #
        for subch_idx in range(self.props['num cycles']):
          subch = subch_idx+1
          self.logger.debug("get_wx_datacubes: WVSR %dth subchannel",subch_idx)
          beam_idx = 0
          for IF_idx in range(self.props['num IFs']):
            self.logger.debug("get_wx_datacubes: %dth IF", IF_idx)            
            #length = len(wx_data['TSYS'][True][:,subch_idx,beam_idx,IF_idx])
            #self.logger.debug("get_wx_datacubes: axis length is %d", length)
            wx_data['TSYS'][True][:,subch_idx,beam_idx,IF_idx] = \
              self.data['TSYS'][goodonweather,IF_idx,0,0,0].flatten()
              #self.data['TSYS'][on_rows[subch_idx::2],IF_idx,0,0,0].flatten()
              #self.data['TSYS'][on_rows[subch_idx:length:2],IF_idx,0,0,0].flatten()
            if len(goodofweather):
              wx_data['TSYS'][False][:,subch_idx,beam_idx,IF_idx] = \
                self.data['TSYS'][goodofweather,IF_idx,0,0,0].flatten()
                #self.data['TSYS'][of_rows[subch_idx::2],IF_idx,0,0,0].flatten()
                #self.data['TSYS'][of_rows[subch_idx:length:2],IF_idx,0,0,0].flatten()
      return wx_data

    def get_first_value(self, column, row):
      """
      Get the first value in a vector cell.
  
      When there is a time axis, some column cells are vectors along the time
      axis of the data cube.  Often, just the first value is needed
      corresponding to the start of the scan.  This relieves the user of having
      to known how to access that datum.
      
      @param column : name of the column
      @type  column : str
      
      @param row : number of the row in the table
      @type  row : int
      """
      try:
        cellshape = self.data[column][row].shape
      except AttributeError:
        # if the table cell has no shape then it is just a single value; for
        # example, self.data['OBJECT'][0] has no shape
        value =  self.data[column][row]
      else:
        if cellshape == ():
          # if the shape is empty the cell is a numpy array of zero dimensions
          # such as self.data['UNIXtime'][0].shape
          value = self.data[column][row]
        else:
          # create an index with the right number of indices, such as
          # tb.data['TSYS'][0].shape which has a shape (2, 1, 1, 1) so that
          # (not clear why I wrote the following this way
          #idx = len(cellshape)*[[0]]
          #return self.data[column][row][idx][0]
          idx = len(cellshape)*[0]
          value = self.data[column][row][idx]
      return value
  
    def get_first_good_value(self, column):
      """
      Get the first good value in the named column.
      
      A row with UNIX time of 0.0 does not have valid data
      """
      row = 0
      while self.data['UNIXtime'][row] == 0.0:
        row += 1
      return self.get_first_value(column, row)
      
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
        for subch in self.cycle_keys:
          subch_idx = subch-1
          spectra[scan_idx][subch_idx] = {} # 2e level dict is for subchannel
          for beam_idx in range(self.props["num beams"]):
            beam = beam_idx+1
            spectra[scan_idx][subch_idx][beam_idx] = {} # 4e levl dict for beam
            for IF_idx in range(self.props["num IFs"]): 
              pol = IF_idx+1
              spectra[scan_idx][subch_idx][beam_idx][IF_idx] = \
                                                      numpy.zeros((num_chans,))
              #self.logger.debug(
              # "prepare_summary_arrays: for scan %d subch %d, beam %d, pol %d",
              # scan, subch, beam, pol)
              #self.logger.debug("prepare_summary_arrays: spectrum shape is %s",
              #            spectra[scan_idx][subch_idx][beam_idx][IF_idx].shape)
      return spectra
      
    def prepare_summary_images(self, num_chans):
      """
      Initiate dict of image arrays for SPECTRUM data
      
      The dicts are indexed by sub-channel, beam and polarization. Each one has
      a single row of zeros.  Each record of each scan for each combination can
      then be appended to form a 2D spectrum-time plot (dynamic spectrum).
      """
      images = {}
      for subch in self.cycle_keys:
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
      all_samples = self.data[self.dataname][good_data]
      return all_samples.min(),  all_samples.max(), \
             all_samples.mean(), all_samples.std()

    def average_records(self, scan, subch, beam, pol):
      """
      average records in a scan
      
      To fill up the images arrays, the very first spectrum 
      
      @param scan : SCAN
      @param subch : CYCLE
      @param beam : 1-based number
      @param pol : 1-based number, i.e., 1 or 2
      """
      self.logger.debug(
                  "average_records: for scan %d subch%d beam%1d pol%d",
                  scan, subch, beam, pol)
       
      IF_idx = pol-1
      cycle_idx = subch-1
      num_records = self.props["num records"][subch]
      first_spectrum = True # first record's spectrum
      for record in range(num_records):
        indices = self.get_indices(scan=scan, cycle=subch, beam=beam, pol=pol,
                                   record=record)
        spectrum = self.data[self.dataname][indices].reshape(
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
   
    def get_spectra(self, rows, remove_warning=False):
      """
      returns each record for each IF as a spectrum
    
      Basically, it just removes the degenerate coordinate axes.  Use
      'get_index_keys' to get the order or the indices of the original array.
      
      Example::
        In [67]: tb0.data['SPECTRUM'].shape
        Out[67]: (12, 2, 22, 2, 1, 1, 32768)
        In [68]: specs = tb0.get_spectra(tb0.acs_rows)
        In [69]: specs.shape
        Out[69]: (12, 2, 22, 2, 32768)
        
      @param rows : a list of row numbers
      @type  rows : list of int
      
      @param remove_warning : no output warning if empty records are removed
      @type  remove_warning : bool
      """
      def remove_empty_records(spectra):
        """
        remove end records if any spectrum is empty
        
        old TAMS spectra may have records with zeros at the end, but any record
        may have a zero in the first 100 channels
        """
        while (spectra[:,:,:,:,100:] == 0).any() == True:
          # remove last records to get rid of empty records
          nrecs = spectra.shape[2]
          spectra = spectra[:,:,:-1,:,:]
          if remove_warning:
            self.logger.warning("nget_spectra: removed last records")
          # compensate integration time
          n_TAMS_recs = self.data['EXPOSURE']/5
          TAMS_intgr = 5*nrecs/n_TAMS_recs
          self.data['EXPOSURE'] -= TAMS_intgr
        return spectra
      
      nrows = len(rows)
      if (self.data['OBSMODE'][rows] == numpy.array(nrows*['LINEPBSW'])).all():
        # this requires a beam axis so...
        spectra = self.data[self.dataname][rows,:,:,:,0,0,:]
      elif (self.data['OBSMODE'][rows] == numpy.array(nrows*['LINEPSSW'])).all():
        if self.num_indices == 6:
          spectra = self.data[self.dataname][rows,:,:,:,0,0,:]
          return remove_empty_records(spectra)
        elif self.num_indices == 5:
          # has a time axis
          spectra = self.data[self.dataname][rows,:,:,0,0,:]
          return spectra
        elif self.num_indices == 4:
          # has pol axis
          spectra = self.data[self.dataname][rows,:,0,0,:]
          return spectra
        else:
          # only one spectrum
          spectra = self.data[self.dataname][rows,0,0,:]
          return spectra
      
    def normalized_beam_diff(self, rows):
      """
      (on-source - off-source)/off_source for each record and each pol
      
      This only works with spectrometers for two or more beams
      
      Example::
        In [70]: tb0.data['SPECTRUM'].shape
        Out[70]: (12, 2, 22, 2, 1, 1, 32768)
        In [71]: norms = tb0.normalized_beam_diff(tb0.acs_rows)
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
      spectra = self.get_spectra(rows)
      if self.num_indices > 4: # 5 or more axes
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
        In [74]: bpsw0 = tb0.BPSW_spectra(tb0.acs_rows)
        In [75]: bpsw0.shape
        Out[75]: (6, 22, 2, 32768)
      The rows were reduced by pairs of one on-source and one off-source.
      
      For the first years of operation, the SAO configuration did not have
      beam 1 system temperatures for the 22H sub-band. Two power meters were
      assigned to the two pols of beam 2.  To construct a Tsys array which
      uses the beam 2 data for both beam difference spectra, do this::
        In [76]: tb0.data['TSYS'][tb0.acs_rows].shape
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
      pm_keys = list(pms.keys())
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
        spectra = self.data[self.dataname][rows]
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
        rows = self.acs_rows
      scanave, Tsys, intgr = self.scans_average(rows)
      subset = self.extract_window(rows, data=scanave, Tsys=Tsys,intgr=intgr,
                                   xlimits=window, frame=frame, source=source)
      if subset:
        x, y, rms = subset.remove_baseline()
        return x, y, rms, Tsys, intgr
      else:
        self.logger.error("reduce_line: no data window found")
        return None
                                      
    def fit_mean_power_to_airmass(self, Tvac_func, first=0, last=None,
                                  replace=False):
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
    
      good_wx_data = self.get_wx_datacubes()
      paramaters = numpy.zeros_like(good_wx_data['TSYS'][True])
      std_devs = numpy.zeros_like(parameters)
      for sig in list(good_wx_data['TSYS'].keys()):
        self.logger.debug('fit_mean_power_to_airmass: processing SIG=%s',
                           sig)
        for beam_idx in range(self.props['num beams']):
          self.logger.debug('fit_mean_power_to_airmass: processing beam %d',
                              beam_idx+1)
          for IFidx in range(self.props['num IFs']):
            self.logger.debug('fit_mean_power_to_airmass: processing IF%d',
                              IFidx+1)
            pol = ['R','L'][IFidx]
            Tvac = Tvac_func(beam=beam_idx, pol=pol)
            msg = "estimated %sCP zenith Tsys in vacuum= %6.2f" % (pol, Tvac)
            self.header.add_history(msg)
            # average records here if needed.
            subch_IDs = list(range(self.props['num cycles']))
            self.logger.debug("fit_mean_power_to_airmass: subchannels: %s",
                              subch_IDs)
            for subch in subch_IDs:
              subchannel = subch+1
              # Get the data for this subchannel.
              mean_power = good_wx_data['TSYS'][sig][:,subch,beam_idx,IFidx]
              elevation  = good_wx_data['ELEVATIO'][sig]
              if last:
                E = elevation[first:last+1]
                P = mean_power[first:last+1]
              else:
                E = elevation[first:]
                P = mean_power[first:]
              # fit the data
              # estimate the slope and intercept:
              ind_max = (1/numpy.sin(numpy.deg2rad(E))).argmax()
              ind_min = (1/numpy.sin(numpy.deg2rad(E))).argmin()
              x_max = (1/numpy.sin(numpy.deg2rad(E)))[ind_max]
              x_min = (1/numpy.sin(numpy.deg2rad(E)))[ind_min]
              y_max = P[ind_max]
              y_min = P[ind_min]
              slope = (y_max-y_min)/(x_max-x_min)
              intercept = y_min - slope*x_min
              self.logger.debug(
               "fit_mean_power_to_airmass: est. intercept = %f and slope = %f",
               intercept,slope)
              #popt, pcov = curve_fit(opacity_fitting, E, P, p0=[intercept, slope])
              #intercept, slope = popt[0], popt[1]
              x = 1/numpy.sin(numpy.deg2rad(E))
              slope, intercept, r_value, p_value, std_err = linregress(x, P)
              errors = (0, std_err)
              self.logger.info("fit_mean_power_to_airmass:" +
                               " B%dP%d sch%d %s intercept, slope: %f, %f",
                               beam_idx+1, IFidx+1, subchannel, sw_state[sig],
                               intercept, slope)
              #self.logger.debug("fit_mean_power_to_airmass: covariance = %s", pcov)
              #if pcov == numpy.inf:
              #  continue
              #errors = numpy.sqrt(numpy.diag(pcov))
              msg = "IF%d, subch%d %s gain=%9.3e +/- %9.3e counts, gain_slope=%9.3e +/- %9.3e counts/am" % \
                    (IFidx+1, subchannel, sw_state[sig],
                     intercept, errors[0], slope, errors[1])
              self.header.add_history(msg)
              gain = Tvac/intercept
              gain_err = (gain/intercept)*errors[0]
              K_per_am = gain*slope
              K_per_am_err = gain*errors[1]
              self.logger.info("fit_mean_power_to_airmass:" +
                               " B%dP%d sch%d %s gain, K/am: %6.2f +/- %5.2f, %4.2f +/- %4.2f",
                               beam_idx+1, IFidx+1, subchannel, sw_state[sig],
                               gain, gain_err, K_per_am, K_per_am_err)
              if replace:
                # there are conversion constants for each switch state (SIG),
                # beam (6th data axis), IF (or pol, 4th data axis) and
                # subchannel (CYCLE)
                if len(self.data['TSYS'].shape) == 5:
                  # WVSR data
                  rowlist = list(self.select({'SIG': sig, 'CYCLE': subchannel}))
                  rowlist.sort()
                  first = rowlist[0]
                  last = rowlist[-1]+1
                  step = rowlist[1]-rowlist[0]
                  self.data['TSYS'][first:last:step, IFidx, 0, 0] *= gain
                else:
                  self.logger.warning(
                 "fit_mean_power_to_airmass: unknown shape; Tsys not computed")

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
      spectrumshape = self.data[self.dataname][0].shape # for common dims
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
          if 'TSYS' in self.data.columns.names:
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
          cycle_indices = list(range(props["num cycles"]))
          for row in cycle_indices: # check each cycle
            cycle = self.data['CYCLE'][row]
            spectrumshape = self.data[self.dataname][row].shape
            # just do the first scan
            if props["num beams"] > 1:
              props["num records"][cycle] = int(spectrumshape[1])
            else:
              props["num records"][cycle] = int(spectrumshape[0])
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
        rows = self.acs_rows
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
            # center freq is middle of channel N/2
            tone_channel = cntr_chan - chan_offset-1 + tone*tone_spacing
          else:
            tone_channel = cntr_chan - chan_offset + tone*tone_spacing
          self.logger.debug("remove_tones: tone %d channel = %d",
                            tone, tone_channel)
          if self.props['full Stokes']:
            pols = list(range(4))
          else:
            pols = self.props['num IFs']
          for pol in pols:
            self.data[self.dataname][row][pol,0,0] = clobber(
                             self.data[self.dataname][row][pol,0,0], tone_channel)

    def get_rows(self, keyword, value):
      """
      Get all the rows where the in the 'keyword' column equals 'value'
      """
      return numpy.where(self.data[keyword] == value)[0]

    def select(self, selections):
      """
      Select rows according to multiple criteria.
      
      'selections' is a dict like 
      """
      keys = list(selections.keys())
      key = keys[0]
      value = selections[key]
      selected = set(self.get_rows(key, value))
      for key in keys[1:]:
        value = selections[key]
        selected = selected.intersection(set(self.get_rows(key, value)))
      return list(selected)
      
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
 
    # these are now obsolete methods of the Table class
      
    def get_good_wx_data(self):
      """
      eliminates rows which have bad data, like 0 for CYCLE, nan, etc
      
      Note that this provides a method for flagging bad data by setting the
      CYCLE value of a row to 0.
      
      The good data flags are True if there are any good data and False if
      there are no good data for that column. They are used internally to
      construct dicts of good data
      
      Good data values are returned in the dict 'good_data'.  For most values,
      there is one key.  If the observing mode is PSSW or BPSW, then this
      key can be True (SIG) or False (REF) according to the True/False
      value in the data['SIG'] column.
      """
      if len(self.obsmodes) > 1:
        raise RuntimeErrror(
                "get_good_wx_data: multiple observing modes not yet supported")
      props, num_indices = self.properties() # table properties
      # check the data to make sure at least some are good and mask for those:
      valid_data, indices =  self.validate_wx_data()
      
      # Initialize for extracting data from simple (single value) columns
      good_data = {}
      good_data['mpltime'] = {}
      if valid_data['ELEVATIO']:
        good_data['elev'] = {}
      if valid_data['TAMBIENT']:
        good_data['Tambient'] = {}
      if valid_data['PRESSURE']:
        good_data['pressure'] = {}
      if valid_data['HUMIDITY']:
        good_data['humidity'] = {}
      if valid_data['WINDSPEE']:
        good_data['windspeed'] = {}
      if valid_data['WINDDIRE']:
        good_data['winddirec'] = {}
      if valid_data['TSYS']:
        good_data['TSYS'] = {}
      # is there any position switching?
      if self.data['SIG'].all():
        # there are no reference (off source) data
        posnsw = False
      else:
        posnsw = True
      self.logger.debug("get_good_wx_data: position switching is %s", posnsw)
      # create empty lists or dicts
      for paramkey in list(good_data.keys()):
        # there are always ONs
        if paramkey == 'TSYS':
          good_data[paramkey][True] = {}
        else:
          good_data[paramkey][True] = []
        if posnsw:
          if paramkey == 'TSYS':
            good_data[paramkey][False] = {}
          else:
            good_data[paramkey][False] = []
      self.logger.debug("get_good_wx_data: initial TSYS dict is %s", good_data['TSYS'])
      # now expand the TSYS dict
      for sig in [True, False]:
        for cycle in self.cycle_keys: # this loops over subchannels
          cycle_idx = self.cycle_keys.index(cycle)
          make_key_if_needed(good_data['TSYS'][sig], cycle_idx, value={})
          for beam_idx in range(self.props["num beams"]):
            make_key_if_needed(good_data['TSYS'][sig][cycle_idx], beam_idx, value={})
            for pol_idx in range(self.props["num IFs"]):
              make_key_if_needed(good_data['TSYS'][sig][cycle_idx][beam_idx], pol_idx, value=[])
      self.logger.debug("get_good_wx_data: expanded TSYS dict is %s",
                        good_data['TSYS'])
      # process all scans
      for row in self.acs_rows:
        # these are simple columns with multiple dimensions
        midnight_unixtime = time.mktime(time.strptime(
                                       self.data['DATE-OBS'][row], "%Y/%m/%d"))
        scan = self.data['SCAN'][row]
        cycle = self.data['CYCLE'][row]
        cycle_idx = cycle - 1
        sig = self.data['SIG'][row]
        if self.props['time axis'] == True:
          # there are time-dependent values for every record
          nrecs = self.props['num records'][cycle]
          for rec in range(nrecs):
            first_time = self.data['CRVAL5'][row]
            rectime = first_time + rec*self.data['CDELT5'][row] # numpy.array
            unixtime = midnight_unixtime + rectime
            datime = datetime.datetime.fromtimestamp(unixtime) # datetime
            # round to roughly 0.1 s
            #good_data['mpltime'][sig].append(round(date2num(datime),6))
            good_data['mpltime'][sig].append(round(datime.toordinal(),6))
            good_data['elev'][sig].append(self.data['ELEVATIO'][row,0,rec,0,0,0,0])
            for beam_idx in range(self.props["num beams"]):
              beam = beam_idx+1
              for pol_idx in range(self.props["num IFs"]):
                pol = pol_idx+1
                indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                           beam=beam, record=rec)
                good_data['TSYS'][sig][cycle_idx][beam_idx][pol_idx].append(
                                                    self.data['TSYS'][indices])
            good_data['Tambient'][sig].append(self.data['TAMBIENT'][row])
            good_data['pressure'][sig].append(self.data['PRESSURE'][row])
            good_data['humidity'][sig].append(self.data['HUMIDITY'][row])
            if valid_data['WINDSPEE']:
              good_data['windspeed'][sig].append(self.data['WINDSPEE'][row])
            if valid_data['WINDDIRE']:
              good_data['winddirec'][sig].append(self.data['WINDDIRE'][row])
        else:
          unixtime = self.data['UNIXtime'][row]
          datime = datetime.datetime.fromtimestamp(unixtime) # datetime object
          # round to roughly 0.1 s
          for beam_idx in range(self.props["num beams"]):
            beam = beam_idx+1
            for pol_idx in range(self.props["num IFs"]):
              pol = pol_idx+1
              indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                         beam=beam)
              good_data['TSYS'][sig][cycle_idx][beam_idx][pol_idx].append(
                                                    self.data['TSYS'][indices])
          if cycle == 1:
            # these are the same for all cycles
            #good_data['mpltime'][sig].append(round(date2num(datime),6))
            good_data['mpltime'][sig].append(round(datime.toordinal(),6))
            good_data['elev'][sig].append(self.data['ELEVATIO'][row])
            good_data['Tambient'][sig].append(self.data['TAMBIENT'][row])
            good_data['pressure'][sig].append(self.data['PRESSURE'][row])
            good_data['humidity'][sig].append(self.data['HUMIDITY'][row])
            if valid_data['WINDSPEE']:
              good_data['windspeed'][sig].append(self.data['WINDSPEE'][row])
            if valid_data['WINDDIRE']:
              good_data['winddirec'][sig].append(self.data['WINDDIRE'][row])
      return good_data
      
class TidTipAnalyzer():
  """
  """
  def __init__(self, extension):
    """
    """
    self.header = extension.header
    self.data = extension.data
    self.logger = logging.getLogger(logger.name+".TidTipAnalyzer")

  def fit_data(self, Tatmos=250, project=None, method="linear"):
    """
    Fits tipping curve data from a specified file
  
    @param tipfiles : dict keyed on filename with file times
    @type  tipfiles : dict
  
    @param index : file in tipfiles to be plotted
    @type  index : int
  
    @param Tatm : average temperature of the atmosphere
    @type  Tatm : float
  
    @return: float (receiver temperature), float (atmospheric opacity)
    """
    # fit the data
    Trx = {}; sigTrx = {}
    tau = {}; sigtau = {}
    Tatm = {}; sigTatm = {}
    for IF in range(4):
      PM = IF+1
      rows = numpy.where(self.data['CYCLE'] == PM)
      elev = self.data['ELEVATIO'][rows]
      tsys = self.data['TSYS'][rows]
      if method == "linear":
        Trx[IF], sigTrx[IF], Tsys_per_am, sigTsys_per_am = \
                    fit_tipcurve_data(elev, tsys, Tatm=250, method="linear")
        tau[IF] = Tsys_per_am/Tatmos
        sigtau[IF] = sigTsys_per_am/Tatmos
      elif method == "quadratic":
        Trx[IF], sigTrx[IF], tau[IF], sigtau[IF], Tatm[IF], sigTatm[IF] = \
                    fit_tipcurve_data(elev, tsys, Tatm=250, method="quadratic")
    return Trx, sigTrx, tau, sigtau, Tatm, sigTatm

