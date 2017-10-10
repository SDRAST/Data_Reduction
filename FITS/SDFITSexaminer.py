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
import astropy.units as u
import datetime
import logging
import numpy
import os
import pyfits
import sys

from astropy.coordinates import FK4, FK5, SkyCoord
from numpy import array, ndarray

# Third party packages
from novas import compat as novas
from novas.compat import eph_manager

from Astronomy import c, MJD, v_sun
from Astronomy.redshift import doppler_radio, doppler_optical, doppler_relat
from Data_Reduction.FITS.DSNFITS import get_indices, session_props
from Data_Reduction.FITS.DSNFITS import get_table_stats
from DatesTimes import UnixTime_to_datetime
from MonitorControl.BackEnds import get_freq_array
from MonitorControl.Configurations.coordinates import DSS
from support import mkdir_if_needed, nearest_index

jd_start, jd_end, number = eph_manager.ephem_open()
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
  def __init__(self, FITSfile):
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

    Extra Attributes of Tables::
      RA
      SB
      data
      dec
      dss
      header
      num_indices
      obsmodes
      props
      row_keys
      scan_keys
    These are in addition to the attributes of a pyfits.hdu.table.BinTableHDU
         
    Example::
      In [2]: run interactive.py --project=67P --DSS=43 --date=2015/204 \
                                 --console_loglevel=debug
      In [3]: E = examiner[0]
      In [4]: N = examiner[0].switched_spectra(0,1)
      In [5]: window = E.extract_window(0, spectrum=N, xlimits=(-100,100))
      In [6]: plot(window.x, window.y[0,:], label="Pol 1")
      In [7]: plot(window.x, window.y[1,:], label="Pol 2")
      In [8]: legend()
      In [9]: grid()
      In [10]: xlabel("V$_{LSR}$ (km/s)")
      In [11]: ylabel("T$_{ant}$/T$_{sys}$")
    
    @param FITSextension : SDFITS file BINTABLE extension
    @type  FITSextension : pyfits HDU
    
    @param dss : DSN station number
    @type  dss : int
    """
    self.logger = logging.getLogger(logger.name+'.DSNFITSexaminer')
    self.hdulist = pyfits.open(FITSfile)
    # hdulist[0] has the file header
    self.header = self.hdulist[0].header
    self.tables = self.get_SDFITS_tables()
    self.logger.info("__init__: %d tables found in %s",
                     len(self.tables), os.path.basename(FITSfile))
       
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
          if tables[index].dataset['CDELT1'][0] > 0:
            tables[index].SB = 'U'
          else:
            tables[index].SB = 'L'
        index += 1
    return tables
  
  class Table(pyfits.BinTableHDU):
    """
    """
    def __init__(self, parent, extension):
      """
      """
      pyfits.BinTableHDU.__init__(self, data=extension.data,
                                        header=extension.header)
      self.logger = logging.getLogger(parent.logger.name+".Table")
      self.logger.debug("__init__: created %s", self)
      # add some attributes for convenience
      self.dss = int(self.header['TELESCOP'].split('-')[1])
      self.logger.debug("__init__: station is DSS-%d", self.dss)
      self.scan_keys, self.cycle_keys, self.row_keys, self.obsmodes \
                                                        = get_table_stats(self)
      self.rows = self.data['CYCLE'].nonzero()[0] # valid data
      self.props, self.num_indices = session_props(self)
      datashape = extension.data['SPECTRUM'].shape
      if datashape == 7:
        # hack until SAO2SDFITS.py is fixed
        self.cycles = [1]
      else:
        self.cycles = list(numpy.unique(self.data['CYCLE'][self.rows]))
      self.sources = numpy.unique(self.data['OBJECT'][self.rows])
      self.report_table()

    def get_rows(self, keyword, value):
      """
      """
      return numpy.where(self.data[keyword] == value)[0]
    
    def report_table(self):
      """
      """
      self.logger.info("report_table: %d rows with valid data", len(self.rows))
      self.logger.info("report_table: cycles:  %s", self.cycles)
      self.logger.info("report_table: frequencies: %s", self.get_obs_freqs())
      self.logger.info("report_table: sources: %s", self.sources)
      self.logger.info("report_table: table properties:  %s", self.props)
  
    def get_obs_freqs(self):
      """
      """
      self.obs_freqs = {}
      for cycle in self.cycles:
        self.obs_freqs[cycle] = self.data['OBSFREQ']\
                                             [numpy.equal(self.rows, cycle)][0]
      return self.obs_freqs
    
    def get_indices(self, scan=1, cycle=1, pol=1, beam=1, IF=1, record=0):
      """
      """
      scan_idx = self.scan_keys.index(scan)
      cycle_idx = self.cycle_keys.index(cycle)
      beam_idx = beam-1
      IF_idx = IF-1
      return get_indices(self.num_indices, self.props, scan_idx=scan_idx, 
                         cycle_idx=cycle_idx, beam_idx=0, IF_idx=0, record=0)
  
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
      self.logger.info("rel_freq_units: reference frequency is %.3f", f_ref)
      rel_freqs = self.freqs(row)-f_ref          # channel frequencies relative
                                                 # to the reference frequency
      v_ref = self.data['VELOCITY'][row] # vel. of object in LSR
      self.logger.info("rel_freq_units: reference velocity is %.3f", v_ref)
      self.logger.info("rel_freq_units: frame velocity is %.3f", v_frame)
      self.frame = frame
      if frame == "CHAN-OBS":
        return range(len(rel_freqs))
      elif frame == "FREQ-OBS":
        return self.freqs(row=row)
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
  
    def compute_X_axis(self, row, frame='RADI-OBS', ref_freq=None,
                       vspline=None, obstime=None):
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
      @type  obstime : UNIX seconds)
    
      @return: numpy.ndarray
      """
      if not obstime:
        # get time from the first row
        unixtime = self.data['UNIXtime'][row]
        obstime = UnixTime_to_datetime(unixtime)
      self.logger.debug("compute_X_axis: using time %s", obstime)
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
          vobj = vspline(time.gmtime(obstime.timetuple()))
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
      self.logger.debug("V_LSR: longitude in degrees = %f", longdeg)
      if longdeg > 180:
        longdeg -= 360
      self.logger.debug("V_LSR: longitude in degrees = %f", longdeg)
    
      DSS43 = novas.make_observer_on_surface(self.header['SITELAT'], longdeg,
                                             self.header['SITEELEV'], 0, 0)
      dt = UnixTime_to_datetime(self.data['UNIXtime'][row])
      self.logger.debug("V_LSR: computing for %s", dt.ctime())
      jd = novas.julian_date(dt.year,dt.month,dt.day,dt.hour+dt.minute/60.)
      mjd = MJD(dt.year,dt.month,dt.day)
      earth = novas.make_object(0, 3, 'Earth', None)
      urthpos,urthvel = novas.ephemeris((jd,0), earth, origin=0)
      (obspos,obsvel) = novas.geo_posvel(jd,0,DSS43,0)
      self.logger.debug("V_LSR: Earth velocity = %s", urthvel)
      self.logger.debug("V_LSR: observer velocity = %s", obsvel)
      totvel = tuple(array(urthvel)+array(obsvel))
      self.logger.debug("V_LSR: total velocity = %s", totvel)
      (srcpos,srcvel) = novas.starvectors(cat_entry)
      self.logger.debug("V_LSR: source velocity = %s", srcvel)
      V = novas.rad_vel(source, srcpos, srcvel, totvel,0,0,0)
      self.logger.debug("V_LSR: velocity of LSR = %.2f km/s", V)
      return V+v_sun(mjd, position.ra.hour, position.dec.deg)
      
    def header_as_dict(self):
      """
      """
      pass
      
    def make_directory(self, dest=sys.stdout):
      """
      """
      labels = "Row Scan ch      Source       Sig Freq"
      flines = "--- ---- -- ---------------- ---- ---------"
      lbform = "%3d  %3d %2d %16s %5s %9.3f"
      print >> dest, labels
      print >> dest, flines
      for row in self.rows:
        print >> dest, lbform % (row, self.data['SCAN'][row],
                                      self.data['CYCLE'][row],
                                      self.data['OBJECT'][row],
                                      self.data['SIG'][row],
                                      self.data['OBSFREQ'][row]/1e6)        

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
      props, num_indices = session_props(self) # session properties
      # check the data:
      #    elevations
      elevation = self.data['ELEVATIO']
      mask = ~(isnan(elevation) | equal(elevation, 0))
      indices = where(mask)
      good_el_data = mask.any()
      if not good_el_data:
        self.logger.warning("get_good_rows: elevation data is bad")
      else:
        self.logger.debug("get_good_rows: elevations: %s", elevation)
      #    system temperature
      if equal(self.data['TSYS'], 0.).all():
        good_tsys_data = False
      else:
        good_tsys_data = True
        self.logger.debug("get_good_rows: system temperatures: %s",
                       self.data['TSYS'])
      #    ambient temperature
      tamb = self.data['TAMBIENT']
      mask = ~(isnan(tamb) | equal(tamb, 0))
      good_tamb_data = mask.any()
      if not good_tamb_data:
        self.logger.warning("get_good_rows: ambient temperature data is bad")
      #    pressure
      pres = self.data['PRESSURE']
      mask = ~(isnan(pres) | equal(pres, 0))
      good_pres_data = mask.any()
      if not good_pres_data:
        self.logger.warning("get_good_rows: pressure data is bad")
      #    humidity
      humi = self.data['HUMIDITY']
      mask = ~(isnan(humi) | equal(humi, 0))
      good_humi_data = mask.any()
      if not good_humi_data:
        self.logger.warning("get_good_rows: humidity data is bad")
      #    windspeed
      wspe = self.data['WINDSPEE']
      mask = ~(isnan(wspe) | equal(wspe, 0))
      good_wspe_data = mask.any()
      if not good_wspe_data:
        self.logger.warning("get_good_rows: windspeed data is bad")
      #    winddir
      wdir = self.data['WINDDIRE']
      mask = ~(isnan(wdir) | equal(wdir, 0))
      good_wdir_data = mask.any()
      if not good_wdir_data:
        self.logger.warning("get_good_rows: wind direction data is bad")
        
      # Initialize for extracting data from simple (single value) columns
      good_data = {}
      good_data['mpltime'] = []
      good_data['elev'] = []
      good_data['Tambient'] = []
      good_data['pressure'] = []
      good_data['humidity'] = []
      good_data['windspeed'] = []
      good_data['winddirec'] = []
      # process all scans
      for row in row_keys:
        # these are simple columns with multiple dimensions
        midnight_unixtime = mktime(strptime(self.data['DATE-OBS'][row],
                                          "%Y/%m/%d"))
        if self.props['time axis'] == True:
          cycle = self.data['CYCLE'][row]
          nrecs = self.props['num records'][cycle]
          for rec in range(nrecs):
            first_time = self.data['CRVAL5'][row]
            rectime = first_time + rec*self.data['CDELT5'][row] # numpy.array
            unixtime = midnight_unixtime + rectime
            datime = datetime.datetime.fromtimestamp(unixtime) # datetime.datetime
            good_data['mpltime'].append(date2num(datime))
            good_data['elev'].append(self.data['ELEVATIO'][row,0,rec,0,0,0,0])
            for beam in [0,1]:
              for pol in [0,1]:
                Tsys[beam][pol].append(self.data['Tsys'][row,beam,rec,pol,0,0,0])
            good_data['Tambient'].append(self.data['TAMBIENT'][row,0,rec,0,0,0,0])
            good_data['pressure'].append(self.data['PRESSURE'][row,0,rec,0,0,0,0])
            good_data['humidity'].append(self.data['HUMIDITY'][row,0,rec,0,0,0,0])
            good_data['windspeed'].append(self.data['WINDSPEE'][row,0,rec,0,0,0,0])
            good_data['winddirec'].append(self.data['WINDDIRE'][row,0,rec,0,0,0,0])
        else:
          good_data['unixtime'] = self.data['UNIXtime'][row]
          datime = datetime.datetime.fromtimestamp(unixtime) # datetime object
          good_data['mpltime'].append(date2num(datime))
          good_data['elev'].append(self.data['ELEVATIO'][row])
          good_data['Tambient'].append(self.data['TAMBIENT'][row])
          good_data['pressure'].append(self.data['PRESSURE'][row])
          good['humidity'].append(self.data['HUMIDITY'][row])
          good_data['windspeed'].append(self.data['WINDSPEE'][row])
          good_data['winddirec'].append(self.data['WINDDIRE'][row])
      return good_data
      
  # convert these to Table methods as needed.
                                      
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
      x_rad = numpy.deg2rad(x)
      x_sec = 1/numpy.sin(x_rad)
      return a + tau*x_sec
    
    for beam_idx in range(self.props['num beams']):
      for IFidx in range(self.props['num IFs']):
        self.logger.debug('fit_mean_power_to_airmass: processing IF%d', IFidx+1)
        Tvac = Tvac_func(beam=beam_idx ,pol=IFidx)
        pol = ['L','R'][IFidx]
        msg = "estimated %sCP zenith Tsys in vacuum= %6.2f" % (pol, Tvac)
        self.logger.debug("fit_mean_power_to_airmass: %s", msg)
        self.header.add_history(msg)
        # average records here if needed.
        subchannels = len(self.props['num cycles'])
        self.logger.debug("fit_mean_power_to_airmass: subchannels: %s",
                          subchannels)
        for subchannel in subchannels:
          subch = subchannels.index(subchannel)
          # Get the data for this subchannel.
          #   the following can be expressed as
          #     mean_power = self.data['avgpower'][subch::2,IFidx,0,0,0]
          #   or
          #     mean_power = self.data['avgpower'][:,IFidx,0,0,0][subch::2]
          #   or
          #     mean_power = self.data[subch::2]['avgpower'][:,IFidx,0,0,0]
          # assuming the first row is the same as the subchannel
          mean_power = tabhdu.data['avgpower'][subch::2,IFidx,0,0,0]
          elevation  = tabhdu.data['ELEVATIO'][subch::2]
          # get an elevation array with the 'nan' and zero values removed
          mask = ~(numpy.isnan(elevation) | numpy.equal(elevation, 0))
          # these are the indices in the data for this subchannel
          #     'where' and its friends return a tuple
          indices = numpy.where(mask)[0]
          self.logger.debug(
             "fit_mean_power_to_airmass: good data rows for subchannel %d: %s",
             subchannel, indices)
          elv = elevation[mask]
          # remove the items with the same indices from mean_power
          pwr = mean_power[mask]
          #self.logger.debug('fit_mean_power_to_airmass: elevations: %s', elv)
          self.logger.debug('fit_mean_power_to_airmass: subch %d pwr shape: %s',
                          subchannel, pwr.shape)
         # fit the data
          popt, pcov = curve_fit(opacity_fitting, elv, pwr, p0=[0, 0])
          intercept, slope = popt[0], popt[1]
          self.logger.debug(
                         "fit_mean_power_to_airmass: intercept, slope: %f, %f",
                         intercept, slope)
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
            new_indices = numpy.where(tabhdu.data['CYCLE'] == subchannel)[0]
            self.logger.debug(
                          "fit_mean_power_to_airmass: table rows for Tsys: %s",
                          new_indices)
            self.logger.debug(
                          "fit_mean_power_to_airmass: destination shape is %s",
                          tabhdu.data['TSYS'][new_indices,IFidx,0,0,0].shape)
            tabhdu.data['TSYS'][new_indices,IFidx,0,0,0] = gain * pwr
          else:
            self.logger.warning(
                        "fit_mean_power_to_airmass: failed; Tsys not computed")
   
  
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
  
  def extract_window(self, row, spectrum=None, xlimits=(0,32767), frame="RADI-LSR"):
    """
    Extracts a subset of a avg_diff spectrum
    """
    # This computes values for the entire 32768 point spectrum
    if frame:
      if frame == "RADI-LSR":
        x = self.compute_X_axis(row, "RADI-LSR")
      elif frame == "RADI-OBJ":
        if re.match('67P', self.dataset['OBJECT'][0]):
          vspline = self.load_ephemeris('67P')
          x = self.compute_X_axis(row, "RADI-OBJ", vspline=vspline)
        else:
          self.logger.error("extract_window: no ephemeris for %s", self.object)
      else:
        x = self.compute_X_axis(row, frame)
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
    if spectrum == None:
      # get the whole 6D spectrum array
      spectrum = self.dataset['SPECTRUM'][row]
    if len(spectrum.shape) == 6:
      s = spectrum[:,:,:,0,0,ch1:ch2]
    elif len(spectrum.shape) == 4:
      s = spectrum[:,:,:,ch1:ch2]
    elif len(spectrum.shape) == 3:
      s = spectrum[:,:,ch1:ch2]
    elif len(spectrum.shape) == 2:
      s = spectrum[:,ch1:ch2]
    d = Data(x[ch1:ch2], s)
    d.frame = frame
    d.channels = (ch1,ch2)
    return d
  

