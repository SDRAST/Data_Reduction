# -*- coding: utf-8 -*-
"""
H5_to_FITS - module FITSfile_from_SAO and TipTableMaker classes

Documentation
=============
For FITS keywords is given in::
  https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html

Notes for further attention::
  * J2000 hard-coded in make_SAO_table() [lines 216, 219] and add_data() [450]
  * obsfreq default is 22 GHz [line 33] for K in add_data() [line 369]
  * mode default is PBSW [line 32]
  * non-Kband passband assumed to be X [line 372]
  * optional columns for CAL and FOFFREF1
"""
import astropy
import astropy.io.fits as pyfits
import astropy.units as u
import calendar
import datetime
import dateutil
import logging
import os
import time
import numpy

from astropy.coordinates import SkyCoord
from time import mktime, strptime

from Automation.NMClogs import NMClogManager,NMC_log_server
from Data_Reduction.FITS.DSNFITS import FITSfile
from DatesTimes import UnixTime_to_datetime
from MonitorControl import ObservatoryError
from MonitorControl.Configurations import station_configuration
from MonitorControl.Configurations.CDSCC.FO_patching import DistributionAssembly
from MonitorControl.Receivers.DSN import DSN_rx
from support.lists import unique

logger = logging.getLogger(__name__)

obsmode = 'LINEPBSW'        # can we get this from the HDF5 file?
obsfreq = 22000000000       # Hz (wrong in HDF5 file)
      
#############################  FITS classes  ###############################

class FITSfile_from_SAO(FITSfile):
  """
  A FITS file object having primary header and binary table extensions.
  
  The header describes the type of FITS file it is and where and when it was
  created.
  
  Each extension consists of a header followed by a column-oriented table
  (columns all have the same type of data).  The cells of the table can
  contain simple values or data arrays.  The array shape is defined in the
  header.  Columns in which all cells have the same value are considered
  'virtual' may be replaced and may be replaced by a keyword, value pair
  in the extension header.

  Public attributes::
    logger - a logging.Logger obhect
    prihdu - a pyfits.PrimaryHDU object
    tables - pyfits.BinTableHDU objects
  """
  def __init__(self, tel):
    """
    Initialize an SDFITS file with no extensions.

    Call this first to initiate an SDFITS file.  It creates the primary HDU
    suitable for a binary table extension.

    Notes
    =====

    The default length of axis 4 (STOKES) is 1 for compatibility with a
    "one spectrum per row" convention.  If this is not to be followed then
    subsequently called methods (probably Observatory.FITS_init_backends)
    must change MAXIS4.

    TDIMxx is calculated from MAXIS1 through MAXISn at the end.
    """
    mylogger = logging.getLogger(logger.name+".FITSfile")
    mylogger.debug(" creating for %s", tel.name)
    FITSfile.__init__(self, tel)
    self.logger = mylogger

  def tables_from_SAO(self, dss, datasets, equipment, tablesize=100,
                            project="GBRA", observer="UNKNOWN"):
    """
    Create extension for each Backend instance
    
    Invoked by the main program (SAO2SDFITS.py) to convert HDF5 datasets to
    SDFITS tables.  It calls 'make_SAO_table()' for each dataset.

    @param dss : DSN station
    @type  dss : Telescope instance

    @param datasets : dict of datasets for each backend
    @type  datasets : dict of SAOdataset objects

    @param tablesize : number of rows
    @type  tablesize : int

    @param project : project ID
    @type  project : str

    @param observer : team leader
    @type  observer : str
    """
    self.tables = {}
    for BE in list(datasets.keys()):
      tabhdu = self.make_SAO_table(dss, BE, datasets[BE], equipment,
                                   tablesize=tablesize,
                                   project=project, observer=observer)
      self.tables[BE] = tabhdu
  
  def make_SAO_table(self, dss, BE, dataset, equipment, tablesize=100,
                           project="GBRA", observer="UNKNOWN"):
    """
    Create extension for one Backend instance from SAOdataset

    The DATA column axes are::
    * CRVAL1 - data axis 1 reference value
    * CRPIX1 - data axis 1 reference pixel
    * CDELT1 - data axis 1 value step
    * CRVAL2 - data axis 2 reference value
    * CRPIX2 - data axis 2 reference pixel
    * CDELT2 - data axis 2 value step
    * CRVAL3 - data axis 3 reference value
    * CRPIX3 - data axis 3 reference pixel
    * CDELT3 - data axis 3 value step
    * CRVAL4 - data axis 4 reference value
    * CRPIX4 - data axis 4 reference pixel
    * CDELT4 - data axis 4 value step
    Optionally a non-degenerate TIME axis may be added but that will break
    compatibilty with ASAP and GBTIDL.  However, neither software package knows
    how to handle dynamic spectra.
    
    @param dss : DSN station
    @type  dss : Telescope instance

    @param BE : backend
    @type  BE : Backend subclass object
    
    @param tablesize : number of rows
    @type  tablesize : int

    @param project : project ID
    @type  project : str

    @param observer : team leader
    @type  observer : str

    @return: pyfits.BinTableHDU instance
    """
    # Get the dimensions of the data cube
    #  Use the first scan to get the axis length
    scans = list(dataset.data.keys())
    numscans = len(scans)
    dims = dataset.data[scans[0]].shape
    self.logger.debug("make_SAO_table: first scan shape is %s", dims)
    nchan = dims[0]
    nlong = dims[1]
    nlat  = dims[2]
    npols = dims[3]
    nrecs = 1e6 # get the short length of records in the data set
    for scan in scans:
      nrecs = min(nrecs, dataset.data[scan].shape[4])
    self.logger.debug("make_SAO_table: min. num. records: %d", nrecs)
    nbeams = dims[5]
    
    self.exthead = self.make_basic_header()
    # start the extension header
    self.exthead['projid']  = project
    self.exthead['observer'] = observer
    self.exthead['FRONTEND'] = equipment['FrontEnd'].name
    self.exthead['RECEIVER'] = equipment['Receiver'].name
    
    # adjust for X frontend and receiver
    self.logger.info("make_SAO_table: receiver is %s", equipment['Receiver'])
    if type(equipment['Receiver']) == DSN_rx:
      nbeams = 1
      self.logger.debug("make_SAO_table: DSN receivers have one beam")
    else:
      self.logger.debug("make_SAO_table: receiver has %d beams", nbeams)
    
    # add the site data
    self.add_site_data(self.exthead)
    
    # make the basic columns that are always needed
    self.make_basic_columns()

    # add the backend data
    self.get_hardware_metadata(BE)
    # override whatever is in the HDF5 file.  sinc(x) response
    self.exthead['FREQRES'] = \
                        1.206*dataset.header['BANDWIDT']/self.exthead['maxis1']
    
    # add multi-dimensioned metadata
    self.add_time_dependent_columns(nrecs)
    
    # beam offsets
    self.make_offset_columns(numrecs=nrecs)

    # column for TSYS, one value for each IF
    if nbeams == 1:
      tsys_dim = "(1,1,1,2,"+str(nrecs)+")"
      unit = "count"
    else:
      tsys_dim = "(1,1,1,2,"+str(nrecs)+",2)"
      unit = "K"
    self.logger.debug("make_SAO_table: TSYS dim is %s", tsys_dim)
    self.columns.add_col(pyfits.Column(name='TSYS', format=str(2*nrecs*2)+'D',
                                       dim=tsys_dim, unit=unit))
    
    # Add columns describing the data matrix
    #   Note that refpix defaults to 0
    axis = 1; self.make_data_axis(self.exthead, self.columns, axis,
                                  nchan, 'FREQ-OBS', 'D', unit='Hz',
                                comment="channel frequency in telescope frame")
                                
    #coordstr =  ' ' .join([dataset.header['R.A.'][0],
    #                       dataset.header['declination'][0]])
    #coord = SkyCoord(coordstr, unit=(u.hourangle, u.deg))
    
    axis +=1; self.make_data_axis(self.exthead, self.columns, axis,
                                  nlong, 'RA---GLS', 'D', unit='deg',
                                  comment="RA J2000.0")
    axis +=1; self.make_data_axis(self.exthead, self.columns, axis,
                                  nlat, 'DEC--GLS','D', unit='deg',
                                  comment="decl. J2000") 
    #   Stokes axis
    #     get the polarizations from the spectrometer input signals
    axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                 npols, 'STOKES',  'I',
                                 comment="polarization code: -8 ... 4")
    #   time axis
    axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                 nrecs, 'TIME', 'E', unit='s',
                                 comment='Python time.time() value')
    #   Beam axis
    if nbeams > 1:
      axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                                nbeams, 'BEAM',    'I',
                                                comment="beam 1 or 2")
    
    # Make the DATA column
    fmt_multiplier = self.exthead['MAXIS1']*self.exthead['MAXIS2']* \
                     self.exthead['MAXIS3']*self.exthead['MAXIS4']* \
                     self.exthead['MAXIS5']
    if nbeams > 1:
      fmt_multiplier *= self.exthead['MAXIS6']
    self.logger.debug("make_SAO_table: format multiplier = %d", fmt_multiplier)
    dimsval = "("+str(self.exthead['MAXIS1'])+"," \
                 +str(self.exthead['MAXIS2'])+"," \
                 +str(self.exthead['MAXIS3'])+"," \
                 +str(self.exthead['MAXIS4'])+"," \
                 +str(self.exthead['MAXIS5'])+")"
    if nbeams > 1:
      dimsval = dimsval[:-1] + ","+str(self.exthead['MAXIS6'])+")"
    self.logger.debug("make_SAO_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("make_SAO_table: data_format = %s", data_format)
    self.columns.add_col(pyfits.Column(name='DATA', format=data_format,
                         dim=dimsval))
    
    # create the table extension
    FITSrec = pyfits.FITS_rec.from_columns(self.columns, nrows=numscans)
    tabhdu = pyfits.BinTableHDU(data=FITSrec, header=self.exthead,
                                name="SINGLE DISH")
    
    # now fill in the rows
    tabhdu = self.add_data(tabhdu, dataset, nrecs, equipment)
    return tabhdu
  
  def add_data(self, tabhdu, dataset, numrecs, equipment):
    """
    Takes data header from an SAOdataset object and puts them into the table

    The dataset has a header which is a dict keyed on FITS-like keywords.  Each
    item is a dict indexed for scans (0-based, not scan number). Most of these
    are 1D arrays with an entry for each record. In the case of EXPOSURE the
    array is length 4, one for each pol and beam (God knows why). Some are 3D
    arrays with  additional indices for polarization and beam.  Header
    items LST, date_obs, time_obs and time are 2-dimensional arrays with one
    dimension of length 1 which must be taken into account when unpacking.
    
    """
    self.logger.debug("add_data: with %d records/scan", numrecs)
    scans = list(dataset.data.keys())
    numscans = len(scans)
    
    # add beam dimension to data shape if needed
    #   dims in C/Python order are: beam, record, pol, dec, RA, freq
    if "MAXIS6" in list(tabhdu.header.keys()):
      time_shp = (1,numrecs,1,1,1,1) # same for both beams
      beam_shp = (2,numrecs,2,1,1,1)
    else:
      time_shp = (numrecs,1,1,1,1)
      beam_shp = (numrecs,2,1,1,1)
    self.logger.debug("add_data: new shape for beam-dependent data is %s",
                      beam_shp)
    for scan in scans:
      # dataset header is keyed on the index of the scan in the set of scans
      index = scans.index(scan)
      tabhdu.data[index]['SCAN'] = scan   # int
      tabhdu.data[index]['CYCLE'] = 1     # int (0 means bad row)
      # use date of first record; see doc string for explanation of extra index
      tabhdu.data[index]['DATE-OBS'] = \
                       dataset.header['date_obs'][index][0][0].replace('-','/')
      # for convenience in Python; not FITS standard
      tabhdu.data[index]['UNIXtime'] = \
               dataset.header['time'][index][:numrecs].reshape(time_shp)
      LSTs = []
      # loop over records because of extra index
      for rec in range(numrecs):
        # see doc string for explanation of extra index
        LSTstr = dataset.header['LST'][index][rec][0]
        LSTtuple = dateutil.parser.parse(LSTstr).timetuple()
        LSTs.append((  LSTtuple.tm_hour*60
                     + LSTtuple.tm_min)*60
                     + LSTtuple.tm_sec    )
      tabhdu.data[index]['LST'] = numpy.array(LSTs).reshape(time_shp)
      tabhdu.data[index]['OBJECT'] = dataset.header['OBJECT'][index]
      tabhdu.data[index]['EXPOSURE'] = dataset.header['EXPOSURE'][0] # for B1P1
      tabhdu.data[index]['BANDWIDT'] = dataset.header['BANDWIDT']
      
      if "MAXIS6" in list(tabhdu.header.keys()):
        tabhdu.data[index]['TSYS'] = \
                  dataset.header['TSYS'][index][:,:numrecs,:].reshape(beam_shp)
        # THIS IS A HACK BECAUSE dataset.header['OBSFREQ'][index] IS WRONG !!!!
        tabhdu.data[index]['OBSFREQ'] = obsfreq # 
        tabhdu.data[index]['OBSMODE'] = obsmode
      else:
        # not K-4ch FE so no power meters and no 'off' beam
        # inspection of spectra suggests channels for averaging
        extracted_passband = dataset.data[scan][300:19300,:,:,:,:numrecs,0]
        averaged_power = extracted_passband.mean(axis=0)
        newshape = tabhdu.data[index]['TSYS'].shape
        tabhdu.data[index]['TSYS'] = \
                                   averaged_power.transpose().reshape(newshape)
        tabhdu.data[index]['OBSFREQ'] = \
                                    equipment['FrontEnd'].data['frequency']*1e6
        tabhdu.data[index]['OBSMODE'] = 'LINEPSSW'
      tabhdu.data[index]['RESTFREQ'] = dataset.header['linefreq']
      tabhdu.data[index]['VELOCITY'] = dataset.header['VELOCITY']
      tabhdu.data[index]['VELDEF'] = dataset.header['VELDEF']
      tabhdu.data[index]['AZIMUTH'] = \
                    dataset.header['ant_az'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['ELEVATIO'] = \
                    dataset.header['ant_el'][index][:numrecs].reshape(time_shp)
      # weather
      tabhdu.data[index]['TAMBIENT'] = \
                  dataset.header['TAMBIENT'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['PRESSURE'] = \
                  dataset.header['PRESSURE'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['HUMIDITY'] = \
                  dataset.header['HUMIDITY'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['WINDSPEE'] = \
                  dataset.header['WINDSPEE'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['WINDDIRE'] = \
                  dataset.header['WINDDIRE'][index][:numrecs].reshape(time_shp)
      # pointing offsets
      tabhdu.data[index]['BEAMXOFF'] = \
                  dataset.header['BEAMXOFF'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['BEAMEOFF'] = \
                  dataset.header['BEAMEOFF'][index][:numrecs].reshape(time_shp)
      
      # for first data axis (frequency)
      # THIS IS A HACK BECAUSE (see above)                               ! !!!!
      tabhdu.data[index]['CRVAL1'] = tabhdu.data[index]['OBSFREQ'] 
      tabhdu.data[index]['CDELT1'] = \
                   tabhdu.data[index]['BANDWIDT']/equipment['Backend'].num_chan
      tabhdu.data[index]['CRPIX1'] = 0
      # GBT SDFITS wants SIDEBAND
      if dataset.header['IFmode'] == 'U':
        tabhdu.data[index]['SIDEBAND'] = +1
        tabhdu.data[index]['CDELT1'] = \
                   tabhdu.data[index]['BANDWIDT']/equipment['Backend'].num_chan
      elif dataset.header['IFmode'] == 'L':
        tabhdu.data[index]['SIDEBAND'] = -1
        tabhdu.data[index]['CDELT1'] = \
                  -tabhdu.data[index]['BANDWIDT']/equipment['Backend'].num_chan
      else:
        self.logger.error("IF mode %s is not supported",
                     tabhdu.data[index]['IFmode'])
      # second and third data axes (coordinates)
      RAstr = dataset.header['R.A.'][0] # same for all records in scan
      decstr = dataset.header['declination'][0]
      coords = ' '.join((RAstr, decstr))
      c = SkyCoord(coords, unit=(u.hourangle,u.deg))
      tabhdu.data[index]['CRVAL2'] = c.ra.hourangle
      tabhdu.data[index]['CRVAL3'] = c.dec.deg
    
      # fourth data axis (polarization)
      refval, delta = equipment['Backend'].polcodes()
      tabhdu.data[index]['CRVAL4'] = refval
      tabhdu.data[index]['CDELT4'] = delta
    
      # fifth data axis (time, from records)
      # UNIX time at midnight
      midnight = time.mktime(dateutil.parser.parse(
                                   tabhdu.data[index]['DATE-OBS']).timetuple())
      seconds_from_midnight = \
               dataset.header['time'][index][:] - midnight
      # we can expect 'midnight' to change during a scan
      tabhdu.data[index]['CRVAL5'] = seconds_from_midnight[0]
      tabhdu.data[index]['CDELT5'] = ( seconds_from_midnight[numrecs-1]
                                      -seconds_from_midnight[0] )/(numrecs-1)
    
      # sixth data axis (beam)
      if "MAXIS6" in list(tabhdu.header.keys()):
        tabhdu.data[index]['CRVAL6'] = 0
        tabhdu.data[index]['CDELT6'] = 1
        # the data in dataset in keyed on scan number
        tabhdu.data[index]['DATA'] = \
                             dataset.data[scan][:,:,:,:,:numrecs,:].transpose()
      else:
        tabhdu.data[index]['DATA'] = \
                             dataset.data[scan][:,:,:,:,:numrecs,0].transpose()
      if dataset.header['nod'][index]:
        tabhdu.data[index]['SIG'] = False
      else:
        tabhdu.data[index]['SIG'] = True
      tabhdu.data[index]['EQUINOX'] = 2000 # should come from Obs. program!
      starttime = datetime.datetime.fromtimestamp(dataset.header['time'][1][0])
      tabhdu.data[index]['VFRAME'] = dataset.V_LSR(starttime, self.tel.number)
      tabhdu.data[index]['RVSYS'] = \
                    tabhdu.data[index]['VELOCITY']+tabhdu.data[index]['VFRAME']
      
      # unnecessary core keywords
      tabhdu.data[index]['TIME'] = tabhdu.data[index]['CRVAL5']

      # these still need attention
      #tabhdu.data[index]['CAL'] =
      #tabhdu.data[index]['FOFFREF1'] = 
    return tabhdu

# -------------------- convert DSS-43 tipping files to FITS -------------------

class TipTableMaker(FITSfile):
  """
  """
  def __init__(self, dss, year, doy, datefiles):
    """
    """
    self.datefiles = datefiles
    # object to get IF routing
    self.da = DistributionAssembly()
    patchdate = self.da.get_sheet_by_date("%4d/%03d" % (year,doy))
    self.IFs = self.da.get_signals('Power Meter')
    # equipment configuration
    obs, equipment = station_configuration('WBDC2_K2')
    tel = equipment['Telescope']
    FITSfile.__init__(self, tel)
    # object for getting NMC log data
    first_file = self.datefiles[0]
    first_time = float(os.path.basename(first_file)[12:-4])
    time = UnixTime_to_datetime(first_time).strftime("%H%M")
    self.NMClogman = NMClogManager(station=dss, year=year,
                            DOY=doy, time=time, use_portal=False)
    # object for providing specific data from the NMC log extracts
    self.logserver = NMC_log_server(dss, year, doy)
    
  def init_header(self):
    """
    """
    # make the table header
    self.header  = pyfits.Header() ## tbhdu.header
    self.header['extname'] = ("SINGLE DISH",          "required keyword value")
    self.header['nmatrix'] = (1,                      "one DATA column")
    self.header['veldef']  = ('FREQ-OBS',             "raw receiver frequency")
    self.header['TIMESYS'] = ('UTC',                       "DSN standard time")
    self.header['datafile'] = ('',            "path and name of original data")
    self.header['UNIXtime'] = ('',                           "UNIX time stamp")
    self.header['DATE-OBS'] = ('',             "day of observation YYYY/MM/DD")
    self.header['TIME'] = ('',                     "seconds UT since midnight")
    self.header['TAMBIENT'] = ('',                   "ambient temperature (K)")
    self.header['PRESSURE'] = ('',                 "atmospheric pressure (mB)")
    self.header['HUMIDITY'] = ('',                              "humidity (%)")
    self.header['WINDSPEE'] = ('',                         "wind speed (km/s)")
    self.header['WINDDIRE'] = ('',              "wind compass direction (deg)")
    # make the columns
    self.columns = pyfits.ColDefs([
      pyfits.Column(name='SCAN',     format='1I'), # one multichannel meas'mnt
      pyfits.Column(name='CYCLE',    format='1I'), # a single channel meas'mnt
      pyfits.Column(name='OBSFREQ',  format='1D', unit='Hz')])
    self.columns.add_col(pyfits.Column(name='BEAM', format='1I'))
    self.columns.add_col(pyfits.Column(name='POL', format='1I'))
    self.columns.add_col(pyfits.Column(name='ELEVATIO', format='1E', unit='deg'))
    self.columns.add_col(pyfits.Column(name='TSYS', format='1E', unit='K'))
    
  def make_table(self, unixtime, el, Tsys):
    """
    """
    self.init_header()
    dt = UnixTime_to_datetime(unixtime)
    nscans, nPMs = Tsys.shape
    # add header items
    self.header['UNIXtime'] = unixtime
    self.header['DATE-OBS'] = dt.strftime("%Y/%m/%d")
    self.header['TIME'] = unixtime - calendar.timegm((dt.year,1,1, 0,0,0,0,0))
    weather = self.logserver.get_weather(unixtime)
    if weather:
      self.header['TAMBIENT'] = weather[0]
      self.header['PRESSURE'] = weather[1]
      self.header['HUMIDITY'] = weather[2]
      self.header['WINDSPEE'] = weather[3]
      self.header['WINDDIRE'] = weather[4]
    
    # create the table extension
    FITSrec = pyfits.FITS_rec.from_columns(self.columns, nrows=nscans*nPMs)
    tabhdu = pyfits.BinTableHDU(data=FITSrec, header=self.header,
                                name="TIPPING CURVE")
    # now fill in the rows
    tabhdu = self.add_data(tabhdu, nscans, nPMs, el, Tsys)
    return tabhdu

  def add_data(self, tabhdu, nscans, nPMs, el, tsys):
    """
    """
    # add additional header data: TAMBIENT, PRESSURE, HUMIDITY
    index = 0
    subch_idx = 0
    for scan_idx in range(nscans):
      for pm in range(nPMs):
        tabhdu.data[index]['SCAN'] = scan_idx+1
        tabhdu.data[index]['ELEVATIO'] = el[scan_idx]
        tabhdu.data[index]['CYCLE'] = pm+1
        tabhdu.data[index]['TSYS'] = tsys[scan_idx][pm]
        sigprops = self.IFs['PM'+str(tabhdu.data[index]['CYCLE'])]
        if sigprops['Pol'] == 'H':
          tabhdu.data[index]['POL'] = -6
        elif sigprops['Pol'] == 'E':
          tabhdu.data[index]['POL'] = -5
        elif sigprops['Pol'] == 'L':
          tabhdu.data[index]['POL'] = -2
        elif sigprops['Pol'] == 'R':
          tabhdu.data[index]['POL'] = -1
        band = sigprops['Band'] = sigprops['Band']
        if sigprops['IF'] == "L":
          tabhdu.data[index]['OBSFREQ'] = (band - 0.5)*1e9
        elif sigprops['IF'] == "H":
          tabhdu.data[index]['OBSFREQ'] = (band + 0.5)*1e9
        tabhdu.data[index]['BEAM'] = sigprops['Receiver']
        index += 1
    return tabhdu

#--------------------------------- obsolete(?) module functions ---------------------------

def make_number_lookup(codes):
  """
  Given a list of codes, return a dict of 1-based ints keyed on the codes
  """
  codes = unique(codes)
  codes.sort()
  lookup = {}
  for index in range(len(codes)):
    code = codes[index]
    lookup[code] = index+1
  return codes, lookup, len(codes)
  
def make_pol_lookup(codes):
  """
  Returns a lookup table from a list of polarization codes
  """
  codes = unique(codes)
  codes.sort()
  lookup = {}
  for code in codes:
    if code == 'X' or code == 'XX' or code == 'H':
      lookup[code] = -5
    elif code == 'Y' or code == 'YY' or code == 'V' or code == 'E':
      lookup[code] = -6
    elif code == 'XY':
      lookup[code] = -7
    elif code == 'YX':
      lookup[code] = -8
    elif code == 'R' or code == 'RR':
      lookup[code] = -1
    elif code == 'L' or code == 'LL':
      lookup[code] = -2
    elif code == 'RL':
      lookup[code] = -3
    elif code == 'LR':
      lookup[code] = -4
    elif code == 'I':
      lookup[code] = 1
    elif code == 'Q':
      lookup[code] = 2
    elif code == 'U':
      lookup[code] = 3
    elif code == 'V':
      lookup[code] = 4
    else:
      # unknown
      lookup[code] = 0
  return codes, lookup, len(codes)
