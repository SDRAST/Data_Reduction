# -*- coding: utf-8 -*-
"""
FITS class - module for DSN SDFITS class

Documentation for keywords is given in::
  https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html

WARNING: DSS=43 hardcoded on line 657, 673
"""
import astropy
import astropy.units as u
import datetime
import dateutil
import logging
#import pyfits
import astropy.io.fits as pyfits
import time
import numpy

from astropy.coordinates import SkyCoord
from time import mktime, strptime

from Data_Reduction.FITS.DSNFITS import FITSfile
from MonitorControl import ObservatoryError
from support.lists import unique

logger = logging.getLogger(__name__)

obsmode = 'LINEPBSW'        # can we get this from the HDF5 file?
obsfreq = 22000000000       # Hz (wrong in HDF5 file)
pols = ['L', 'L', 'R', 'R'] # replaces P1, P2
      
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
    #self.make_prihdu()
    #self.add_site_data(self.prihdu.header, tel)
    #self.logger.debug(" Initial header: %s", self.prihdu.header)
    FITSfile.__init__(self, tel)
    self.logger = mylogger

  def tables_from_SAO(self, dss, datasets, equipment, tablesize=100,
                            project="GBRA", observer="UNKNOWN"):
    """
    Create extension for each Backend instance

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
    for BE in datasets.keys():
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
    scans = dataset.data.keys()
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
    nbeam = dims[5]
    
    self.exthead = self.make_basic_header()
    # start the extension header
    #tbhead  = pyfits.Header() ## tbhdu.header
    #tbhead['extname'] = ("SINGLE DISH", "required keyword value")
    #tbhead['nmatrix'] = (1,"one DATA column")
    #tbhead['veldef']  = ('FREQ-OBS', "raw receiver frequency")
    self.exthead['projid']  = project
    self.exthead['observer'] = observer
    self.exthead['FRONTEND'] = equipment['FrontEnd'].name
    self.exthead['RECEIVER'] = equipment['Receiver'].name
    #tbhead['TIMESYS'] = ('UTC', "DSN standard time")
    
    # add the site data
    self.add_site_data(self.exthead)
    
    # make the basic columns that are always needed
    self.make_basic_columns(numrecs=nrecs)

    # add the backend data
    self.get_hardware_metadata(BE)
    # override whatever is in the HDF5 file.  sinc(x) response
    self.exthead['FREQRES'] = \
                        1.206*dataset.header['BANDWIDT']/self.exthead['maxis1']
    
    # add multi-dimensioned metadata
    self.add_time_dependent_columns(nrecs)

    # column for TSYS, one value for each IF
    beam_dim = "(1,1,1,2,"+str(nrecs)+",2)"
    col = pyfits.Column(name='TSYS',
                         format=str(2*nrecs*2)+'D', dim=beam_dim)
    self.logger.debug("make_SAO_table: col is type %s", type(col))
    self.columns.add_col(col)
    
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
    #     Elapsed seconds since the start of UT day
    #     subtract seconds at midnight
    first_ds = dataset.header['date_obs'][0]
    start = dataset.header['start'] - time.mktime(time.strptime(
                                                   first_ds[0][0], "%Y-%m-%d"))
    #hdrindex = scans.index(scans[-1])
    #last_ds = dataset.header['date_obs'][hdrindex]
    #end   = dataset.header['end'] - time.mktime(time.strptime(
    #                                                last_ds[0][0], "%Y-%m-%d"))
    #interval = (end-start)/nrecs
    
    axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                 nrecs, 'TIME', 'E', unit='s',
                                 comment='Python time.time() value')
    #   Beam axis
    axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                                nbeam, 'BEAM',    'I',
                                                comment="beam 1 or 2")
    
    # Make the DATA column
    fmt_multiplier = self.exthead['MAXIS1']*self.exthead['MAXIS2']* \
                     self.exthead['MAXIS3']*self.exthead['MAXIS4']* \
                     self.exthead['MAXIS5']*self.exthead['MAXIS6']
    self.logger.debug("make_SAO_table: format multiplier = %d", fmt_multiplier)
    dimsval = "("+str(self.exthead['MAXIS1'])+"," \
                 +str(self.exthead['MAXIS2'])+"," \
                 +str(self.exthead['MAXIS3'])+"," \
                 +str(self.exthead['MAXIS4'])+"," \
                 +str(self.exthead['MAXIS5'])+"," \
                 +str(self.exthead['MAXIS6'])+")"
    self.logger.debug("make_SAO_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("make_SAO_table: data_format = %s", data_format)
    self.columns.add_col(pyfits.Column(name='SPECTRUM', format=data_format,
                         dim=dimsval))
    
    # create the table extension
    # tabhdu   = pyfits.new_table(cols, header=self.exthead, nrows=numscans)
    FITSrec = pyfits.FITS_rec.from_columns(self.columns, nrows=numscans)
    tabhdu = pyfits.BinTableHDU(data=FITSrec, header=self.exthead, name="SINGLE DISH")
    
    # now fill in the rows
    tabhdu = self.add_data(tabhdu, dataset, nrecs, equipment)
    return tabhdu
    
  def get_hardware_metadata(self, BE):
    """
    Initializes columns for the backends.

    Creates SDFITS columns for the backends.  Each distinct backend gets its
    own extension so BACKEND doesn't need a column.

    Notes
    =====
    Data Axes
    ---------
    Backend data are represented by a multi-dimensional cube.  For DSN FITS::
     - axis 1 = frequency
     - axis 2 = right ascension
     - axis 3 = declination
     - axis 4 = Stokes code
    Axis 2, 3 and 4 have only a single value, so that a data table is only a
    spectrum compatible with ASAP and GBTIDL.

    Beyond these there is an optional axes::
     - axis 5 = time (as for dynamic spectra)
    This makes the table incompatible with ASAP and GBTIDL but those programs
    can't handle dynamic spectra anyway.

    @param hdr : pyFITS header object

    @param cols : pyFITS columns object

    @return: tuple
      (hdr, cols) a revised header object and a revised column object
    """
    self.exthead['backend'] = BE.name
    try:
      self.exthead['firmware'] = BE.firmware
    except AttributeError:
      pass
    try:
      self.exthead['boffile']= str(BE.roach.bitstream)
    except AttributeError:
      pass
    self.exthead['maxis1'] =  (BE['num_chan'], "length of DATA axis 1")
    self.logger.debug("get_hardware_metadata: MAXIS1 = %d",
                      self.exthead['maxis1'])
    self.exthead['freqres'] =  BE['freqres']

    # GBTIDL needs a 0-based integer sequence for polarization, beam and IF
    # How many polarizations, beams and IFs?
    #pols = []
    #beams = []
    #IFs = []
    #for key in BE.inputs.keys():
    #  pols.append(BE.inputs[key].signal['pol'])
    #  beams.append(BE.inputs[key].signal['beam'])
    #  IFs.append(BE.inputs[key].name)
    #if len(pols) < 1:
    #  raise ObservatoryError(
    #      "get_hardware_metadata: ", "There are no backend outputs defined")
    #return hdr, cols
  
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
    scans = dataset.data.keys()
    numscans = len(scans)
    time_shp = (1,numrecs,1,1,1,1)
    beam_shp = (2,numrecs,2,1,1,1)
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
      tabhdu.data[index]['OBSMODE'] = obsmode
      tabhdu.data[index]['EXPOSURE'] = dataset.header['EXPOSURE'][0] # for B1P1
      tabhdu.data[index]['BANDWIDT'] = dataset.header['BANDWIDT']
      tabhdu.data[index]['TSYS'] = \
                  dataset.header['TSYS'][index][:,:numrecs,:].reshape(beam_shp)
      tabhdu.data[index]['RESTFREQ'] = dataset.header['linefreq']
      tabhdu.data[index]['VELOCITY'] = dataset.header['VELOCITY']
      # THIS IS A HACK BECAUSE dataset.header['OBSFREQ'][index] IS WRONG !!!!!!
      tabhdu.data[index]['OBSFREQ'] = obsfreq # 
      tabhdu.data[index]['VELDEF'] = dataset.header['VELDEF']
      tabhdu.data[index]['AZIMUTH'] = \
                        dataset.header['ant_az'][index][:numrecs].reshape(time_shp)
      tabhdu.data[index]['ELEVATIO'] = \
                        dataset.header['ant_el'][index][:numrecs].reshape(time_shp)
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
        logger.error("IF mode %s is not supported",
                     tabhdu.data[index]['IFmode'])
      # second and third data axes (coordinates)
      RAstr = dataset.header['R.A.'][0] # same for all records in scan
      decstr = dataset.header['declination'][0]
      coords = ' '.join((RAstr, decstr))
      c = SkyCoord(coords, unit=(u.hourangle,u.deg))
      tabhdu.data[index]['CRVAL2'] = c.ra.hourangle
      tabhdu.data[index]['CRVAL3'] = c.dec.deg
    
      # fourth data axis (polarization)
      refval, delta = polcodes(equipment['Backend'])
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
      tabhdu.data[index]['CRVAL6'] = 0
      tabhdu.data[index]['CDELT6'] = 1
      
      # the data in dataset in keyed on scan number
      tabhdu.data[index]['SPECTRUM'] = dataset.data[scan][:,:,:,:,:numrecs,:].transpose()
      
      if dataset.header['nod'][index]:
        tabhdu.data[index]['SIG'] = False
      else:
        tabhdu.data[index]['SIG'] = True
      tabhdu.data[index]['EQUINOX'] = 2000 # should come from Obs. program!
      starttime = datetime.datetime.fromtimestamp(dataset.header['time'][1][0])
      tabhdu.data[index]['VFRAME'] = dataset.V_LSR(starttime, 43) # get from obs. program!!
      tabhdu.data[index]['RVSYS'] = tabhdu.data[index]['VELOCITY']+tabhdu.data[index]['VFRAME']
      
      # unnecessary core keywords
      tabhdu.data[index]['TIME'] = tabhdu.data[index]['CRVAL5']
      
      # these still need attention
      #tabhdu.data[index]['CAL'] =
      #tabhdu.data[index]['FOFFREF1'] = 
    return tabhdu

# ------------------------------ removed from class ---------------------------
    
  def old_make_data_axis(self, header, columns,
                             daxis,
                             length,
                             dtype,
                             dformat,
                             unit=None, comment=None):
    """
    create header item and columns for a data axis

    @param tbhdu : extension (table HDU)
    @type  tbhdu : pyfits.TableHDU objetc
    
    @param daxis : axis number
    @type  daxis : int
    
    @param length : length of the axis (number of 'pixels')
    @type  length : int

    @param dtype : axis data type
    @type  dtype : str

    @param dformat : "D", "E", "I"
    @type  dformat : str
    
    @param unit : unit of measurement (defaults to SI or None)
    @type  unit : str
    """
    if comment == None:
      comment = "SPECTRUM axis "+str(daxis)
    header['ctype'+str(daxis)] = (dtype, comment)
    if unit:
      newcols = pyfits.ColDefs(
            [pyfits.Column(name='CRVAL'+str(daxis), format='1'+dformat, unit=unit),
             pyfits.Column(name='CRPIX'+str(daxis), format='1I'),
             pyfits.Column(name='CDELT'+str(daxis), format='1'+dformat, unit=unit)])
    else:
      newcols = pyfits.ColDefs(
            [pyfits.Column(name='CRVAL'+str(daxis), format='1'+dformat),
             pyfits.Column(name='CRPIX'+str(daxis), format='1I'),
             pyfits.Column(name='CDELT'+str(daxis), format='1'+dformat)])
    for col in newcols:
      columns.add_col(col)
    header['maxis'+str(daxis)] = length
    self.logger.debug("make_data_axis: MAXIS%d = %d", daxis, length)
    return header, columns

  def old_make_basic_columns(self, numrecs=1):
    """
    Make the minimum set of columns needed by SDFITS

    This make the REQUIRED columns for an SDFITS binary table::
    * SCAN     - scan number
    * CYCLE    - subscan number; increments by 1 within a row
    * DATE-OBS - ISO format date and time
    * OBJECT   - source name
    * OBSMODE  - observing mode
    * SIG      - True if on-source
    * CAL      - True if noise diode is on
    * TCAL     - effective noise diode temperature
    * EXPOSURE - integration time, sec
    * BANDWIDT - spectrometer bandwidth, Hz
    * TSYS     - system temperature
    * RESTFREQ - frequency of spectral line in the local frame
    * OBSFREQ  - center frequency of the receiver
    * VELDEF   - definition of the reference frame
    * RVSYS    - radial velocity of source w.r.t. telescope
    * VFRAME   - radial velocity of rest frame w.r.t. telescope
    * VELOCITY - radial velocity of source w.r.t. rest frame
    * EQUINOX  - epoch for source coordinates
    * FOFFREF1 - frequency offset for frequency switching
    * FRONTEND -
    * RECEIVER -
    * SAMPLER  -
    * UTIME    - a required FITS keyword

    To these are usually appended the various beam offsets available at DSN
    antennas.  See 'make_offset_columns'.

    @param numrecs : minimum of the number of records in each scan
    @type  numrecs : 1
    """
    # create empty column data arrays
    # create required columns.

    cols = pyfits.ColDefs([
      pyfits.Column(name='SCAN',     format='1I'),
      pyfits.Column(name='CYCLE',    format='1I'),
      pyfits.Column(name='DATE-OBS', format='16A'),
      pyfits.Column(name='OBJECT',   format='16A'),
      pyfits.Column(name='OBSMODE',  format='8A'),
      pyfits.Column(name='SIG',      format='1L'),
      pyfits.Column(name='CAL',      format='1L'),
      pyfits.Column(name='TCAL',     format='1E'),
      pyfits.Column(name='EXPOSURE', format='1E', unit='s'),
      pyfits.Column(name='TIME',     format='1E', unit='s'),
      pyfits.Column(name='BANDWIDT', format='1E', unit='Hz'),
      pyfits.Column(name='SIDEBAND', format='1A'),
      pyfits.Column(name='RESTFREQ', format='1D', unit='Hz'),
      pyfits.Column(name='OBSFREQ',  format='1D', unit='Hz')])
    # Velocity data
    cols += pyfits.ColDefs(
                [pyfits.Column(name='VELDEF',   format='8A'),
                 pyfits.Column(name='RVSYS',    format='1E', unit='m/s'),
                 pyfits.Column(name='VFRAME',   format='1E', unit='m/s'),
                 pyfits.Column(name='VELOCITY', format='1E', unit='m/s'),
                 pyfits.Column(name='EQUINOX',  format='1E')])

    # frequency switching offset
    cols += pyfits.ColDefs(
                     [pyfits.Column(name='FOFFREF1',  format='1E', unit='Hz')])
    
    cols += self.make_offset_columns(numrecs)

    # MORE TO ADD !!!!!!!!!!!!!!!
    #  for consideration
    # NSCANS
    return cols

  def old_add_time_dependent_columns(self, numrecs, cols):
    """
    create columns for the time-dependent metadata
    """
    time_dim = "(1,1,1,1,"+str(numrecs)+",1)" # FORTRAN order means reversed
    #cols.add_col(pyfits.Column(name='TIMESTAMP',
    #                           format=str(numrecs)+'D',
    #                           dim=time_dim))
    cols.add_col(pyfits.Column(name='LST',
                               format=str(numrecs)+'D',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='UNIXtime',
                               format=str(numrecs)+'D',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='AZIMUTH',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='ELEVATIO',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='TAMBIENT',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='PRESSURE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='HUMIDITY',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='WINDSPEE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='WINDDIRE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    return cols
    
  def old_make_offset_columns(self, numrecs, 
                          Aoff=False,Xoff=True,Eoff=True,
                          equatorial=False,
                          galactic=False,
                          refpos=False):
    """
    beam offsets such as those used to make maps

    Notes
    =====
    
    Pointing
    ~~~~~~~~
    These are not pointing offsets used by the antenna control computer to
    align the beam with the source.
    
    Columns
    ~~~~~~~
    Unused columns, that is, all zeros, should not be included so I think this
    method needs some arguments.  Do we need the REF offsets?
    These are the offsets used by DSN antennas::
      BEAMAOFF - azimuth offset
      BEAMXOFF - cross-elevation offset
      BEAMEOFF - elevation offset
      BEAMHOFF - hour angle offset
      BEAMCOFF - cross-declination offset
      BEAMDOFF - declination offset
    On 2016 Dec 22 these were added::
      BEAMLOFF - galactic longitude offset
      BEAMBOFF - galactic latitude offset
      BEAMGOFF - galactic cross-latitude offset\end{verbatim}
    """
    # always required
    if numrecs > 1:
      cols = pyfits.ColDefs([pyfits.Column(name='BEAMXOFF',
                                           format=str(numrecs)+'E',
                                          dim="(1,1,1,1,"+str(numrecs)+",1)"),
                             pyfits.Column(name='BEAMEOFF',
                                           format=str(numrecs)+'E',
                                          dim="(1,1,1,1,"+str(numrecs)+",1)")])
    else:
      cols = pyfits.ColDefs([pyfits.Column(name='BEAMXOFF',
                                           format='E', unit='deg'),
                             pyfits.Column(name='BEAMEOFF',
                                           format='E', unit='deg')])
                                           
    # the following are traditional columns.  Use above format for dimensioned
    # columns
    if Aoff: # if needed
      cols.add_col(pyfits.Column(name='BEAMAOFF', format='E', unit='deg'))
    if equatorial: # equatorial offsets
      cols.add_col(pyfits.Column(name='BEAMHOFF', format='E', unit='deg'))
      cols.add_col(pyfits.Column(name='BEAMCOFF', format='E', unit='deg'))
      cols.add_col(pyfits.Column(name='BEAMDOFF', format='E', unit='deg'))
    if galactic:# find this code
      pass 
    if refpos: # reference position for position switching
      cols.add_col(pyfits.Column(name='REF_HOFF', format='E', unit='deg'))
      cols.add_col(pyfits.Column(name='REF_DOFF', format='E', unit='deg'))
    return cols


#--------------------------------- module functions ---------------------------

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

def polcodes(BE):
  """
  Extract polarization codes from signals
  """
  pols = []
  for inpt in BE.inputs.keys():
    pols.append(BE.inputs[inpt].signal['pol'])
  pols = unique(pols)
  pols.sort()
  # convert to NRAO code
  if pols == ['E', 'H']:
    refval = -5; delta = +1
  elif pols == ['L', 'R']:
    refval = -2; delta = -1
  return refval, delta
  
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
