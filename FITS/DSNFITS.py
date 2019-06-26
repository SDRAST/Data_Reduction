"""
Classes and functions for working with DSN SDFITS files

Could also be useful for other SDFITS files without various structures
"""
import astropy.io.fits as pyfits
import logging
import numpy
import time

from support.lists import unique

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
###################################### classes ################################

class FITSfile(object):
  """
  A FITS file object having primary header and binary table extensions.
  
  This class is a superclass with elements common to all DSN SDFITS files.
  
  The file header describes the type of FITS file it is and where and when it
  was created.
  
  Each extension consists of a header followed by a column-oriented table
  (columns all have the same type of data).  The cells of the table can
  contain simple values or data arrays.  The array shape is defined in the
  header.  Columns in which all cells have the same value are considered
  'virtual' may be replaced and may be replaced by a keyword, value pair
  in the extension header.
  
  When multiple table are created, the completed extensions are stored in
  the attribute 'tables'

  Public attributes::
    columns - columns of the extension being created
    exthead - header of the extension being created
    logger  - a logging.Logger obhect
    prihdu  - a pyfits.PrimaryHDU object
    tel     - Telescope object
    tables  - pyfits.BinTableHDU objects
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
    
    @param tel : DSN station
    @type  tel : Telescope object
    """
    self.logger = logging.getLogger(logger.name+".FITSfile")
    self.tel = tel
    self.logger.debug(" creating for %s", self.tel.name)
    self.make_prihdu()
    self.add_site_data(self.prihdu.header)
    self.tables = {}

  def make_prihdu(self):
    """
    Creates the FITS file's primary header (and data) unit
    """
    self.prihdu = pyfits.PrimaryHDU()
    self.prihdu.header['BLOCKED'] = 'T'
    self.prihdu.header['DATE'] = time.strftime("%Y/%m/%d",time.gmtime())
    self.prihdu.header['ORIGIN'] = 'FITSfile.__init__'
    
  def add_site_data(self, hdr):
    """
    Modifies 'hdr' by adding telescope data to header
    
    We pass the attribute header explicitly since so this method can be used
    to set values in the primary header and the extension header

    @param hdr : the header to be modified
    @type  hdr : pyfits Header instance
    """
    hdr['telescop'] = self.tel.name
    hdr['sitelong'] = (self.tel['longitude'], "degrees west of Greenwich")
    hdr['sitelat']  = (self.tel['latitude'],  "degrees")
    hdr['siteelev'] = (self.tel['elevation'], "meters")
    hdr['obsgeo-x'] = (self.tel['geo-x'],     "meters")
    hdr['obsgeo-y'] = (self.tel['geo-y'],     "meters")
    hdr['obsgeo-z'] = (self.tel['geo-z'],     "meters")

  # the following methods are invoked by the subclass which inherits from this
  # superclass after the subclass has initialized a new table.
    
  def make_basic_header(self):
    """
    Starts an extension header with the required values
    
    This includes values that are applicable to DSN SDFITS as well as SDFITS
    """
    header  = pyfits.Header() ## tbhdu.header
    header['extname'] = ("SINGLE DISH",               "required keyword value")
    header['nmatrix'] = (1,                           "one DATA column")
    header['veldef']  = ('FREQ-OBS',                  "raw receiver frequency")
    header['TIMESYS'] = ('UTC', "DSN standard time")
    return header
  
  def make_basic_columns(self):
    """
    Make the minimum set of columns needed by SDFITS

    This make the REQUIRED columns for an SDFITS binary table::
    * SCAN     - scan number
    * CYCLE    - subscan number; increments by 1 for every row
    * DATE-OBS - ISO format date and time
    * OBJECT   - source name
    * OBSMODE  - observing mode
    * SIG      - True if on-source
    * CAL      - True if noise diode is on
    * TCAL     - effective noise diode temperature
    * EXPOSURE - integration time, sec
    * TIME     - a required FITS keyword
    * BANDWIDT - spectrometer bandwidth, Hz
    * SIDEBAND - lower or upper with respect to the LO
    * RESTFREQ - frequency of spectral line in the local frame
    * OBSFREQ  - center frequency of the receiver
    * VELDEF   - definition of the reference frame
    * RVSYS    - radial velocity of source w.r.t. telescope
    * VFRAME   - radial velocity of rest frame w.r.t. telescope
    * VELOCITY - radial velocity of source w.r.t. rest frame
    * EQUINOX  - epoch for source coordinates    
    * FOFFREF1 - frequency offset for frequency switching

    To these are usually appended the various beam offsets available at DSN
    antennas.  See 'make_offset_columns'.
    """
    # create empty column data arrays
    # create required columns.

    self.columns = pyfits.ColDefs([
      pyfits.Column(name='SCAN',     format='1I'),
      pyfits.Column(name='CYCLE',    format='1I'), # used for subchannel
      pyfits.Column(name='DATE-OBS', format='16A'),
      pyfits.Column(name='OBJECT',   format='16A'),
      pyfits.Column(name='OBSMODE',  format='8A'),
      pyfits.Column(name='SIG',      format='1L'),
      pyfits.Column(name='CAL',      format='1L'),
      pyfits.Column(name='TCAL',     format='1E'),
      pyfits.Column(name='EXPOSURE', format='1E', unit='s'),
      pyfits.Column(name='TIME',     format='1E', unit='s'),
      pyfits.Column(name='BANDWIDT', format='1E', unit='Hz'),
      pyfits.Column(name='SIDEBAND', format='1I'),
      pyfits.Column(name='RESTFREQ', format='1D', unit='Hz'),
      pyfits.Column(name='OBSFREQ',  format='1D', unit='Hz')])
    # Velocity data
    self.columns += pyfits.ColDefs(
                [pyfits.Column(name='VELDEF',   format='8A'),
                 pyfits.Column(name='RVSYS',    format='1E', unit='m/s'),
                 pyfits.Column(name='VFRAME',   format='1E', unit='m/s'),
                 pyfits.Column(name='VELOCITY', format='1E', unit='m/s'),
                 pyfits.Column(name='EQUINOX',  format='1E')])

    # frequency switching offset
    self.columns += pyfits.ColDefs(
                     [pyfits.Column(name='FOFFREF1',  format='1E', unit='Hz')])
    
  def make_offset_columns(self, numrecs=1, 
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
      BEAMGOFF - galactic cross-latitude offset
    
    @param numrecs : number of records in a scan
    @type  numrecs : int
    
    @param Aoff : include azimuth offsets
    @type  Aoff : bool
    """
    # always required
    if numrecs > 1:
      self.columns += pyfits.ColDefs([pyfits.Column(name='BEAMXOFF',
                                           format=str(numrecs)+'E',
                                          dim="(1,1,1,1,"+str(numrecs)+",1)"),
                             pyfits.Column(name='BEAMEOFF',
                                           format=str(numrecs)+'E',
                                          dim="(1,1,1,1,"+str(numrecs)+",1)")])
    else:
      self.columns += pyfits.ColDefs([pyfits.Column(name='BEAMXOFF',
                                           format='E', unit='deg'),
                             pyfits.Column(name='BEAMEOFF',
                                           format='E', unit='deg')])
                                           
    # the following are traditional columns.  Use above format for dimensioned
    # columns
    if Aoff: # if needed
      self.columns.add_col(pyfits.Column(name='BEAMAOFF', format='E',
                                         unit='deg'))
    if equatorial: # equatorial offsets
      self.columns.add_col(pyfits.Column(name='BEAMHOFF', format='E',
                                         unit='deg'))
      self.columns.add_col(pyfits.Column(name='BEAMCOFF', format='E',
                                         unit='deg'))
      self.columns.add_col(pyfits.Column(name='BEAMDOFF', format='E',
                                         unit='deg'))
    if galactic:# find this code
      pass 
    if refpos: # reference position for position switching
      self.columns.add_col(pyfits.Column(name='REF_HOFF', format='E',
                                         unit='deg'))
      self.columns.add_col(pyfits.Column(name='REF_DOFF', format='E',
                                         unit='deg'))

  def add_time_dependent_columns(self, numrecs):
    """
    create columns for the time-dependent metadata
    
    """
    if numrecs == 1:
      time_dim = None
    else:
      # index order is spectrum, RA, dec, pol, time, beam
      time_dim = "(1,1,1,1,"+str(numrecs)+",1)" # FORTRAN order to C reverses
      # for a single beam the last index is dropped
    self.columns.add_col(pyfits.Column(name='LST',
                               format=str(numrecs)+'D',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='UNIXtime',
                               format=str(numrecs)+'D',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='AZIMUTH',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='ELEVATIO',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='TAMBIENT',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='PRESSURE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='HUMIDITY',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='WINDSPEE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    self.columns.add_col(pyfits.Column(name='WINDDIRE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    
  def make_data_axis(self, header, columns,
                     daxis, length, dtype, dformat, unit=None, comment=None):
    """
    create header item and columns for a data axis
    
    @param header : header to be modified
    @type  header : pyfits.Header object
    
    @param columns : table columns to be added
    @type  columns : pyfits.ColDefs object
    
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
      self.exthead['boffile']= str(BE.bitstream)
    except AttributeError:
      pass
    self.exthead['maxis1'] =  (BE['num_chan'], "length of DATA axis 1")
    self.logger.debug("get_hardware_metadata: MAXIS1 = %d",
                      self.exthead['maxis1'])
    self.exthead['freqres'] =  BE['freqres']

##################################### functions ###############################

def get_SDFITS_tables(hdulist):
  """
  Finds the SINGLE DISH extensions in a FITS file
  
  This is used by make_session_index.py
  """
  tables = []
  # hdulist[0] is the file header
  for extension in hdulist[1:]:
    if extension.header['extname'] == 'SINGLE DISH':
      tables.append(extension)
  return tables

def get_row(tabletype, scans, scan=0,
                       num_cycles=1, cycle=1,
                       num_tones=0, tone=0):
  """
  get the row number in a table of the specified type
  
  The bandwidth could be so narrow that there are no tones within it.
  
  @param tabletype : required 'SINGLE DISH' or 'TONES PCG', case insensitive
  @type  tabletype : str
  
  @param scans : required ordered list of scan numbers in the table
  @type  scans : list of int
  
  @param scan : SCAN number, must be in 'scans'
  @type  scan : int
  
  @param num_cycles : number of CYCLEs in a scan; default 1
  @type  num_cycles : int
  
  @param cycle : CYCLE number within a scan, starting with 1; default 1
  @type  cycle : int
  
  @param tone : tone number starting with 0 for the lowest frequency tone
  @type  tone : int
  """
  scan_idx = scans.index(scan)
  cycle_idx = cycle - 1
  if tabletype.upper() == 'SINGLE DISH':
    return num_cycles*scan_idx + cycle_idx
  elif tabletype.upper() == 'TONES PCG':
    if num_tones:
      return (scan_idx*num_cycles + cycle_idx)*num_tones + tone
    else:
      return -1 # no tone, no row
  else:
    logger.error("get_row: table type %s is not recognized", tabletype)
  
def get_indices(num_indices, props,
                scan_idx=0, cycle_idx=0, beam_idx=0, IF_idx=0, record=0,
                trimmed=False):
  """
  get the row and construct a tuple for a FITS data array index
  
  The maximum is (beam,record,IF,0,0) which will return a spectrum or whatever
  is in the first column of the multi-dimensional array.
  
  In FITS convention, lists start with 1. For Python, subtract 1 (if the list
  is continuous)
  
  This is used by class DSNFITSexaminer.
  
  @param num_indices : number of indices in the SPECTRUM array
  @type  num_indices : int
  
  @param props : table properties
  @type  props : dict
  
  @param scan_idx : index in list of SCAN numbers
  @type  scan_idx : int
  
  @param cycle_idx : index in list of CYLE numbers
  @type  cycle_idx : int
  
  @param beam_idx : index in list of BEAM axis values
  @type  beam_idx : int
  
  @param IF_idx : index in a list of IF numbers
  @type  IF_idx : int
  
  @param record : index in list of TIME axis values
  @type  record : int
  
  @param trimmed : return tuple with 'RA' and 'dec' indices removed (always 0)
  @type  trimmed : bool
  """
  row = props["num cycles"]*scan_idx + cycle_idx
  if num_indices < 3:
    raise RuntimeError("minimum data array index is 'channel, RA, dec'")
  elif num_indices == 3:
    # returns row, dec, RA
    index_tuple =  (row,0,0)
  elif num_indices == 4:
    # returns s row, IF, dec, RA
    index_tuple =  (row,IF_idx,0,0)
  elif num_indices == 5:
    index_tuple =  (row,record,IF_idx,0,0)
  elif num_indices == 6:
    index_tuple =  (row,beam_idx,record,IF_idx,0,0)
  else:
    raise RuntimeError("cannot have more that 6 axes in SPECTRUM array")
    index_tuple = (0,0)
  if trimmed:
    return index_tuple[:-2]
  else:
    return index_tuple

