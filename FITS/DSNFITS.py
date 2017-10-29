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
    self.make_offset_columns()
    
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
      time_dim = "(1,)"
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


##################################### functions ###############################

def get_SDFITS_tables(hdulist):
  """
  Finds the SINGLE DISH extensions in a FITS file
  """
  tables = []
  # hdulist[0] is the file header
  for extension in hdulist[1:]:
    if extension.header['extname'] == 'SINGLE DISH':
      tables.append(extension)
  return tables

def get_table_stats(table):
  """
  get number of scans, cycles, rows and observing mode
  
  Notes
  =====
   * each subchannel gets its own cycle
   * each on-source or off-source position gets its own scan 
  """
  logger.debug("get_table_stats: for %s", table)
  # I suspect the following works because all are the same
  indices_for_nonzero_scans = table.data['SCAN'].nonzero()
  indices_for_nonzero_cycles = table.data['CYCLE'].nonzero()
  indices_for_nonzero_freqs = table.data['OBSFREQ'].nonzero()
  # the above should probably be made into sets taking 0th element then
  # row_indices = set(indices_for_nonzero_scans).intersection(
  #     set(indices_for_nonzero_cycles).intersection(indices_for_nonzero_freqs)
  indices_for_OK_obsmodes = table.data['OBSMODE'][table.data['OBSMODE'] != ""] # not used
  # (might need to check other columns)
  row_indices = (indices_for_nonzero_scans and indices_for_nonzero_cycles and \
                 indices_for_nonzero_freqs)[0]
  # the above returns a tuple with the array of indices as the first item
  logger.debug("get_table_stats: rows without nan or zero: %s", row_indices)
  
  cycle_keys = list(numpy.unique(table.data['CYCLE'][row_indices]))
  logger.debug("get_table_stats: nonzero cycles: %s", cycle_keys)
  num_cycles = len(cycle_keys)
  logger.debug("get_table_stats: %d cycles per scan", num_cycles)
  
  nonzero_subch = list(numpy.unique(table.data['OBSFREQ'][row_indices]))
  num_subch = len(nonzero_subch)
  
  if num_cycles != num_subch:
    raise RuntimeError(
        "get_table_stats: number of cycles not equal to number of subchannels")
  
  scan_keys = numpy.unique(table.data['SCAN'][row_indices])
  logger.debug("get_table_stats: good scans: %s", scan_keys)
  num_scans = len(scan_keys)
  
  n_rows = num_scans * num_cycles
  logger.debug("get_table_stats: %d scans of %d cycles", num_scans, num_cycles)
  if n_rows != len(table.data):
    logger.info("get_table_stats: there is a scan without all its cycles")
    complete_scans = []
    # select the rows with non-zero CYCLE for each of the scans
    for scan in scan_keys:
      num_rows = len(table.data['CYCLE'][table.data['SCAN'] == scan])
      if num_rows == num_cycles:
        complete_scans.append(scan)
    logger.debug("get_table_stats: complete scans: %s", complete_scans)
    good_rows = []
    for row in row_indices:
      logger.debug("get_table_stats: processing row %s", row)
      this_scan = table.data['SCAN'][row]
      logger.debug("get_table_stats: checking scan: %s", this_scan)
      if this_scan in complete_scans:
       good_rows.append(row)
    #logger.debug("get_table_stats: good rows: %s", good_rows)
    scan_keys = list(complete_scans)
    row_keys = good_rows
  else:
    scan_keys = list(scan_keys)
    row_keys = list(row_indices)
  
  obsmodes = unique(list(table.data['OBSMODE'][row_keys]))
  return scan_keys, cycle_keys, row_keys, obsmodes

def session_props(table):
  """
  get session properties
  
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
  spectrumshape = table.data['SPECTRUM'][0].shape # for common dims
  props["num cycles"] = \
               len(numpy.unique(table.data['CYCLE'][table.data['CYCLE'] != 0]))
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
      if 'IFSPECTR' in table.data.columns.names:
        props["num IFs"] = 2
        IFspecshape = table.data['IFSPECTR'][0].shape
        props["num IFspec chans"] = int(IFspecshape[-1])
      else:
        props["num _IFs"] = 1
        logger.warning("session_props: no IF data; will use Stokes I for monitoring")
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
        cycle = table.data['CYCLE'][cycle_idx]
        spectrumshape = table.data['SPECTRUM'][cycle_idx].shape
        # just do the first scan
        props["num records"][cycle] = int(spectrumshape[1])
    else:
      props["num records"] = {1: 1}
    logger.debug("session props:\n %s", props)
    return props, len(spectrumshape)

def get_indices(num_indices, props,
                scan_idx=0, cycle_idx=0, beam_idx=0, IF_idx=0, record=0):
  """
  get the row and construct a tuple for a data array index
  
  The maximum is (beam,record,IF,0,0) which will return a spectrum or whatever
  is in the first column of the multi-dimensional array.
  
  """
  row = props["num cycles"]*scan_idx + cycle_idx
  if num_indices < 3:
    raise RuntimeError("minimum data array index is 'channel, RA, dec'")
  elif num_indices == 3:
    return (row,0,0)
  elif num_indices == 4:
    return (row,IF_idx,0,0)
  elif num_indices == 5:
    return (row,record,IF_idx,0,0)
  elif num_indices == 6:
    return (row,beam_idx,record,IF_idx,0,0)
  else:
    raise RuntimeError("cannot have more that 6 axes in SPECTRUM array")

