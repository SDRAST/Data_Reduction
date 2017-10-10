"""
functions for working with DSN SDFITS files

Could also be useful for other SDFITS files without various structures
"""
import logging
import numpy

from support.lists import unique

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

def get_first_value(table, column, row):
  """
  Get the first value in a vector cell.
  
  When there is a time axis, some column cells are vectors along the time axis
  of the data cube.  Often, just the first value is needed corresponding to the
  start of the scan.  This relieves the user of having to known how to access
  that datum.
  """
  try:
    cellshape = table.data[column][row].shape
  except AttributeError:
    return table.data[column][row]
  else:
    if cellshape == ():
      return table.data[column][row]
    else:
      idx = len(cellshape)*[[0]]
      return table.data[column][row][idx][0]

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
        logger.warning(" no IF data; will use Stokes I for monitoring")
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
    logger.debug(" session properties:\n %s", props)
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

