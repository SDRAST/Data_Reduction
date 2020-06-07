"""
module for reducing SAO spectrometer data
"""
import dill as cPickle
import logging
import re

from numpy import array, ndarray, zeros
from os.path import basename, exists

from DatesTimes import logger as dtl
dtl.setLevel(logging.WARNING)
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
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
    if type(y) == dict:
      if list(y.keys()) == [0, 1]:
        self.y = y
      else:
        self.logger.error("__init__: pols keys must be [0,1]")
    elif type(y) == list:
      if type(y[0]) == list:
        self.y = array(y)
      else:
        self.logger.error("__init__: %s is not a dict of valid Y-values")
    self.frame = None
    self.channels = None

class SAOhdf5:
  """
  Container for data in SAOspec HDF5 data files

  Notes
  =====
  **SAO HDF5 File Structure**
  
  SAO spectra are acquired from each ROACH once every five seconds.  Each of
  these is a ``record``.  Records are grouped into ``scans``.  So a 2-min
  integration will have 24 records for each of the four ROACH boards.

  A typical dataset looks like this::
  
    In [8]: data.values()
    Out[8]:
    [<HDF5 dataset "NSCANS":          shape (1, 1), type "<f4">,
     <HDF5 dataset "SITEELEV":        shape (1, 1), type "|S13">,
     <HDF5 dataset "SITELAT":         shape (1, 1), type "|S13">,
     <HDF5 dataset "SITELONG":        shape (1, 1), type "|S13">,
     <HDF5 dataset "obs_freq":        shape (1, 1), type "<f4">,
     <HDF5 dataset "rest_freq":       shape (1, 1), type "<f4">,
     <HDF5 dataset "v_ref":           shape (1, 1), type "|S10">,
     <HDF5 dataset "vsys":            shape (1, 1), type "<f4">,
     <HDF5 dataset "integ_time":      shape (1, 4), type "<f4">,
     <HDF5 dataset "scan_number":     shape (326, 1), type "<f4">,
     <HDF5 dataset "date_obs":        shape (326, 1), type "|S10">,
     <HDF5 dataset "time_obs":        shape (326, 1), type "|S10">,
     <HDF5 dataset "timestamp":       shape (326, 1), type "<i8">,
     <HDF5 dataset "LST":             shape (326, 1), type "|S13">,
     <HDF5 dataset "scan_duration":   shape (326, 1), type "<f4">,
     <HDF5 dataset "observer":        shape (326, 1), type "|S10">,
     <HDF5 dataset "source_name":     shape (326, 1), type "|S10">,
     <HDF5 dataset "onsource":        shape (326, 1), type "|b1">,
     <HDF5 dataset "current_azel":    shape (326, 2), type "<f4">,
     <HDF5 dataset "offsets":         shape (326, 2), type "<f4">,
     <HDF5 dataset "source_azel":     shape (326, 2), type "|S12">,
     <HDF5 dataset "source_long_lat": shape (326, 2), type "|S13">,
     <HDF5 dataset "source_radec":    shape (326, 2), type "|S13">,
     <HDF5 dataset "weather":         shape (326, 6), type "<f4">]
     <HDF5 dataset "Tsys":            shape (326, 4), type "<f4">,
     <HDF5 dataset "bandwidth":       shape (326, 4), type "<f4">,
     <HDF5 dataset "mode":            shape (326, 4), type "|S10">,
     <HDF5 dataset "pol":             shape (326, 4), type "|S3">,
     <HDF5 dataset "spectraCh1":      shape (326, 32768), type "<f4">,
     <HDF5 dataset "spectraCh2":      shape (326, 32768), type "<f4">,
     <HDF5 dataset "spectraCh3":      shape (326, 32768), type "<f4">,
     <HDF5 dataset "spectraCh4":      shape (326, 32768), type "<f4">,

  In order, the first set are parameters which do not change for the entire
  dataset. ``integ_time`` is specified for each backend.  The remaining
  parameters are different for each integration. Most of the remaining
  arrays are arrays with multiple values in an obvious way such as five
  sec integrations, or parameters
  for each DSProc (ROACH).  The data associated with each nominally 5 sec are
  associated with scans like this::
  
    In [49]: hdf.data['weather'].value
    Out[49]:
    array([[ 28.05555725, 946.18383789,  21.9917202 ,  9.88373375, 187.,  0.],
           [ 28.05555725, 946.18383789,  21.24612045,  9.26599979, 157.,  0.],
           [ 28.05555725, 946.18383789,  21.24612045,  9.26599979, 157.,  0.],
           ...,
           [ 28.38888931, 945.9130249 ,  19.21655655,  5.55959988,  42.,  0.],
           [ 28.38888931, 945.9130249 ,  18.9209938 ,  5.55959988,  44.,  0.],
           [  0.        ,   0.        ,   0.        ,  0.        ,   0.,  0.]],
          dtype=float32)

  The records are associated with scans like this::
  
    In [48]: hdf.data['scan_number'].value
    Out[48]:
    array([[ -1.],[ -1.],[ -1.],[ -1.],[ -1.],[ -1.],[ -1.],[ -1.],[-1.],[-1.],
           [  1.],[  1.],[  1.],[  1.],[  1.],[  1.],[  1.],[  1.],[ 1.],[ 1.],
           [  1.],[  1.],
           [ -1.],[ -1.],
           [  2.],[  2.],[  2.],[  2.],[  2.],[  2.],[  2.],[  2.],[ 2.],[ 2.],
           [  2.],[  2.],
           [ -1.],[ -1.],[ -1.],
            ...
           [ -1.],[ -1.],
           [ 19.],[ 19.],[ 19.],[ 19.],[ 19.],[ 19.],[ 19.],[ 19.],[19.],[19.],
           [ 19.],[ 19.],[ 19.],
           [ -1.],[ -1.],
           [ 20.],[ 20.],[ 20.],[ 20.],[ 20.],[ 20.],[ 20.],[ 20.],[20.],[20.],
           [ 20.],[ 20.],
           [ -1.],
           [  0.]], dtype=float32)

  This class has an attribute 'container' similar to an SDFITS table.
  The 'spec' axes are::
  
    frequency   - 32768 equally spaced frequencies in the topocentric frame
    R.A.        -     1 J2000 coordinate
    declination -     1 J2000 coordinate
    pol         -     2 for E,H or L,R
    time        -     N seconds since midnight for each record in a scan
    beam        -     2 for 1,2

  Attributes
  ==========
    clobber_0      - set 0th value to 1st value if True
    clobber_center - set 16384th value to mean of 16383th and 16385th if True
    container      - SAOdataset object in which the reformatted data are stored
    data           - HDF5 file contents (temporary)
    filename       - name of HDF5 file
    logger         - logging.Logger instance
    meta           - metadata from HDF5 file
    num_DSP        - number of signal processors
    num_uniq       - number of unique records for each scan and beam and pol
    PROJID         - identifier string for project
    SITEELEV       - elevation of telescope (m)
    SITELAT        - latitude of telescope (deg)
    SITELONG       - east longitude of telescope (deg)
    spec           - nparray of spectra
    uniq_recs      - dict with number of unique records for each pol and beam

  Methods
  =======
    extract_valid_records - gets non-duplicate, non-zero records
    to_dataset            - returns an SAOdataset instance
  """
  def __init__(self, filename, clobber_0=True, clobber_center=True):
    """
    Create as SAOspec data container

    @param filename : SAO dataset, original hdf5 (.h5) or Python pickle (.pkl)
    @type  filename : str

    @param clobber_0 : replace ch 0 with ch 1
    @type  clobber_0 : bool

    @param clobber_center : replace ch 16384 with average of adjacents
    @type  clobber_center : bool
    """
    self.logger =  logging.getLogger(logger.name+".SAOhdf5")
    self.filename = filename
    self.logger.debug("__init__: processing file %s", self.filename)
    self.clobber_0 = clobber_0
    self.clobber_center = clobber_center
    self.container = SAOdataset()
    if self._load_() == None:
      raise RuntimeError("could not create SAOhdf5 object")
    self.uniq_recs, self.rec_time, self.num_uniq = self.extract_valid_records()
    self.container.scan_keys = list(self.uniq_recs.keys())
    self.container.scan_keys.sort()
    self.container.first_scan = self.container.scan_keys[0]

  def _load_(self):
    """
    Initialize an SAOhdf5 object from an HDF5 file

    @param reload : if False, find and load the associated pickle file
    @type  reload : bool
    """
    filetype = splitext(basename(self.filename))[1]
    self.logger.debug("_load_: opening %s file %s", filetype, self.filename)
    try:
      self.data = h5py.File(self.filename,'r') # HDF5 group
    except IOError as details:
      if re.search("unable to open file",str(details)):
        self.logger.error("_load_: missing or badly closed file")
      else:
        self.logger.error("_load_: %s", str(details))
      raise RuntimeError("_load_: failed")
    self.logger.debug("_load_: keys: %s", list(self.data.keys()))
    self.scan_number = self.data['scan_number'].value # array for each record
    self.num_recs = len(self.scan_number) # scan number of each record
    self.numDSP = self.data['integ_time'].value.shape[1] # should be 4
    self._read_metadata_file_()
    self._unpack_spectra_()
    self.logger.debug("_load_: loaded %s", self.filename)
    return True

  def _unpack_spectra_(self, beamsw=False):
    """
    Unpacks the data in an SAO HDF5 file.

    A scan is a series of integrations. Scans are numbered sequentially from 1.
    A scan number -1 indicates that the antenna is moving between feeds.
    The first feed selection is considered 'nod 0', the other 'nod 1'.
    (The feed selection can also be obtained from offsets.)

    Note that this is not for double-Dicke switching as described in the
    module docstring. Nothing happens with the antenna or receiver from one
    record (integration) to the next.  All the records in a scan for one
    ROACH may be considered parts of one long integration.  Such an integration
    is called a 'spec' here.

    This still needs code added for the case 'beamsw=True'.  This is used if
    the cross-over switch is cycled between every pair of records.

    @param beamsw - True for beam switching, False for position switching
    @type  beamsw - bool
    """
    def _apply_metadata_(self, groupIndex, records):
      """
      Add metadata to scan spectra.

      @param groupIndex - 0 for first scan, incrementing by 1
      @type  groupIndex - int

      @param records - list with indices for the records in each scan
      @type  records - list of int
      """
      group = groupIndex+1
      self.logger.debug(
      "_unpack_spectra_._apply_metadata_: Getting group (scan) %d metadata with records %s",
       group, records)
      self.meta[groupIndex]['LST'] = self.data['LST'][records[0]:records[-1]]
      self.SITEELEV = float(self.data['SITEELEV'][0,0])
      self.SITELAT = float(self.data['SITELAT'][0,0])
      self.SITELONG = float(self.data['SITELONG'][0,0])

      self.meta[groupIndex]['BANDWIDT'] = \
                                       self.data['bandwidth'][records[0],0]*1e6
      self.meta[groupIndex]['ant_az'] = \
                            self.data['current_azel'][records[0]:records[-1],0]
      self.meta[groupIndex]['ant_el'] = \
                            self.data['current_azel'][records[0]:records[-1],1]
      self.meta[groupIndex]['date_obs'] = \
                                  self.data['date_obs'][records[0]:records[-1]]
      self.meta[groupIndex]['time_obs'] = \
                                  self.data['time_obs'][records[0]:records[-1]]
      self.meta[groupIndex]['time'] = \
                                 self.data['timestamp'][records[0]:records[-1]]
      self.meta[groupIndex]['EXPOSURE'] = \
                                      self.data['integ_time'][0,:]*len(records)

      self.meta[groupIndex]['linefreq'] = self.data['rest_freq'][0,0]
      # now we make some intelligent guesses about parameters not in HDF5 file
      if 'band' in self.meta:
        self.meta[groupIndex]['band'] = self.meta['band']
        self.logger.warning(
         "_unpack_spectra_._apply_metadata_: override of %s applied to 'band'",
         self.meta['band'])
      else:
        self.meta[groupIndex]['band'] = \
                                     round(self.meta[0]['linefreq']/2e9)*2
      if self.meta[groupIndex]['linefreq'] < self.meta[groupIndex]['band']*1e9:
        self.meta[groupIndex]['IFmode'] = "L"
      else:
        self.meta[groupIndex]['IFmode'] = "U"
      if self.meta[groupIndex]['IFmode'] == "U":
        obsfreq = self.meta[groupIndex]['band']*1e9 + 500e6
      elif self.meta[groupIndex]['IFmode'] == "L":
        obsfreq = self.meta[groupIndex]['band']*1e9 - 500e6
      else:
        obsfreq = self.meta[groupIndex]['band']*1e9
      self.meta[groupIndex]['OBSFREQ'] = obsfreq
      self.PROJID = self.data['observer'][records[0],0]

      self.meta[groupIndex]['BEAM'] = []
      self.meta[groupIndex]['pol'] = []
      self.logger.debug("_unpack_spectra_._apply_metadata_: processing %s",
                                                 self.data['mode'][records[0]])
      num_IFs = len(self.data['mode'][records[0]])
      self.meta[groupIndex]['TSYS'] = empty((2,len(records),2))
      self.logger.debug("Created TSYS for %d with shape %s", groupIndex,
                                           self.meta[groupIndex]['TSYS'].shape)
      for key in range(num_IFs):
        # beam and pol information for the 4 ROACH boards are encoded as 'mode'
        # like this:  ['A1P1' 'A2P1' 'A1P2' 'A2P2']
        beam_code = self.data['mode'][records[0]][key][:2]
        pol_code  = self.data['mode'][records[0]][key][2:]
        self.meta[groupIndex]['BEAM'].append(beam_code)
        self.meta[groupIndex]['pol'].append(pol_code)
        pol_index  = int(self.meta[groupIndex]['pol'][key][1])-1
        beam_index = int(self.meta[groupIndex]['BEAM'][key][1])-1
        # format Tsys into an array indexed by pol, record, beam; same order as
        # main data array
        for record in records:
          index = records.index(record)
          self.meta[groupIndex]['TSYS'][pol_index, index, beam_index] = \
                                         self.data['Tsys'].value[record,key]
      self.logger.debug(
                  "_unpack_spectra_._apply_metadata_: %d'th scan beams are %s",
                  groupIndex, self.meta[groupIndex]['BEAM'])
      self.logger.debug(
                   "_unpack_spectra_._apply_metadata_: %d'th scan pols are %s",
                   groupIndex, self.meta[groupIndex]['pol'])
      # antenna status
      self.meta[groupIndex]['BEAMXOFF'] = \
                                 self.data['offsets'][records[0]:records[-1],0]
      self.meta[groupIndex]['BEAMEOFF'] = \
                                 self.data['offsets'][records[0]:records[-1],1]
      self.meta[groupIndex]['tracking'] = \
                                self.data['onsource'][records[0]:records[-1],0]
      # what is self.data['scan_duration']?
      self.meta[groupIndex]['SCAN'] = self.data['scan_number'][records[0],0]
      self.meta[groupIndex]['src_az'] = \
                             self.data['source_azel'][records[0]:records[-1],0]
      self.meta[groupIndex]['src_el'] = \
                             self.data['source_azel'][records[0]:records[-1],1]
      # what is self.data['source_long_lat']?
      self.meta[groupIndex]['R.A.'] = \
                            self.data['source_radec'][records[0]:records[-1],0]
      self.meta[groupIndex]['declination'] = \
                            self.data['source_radec'][records[0]:records[-1],1]
      if 'source' in self.meta:
        self.meta[groupIndex]['OBJECT'] = self.meta['source']
        self.logger.warning(
        "_unpack_spectra_._apply_metadata_: override of %s applied to 'source'",
                            self.meta['source'])
      else:
        self.meta[groupIndex]['OBJECT'] = \
                                        self.data['source_name'][records[0]][0]
      self.meta[groupIndex]['start'] = self.data['timestamp'][records[0],0]
      self.meta[groupIndex]['end'] = ( self.data['timestamp'][records[-1],0]+
                                       self.data['integ_time'][0,0])
      self.meta[groupIndex]['VELDEF'] = self.data['v_ref'][0,0]
      self.meta[groupIndex]['VELOCITY'] = self.data['vsys'][0,0]
      self.meta[groupIndex]['nod'] = nod #                        SUBREF_STATE?
      self.meta[groupIndex]['TAMBIENT'] = \
                                 self.data['weather'][records[0]:records[-1],0]
      self.meta[groupIndex]['PRESSURE'] = \
                                 self.data['weather'][records[0]:records[-1],1]
      self.meta[groupIndex]['HUMIDITY'] = \
                             self.data['weather'][records[0]:records[-1],2]/100
      self.meta[groupIndex]['WINDSPEE'] = \
                                 self.data['weather'][records[0]:records[-1],3]
      self.meta[groupIndex]['WINDDIRE'] = \
                                 self.data['weather'][records[0]:records[-1],4]
      if 'IFmode' in self.meta:
        self.meta[groupIndex]['IFmode'] = self.meta['IFmode']
        self.logger.warning(" override of %s applied to property 'IFmode'",
                            self.meta['IFmode'])
      else:
        self.meta[groupIndex]['IFmode'] = \
          IFmodes[int(        self.meta[groupIndex]['linefreq']/1e9
                      - round(self.meta[groupIndex]['linefreq']/2e9)*2 > 0)]
      self.logger.debug(" got metadata for %d'th scan", groupIndex)

    def _load_spectra_(groupIndex):
      """
      Pack the four spectra from an HDF5 record into the data matrix.
      The matrix axes are freq, ra, dec, pol, time, beam.
      """
      # AD HOC COMPUTING OF GROUPINDEX FROM GROUP AND VICE VERSA ASSUMES
      # THE FIRST GROUP/SCAN is 1.
      group = groupIndex+1
      self.spec[group] = zeros((32768,1,1,2,num_records,2))

      self.logger.debug(
                "_unpack_spectra_._load_spectra_: 'spec' group %d shape is %s",
                group, self.spec[group].shape)
      self.logger.debug("_unpack_spectra_._load_spectra_: records: %s",
                            records[groupIndex])
      for record in records[groupIndex]:
        index = records[groupIndex].index(record)
        self.logger.debug(
             "_unpack_spectra_._load_spectra_: doing %dth record of %dth scan",
                        index, groupIndex)
        self.logger.debug(
              "_unpack_spectra_._load_spectra_: copying spectra from record %d",
              record)
        self.logger.debug(
                         "_unpack_spectra_._load_spectra_: taking pol from %s",
                         self.meta[groupIndex]['pol'])
        self.logger.debug(
                        "_unpack_spectra_._load_spectra_: taking beam from %s",
                                                 self.meta[groupIndex]['BEAM'])
        # spectraCh1 has pol[0], beam[0]; second char has numeric value 1 or 2
        for specCh in [1,2,3,4]:
          specChIndex = specCh - 1
          pol_index  = int(self.meta[groupIndex]['pol'][specChIndex][1])-1
          beam_index = int(self.meta[groupIndex]['BEAM'][specChIndex][1])-1
          self.logger.debug("_unpack_spectra_._load_spectra_: storing IF%d"
                            +" spectra for %dth record, pol %s and beam %s"
                            +" for the %dth scan",
                            specCh, index, pol_index, beam_index, groupIndex)
          self.spec[group][:,0,0, pol_index, index, beam_index] = \
                             self.data['spectraCh'+str(specCh)].value[record,:]

    # The data array contains all the integrations for one scan.

    group = 0            # first group is 1
    nod = 1              # so the first nod will be 0
    num_records = 0      # to note the largest number of records in a scan
    records = {}         # the records comprising a scan

    # The next loop over 'num_recs' will fail if the file was terminates early,
    # that is, the number of records in the actual file is less than 'num_recs'
    self.spec = {}
    for index in range(self.num_recs):
      # 'scan_number' is the number of the scan to which the record belongs.
      # Each entry is a 1-item list in [1],[2],[3],[4],... or [-1] for no scan.
      # THIS MAKES THE SAFE ASSUMPTION THAT THERE IS ALWAYS AT LEAST ONE
      # type -1 record at the beginning of the file
      if self.scan_number[index] > 0:
        # Antenna on point
        if self.scan_number[index][0] == group:
          # Additional record of this group
          self.logger.debug("_unpack_spectra_: record %3d, scan %2d, nod %d",
                    index, int(self.scan_number[index][0]), nod)
          records[groupIndex].append(index)
        else:
          # First record of a group (scan)
          if group:
            # get the metadata for the group (scan)
            self.logger.debug("_unpack_spectra_: records for scan %d: %s",
                              group, records[groupIndex])
            records[groupIndex].append(index)
            # The order of the next three lines is important.  The first
            # creates the meta dict for the scan.  The second adds an item.
            _apply_metadata_(self, groupIndex, records[groupIndex])
            num_records = max(num_records, len(records[groupIndex]))
          self.logger.debug("_unpack_spectra_: record %3d, scan %2d, nod %d",
                    index, int(self.scan_number[index][0]), nod)
          # start a new group (scan)
          group += 1
          groupIndex = group-1
          records[groupIndex] = []
          self.meta[groupIndex] = {} # because the group numbers start with 1
          self.logger.debug(" Starting group (scan) %d", group)

          self.logger.debug(" spectra initialized for group %d", group)
          nod = 1-nod
    self.logger.debug(" have groups %s", list(self.spec.keys()))
    self.logger.debug(" have metadata for group indices %s", list(self.meta.keys()))
    self.logger.debug(" the largest number of records in a scan is %d", num_records)

    # normalize the last group
    _apply_metadata_(self, groupIndex, records[groupIndex])
    # this should not be needed

    # get spectra and optionally remove channel 0 and center channel
    self.logger.debug("_unpack_spectra_: Getting spectra for groups indices %s"
                      ,list(self.meta.keys()))
    self.logger.debug("_unpack_spectra_: Using records %s", records)
    for key in list(self.meta.keys()):
      # scan is the same as group number
      _load_spectra_(key)
      scan = key+1
      n_freqs, n_ra, n_dec, n_pols, n_recs, n_beams = self.spec[scan].shape
      for pol in range(n_pols):
        for beam in range(n_beams):
          for rec in range(n_recs):
            if self.clobber_0:
              self.spec[scan][0,0,0,pol,rec,beam] = \
                                            self.spec[scan][1,0,0,pol,rec,beam]
            if self.clobber_center:
              self.spec[scan][16384,0,0,pol,rec,beam] = \
                            (self.spec[scan][16383,0,0,pol,rec,beam] +
                             self.spec[scan][16385,0,0,pol,rec,beam])/2.
    # we don't need this anymore
    del(self.scan_number)
    del(self.num_recs)

  def _read_metadata_file_(self):
    """
    Some data are not included in the HDF5 file, like the state of the band
    switch, the source name and the sideband cabling.  These are in a separate
    file.
    """
    metaFile = self.filename.split('.')[0]+".meta"
    self.meta = {}
    if exists(metaFile):
      # for overrides metadata (header)
      fd = open(metaFile,'r')
      md = fd.readlines()
      fd.close()
      for line in md:
        if line.isspace():
          continue
        elif line[0] == "#":
          continue
        else:
          items = line.strip().split()
          if items[0].lower() == "band":
            self.meta["band"] = array(items[1].split(','), dtype=int)
          elif items[0].lower() == "ifmode":
            self.meta["IFmode"] = array(items[1].split(','))
          elif items[0].lower() == "source":
            self.meta["source"] = items[1]
          elif items[0].lower() == "opacity":
            self.opacity = float(items[1])
          elif items[0].lower() == "trec":
            self.Trec = float(items[1])

  def extract_valid_records(self):
    """
    Remove duplicate and empty records from each scan.

    Returns nparray of valid records, number of valid records for each pol,beam
    """
    valid_records = {}
    valid_times = {}
    newrec = {}
    for scannumber in list(self.spec.keys()):
      scan_index = scannumber-1
      num_recs = self.spec[scannumber].shape[4]
      valid_records[scannumber] = empty((self.spec[scannumber].shape))
      max_valid = 0 # largest number of valid records found, to size the array
      newrec[scannumber] = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
      for pol in [0,1]:
        for beam in [0,1]:
          newrec[scannumber][pol,beam] = 0
          for record in range(num_recs):
            if newrec[scannumber][pol,beam] == 0:
              if self.spec[scannumber][0,0,0,pol,record,beam]:
                # This is the first valid record; start new list of records
                rec = newrec[scannumber][pol,beam]
                valid_records[scannumber][:,0,0,pol,rec,beam] = \
                                   self.spec[scannumber][:,0,0,pol,record,beam]
                valid_times[scannumber] = self.meta[scan_index]['time']
                newrec[scannumber][pol,beam] += 1
              else:
                # try next record
                continue
            else:
              if self.spec[scannumber][0,0,0,pol,record,beam] != 0.0:
                # non -zero record
                rec = newrec[scannumber][pol,beam]
                if self.spec[scannumber][0,0,0,pol,record,beam] \
                   != valid_records[scannumber][0,0,0,pol,rec-1,beam]:
                  valid_records[scannumber][:,0,0,pol,rec,beam] = \
                    self.spec[scannumber][:,0,0,pol,record,beam]
                  valid_times[scannumber] = self.meta[scan_index]['time']
                  newrec[scannumber][pol,beam] += 1
                else:
                  # duplicate record
                  self.logger.debug(
                         "extract_valid_records: record %d of %d is duplicate",
                         record, num_recs)
              else:
                self.logger.debug(
                             "extract_valid_records: record %d of %d is empty",
                             record, num_recs)
          self.logger.info(
            "extract_valid_records: scan %d: %d records found for pol %d, beam %d",
                           scannumber, newrec[scannumber][pol,beam], pol, beam)
          max_valid = max(max_valid, newrec[scannumber][pol,beam])
      # Note that 'max_valid' is one higher than the highest index for valid
      # data because of the advance of the index after each record found.
      valid_records[scannumber] = valid_records[scannumber][:,:,:,:,:max_valid,:]
    return valid_records, valid_times, newrec

  def to_dataset(self):
    """
    Fills container with data in format suitable for pickling (saving to disk)

    This method is used by the calling program which initializes this class to
    put the data in a FITS-like format that can be put in a pickle file.
    """
    # TSYS is in here to conform to the DSN SDFITS convention with
    # shape (2,nrecs,2)
    static_keys = ['BANDWIDT','BEAM','EXPOSURE','IFmode','R.A.', 'VELDEF',
                   'VELOCITY','band','declination','end','linefreq','pol',
                   'start']
    self.container.header = {}
    metakeys = list(self.meta[0].keys())
    metakeys.sort()
    scanindexkeys = list(self.meta.keys()) # all scans have the same keys
    scanindexkeys.sort()
    for key in metakeys:
      if key in static_keys:
        self.logger.debug("to_dataset: %s has no keys", key)
        self.container.header[key] = self.meta[0][key]
      else:
        self.container.header[key] = {}
        self.logger.debug("to_dataset: %s has keys %s", key, list(self.meta.keys()))
        for index in scanindexkeys:
          self.container.header[key][index] = self.meta[index][key]
    self.container.data = self.uniq_recs # This puts in the spectral line data
    self.container.origin = self.filename
    self.container.header['num_records'] = self.num_uniq
    self.container.header["SITEELEV"] = float(self.data['SITEELEV'][0,0])
    self.container.header["SITELAT"]  = float(self.data['SITELAT'][0,0])
    self.container.header["SITELONG"] = float(self.data['SITELONG'][0,0])
    self.container.header["PROJID"]   = self.PROJID
    try:
      self.container.header['TAUZENIT'] = self.opacity
    except:
      pass
    return self.container
    
class SAOexaminer(object):
  """
  Class for examining SAOdataset objects
  """
  def __init__(self, project=None, dss=None, date=None):
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
      get_datasets          - get datasets from cPickle files
      get_obs_dirs          - get the directories used in this obs. session
      load_HDF5_files       - get datasets from HDF5 files
      load_ephemeris        - get the ephemeris for the designated source
      plot_Tsys             -
      plot_spectra          - plots avg_diff
      save_average          -
      smoothed              -
    
    @param datafiles : optional list of cPickle data files
    @type  datafiles : list of str
    """
    self.logger = logging.getLogger(logger.name+'.SAOexaminer')
    self.avg_diff = {}
    self.datafiles = None
    self.dataset = {}
    self.tau = {}      # cumulative integration time
    try:
      self.logger.debug("__init__: getting pickle files for %s at DSS-%s on %s",
                        project, dss, date)
      self.get_datasets(project=project, dss=dss, date=date)
    except TypeError as details:
      if re.search("NoneType", str(details)):
        self.logger.debug("__init__: getting HDF5 files instead")
        self.load_HDF5_files(project=project, dss=dss, date=date)
      else:
        self.logger.error("__init__: failed with TypeError: %s", details)
        raise RuntimeError("__init__: failed with TypeError: %s", details)
    except RuntimeError as details:
      if re.search("No data files found", str(details)):
        self.logger.debug("__init__: getting HDF5 files")
        self.load_HDF5_files(project=project, dss=dss, date=date)
      else:
        self.logger.error("__init__: failed with RuntimeError: %s", details)
        raise RuntimeError("__init__: failed with RuntimeError: %s", details)
    self.beam_switched_average()
    self.frame = None
    
  def get_obs_dirs(self, project=None, dss=None, date=None, rawfmt="HDF5"):
    """
    Gets the directories used for an observing session.
    
    If the session details are not provided, a menu-based inquiry is made
    
    The attribute rawdatapath is also computed but not returned.
    
    @param project : the project for which the observations were made
    @type  project : str
    
    @param dss : DSN station number
    @type  dss : int
    
    @param date : year and day of year as YYYY/DDD
    @type  data : str
    
    @param raw : get HDF5 data from RA_data store if True
    @type  raw : bool
    """
    self.datapath, self.projworkpath, self.rawdatapath = get_obs_dirs(
         *get_obs_session(dss=dss, project=project, date=date, datafmt=rawfmt))
    self.logger.debug("get_obs_dirs: found %s, %s, %s",
                            self.datapath, self.projworkpath, self.rawdatapath)
    return self.datapath, self.projworkpath, self.rawdatapath
  
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

  def get_datasets(self, project=None, dss=None, date=None):
    """
    Loads selected datasets and defines additional attributes.
    
    @param project : the project for which the observations were made
    @type  project : str
    
    @param dss : DSN station number
    @type  dss : int
    
    @param date : year and day of year as YYYY/DDD
    @type  data : str
    """
    # get the datafiles to be processed
    self.logger.debug("get_datasets: for %s at %s on %s", project, dss, date)
    dirs = self.get_obs_dirs(project=project, dss=dss, date=date)
    self.logger.debug("get_datasets: dirs are %s", dirs)
    self.datafiles = select_data_files(self.datapath)
    if self.datafiles:
      self.logger.debug("get_datasets: pickle files found")
      for datafile in self.datafiles:
        self.logger.info("get_datasets: processing %s", datafile)
        index = self.datafiles.index(datafile)
        self.savebasename = self.projworkpath+basename(datafile).split('.')[0]
        self.logger.info(" will save work as %s...", self.savebasename)
        # get the data
        fd = open(datafile, "rb")
        self.dataset[index] = cPickle.load(fd)
        fd.close()
        self.dataset[index].file    = datafile
        datapath_parts = datafile.split('/')
        self.dataset[index].project = datapath_parts[-5]
        self.dataset[index].dss     = int(datapath_parts[-4][-2:])
        self.dataset[index].year    = int(datapath_parts[-3])
        self.dataset[index].doy     = int(datapath_parts[-2])
        self.dataset[index].source  = self.dataset[index].header['OBJECT'][0]
      return len(self.dataset)
    else:
      self.logger.error("get_datasets: getting HDF5 files")
      return self.load_HDF5_files(project=project, dss=dss, date=date)

  def load_HDF5_files(self, project=None, dss=None, date=None):
    """
    """
    try:
      self.datafiles = select_data_files(self.rawdatapath, load_hdf=True)
    except AttributeError:
      self.get_obs_dirs(project=project, dss=dss, date=date)
      self.datafiles = select_data_files(self.rawdatapath, load_hdf=True)
    if self.datafiles:
      for datafile in self.datafiles:
        self.logger.info("load_HDF5_files: processing %s", datafile)
        index = self.datafiles.index(datafile)
        hdf = SAOhdf5(datafile)
        self.dataset[index] = hdf.to_dataset() # loads the data into the dataset
        mkdir_if_needed(self.datapath)
        self.dataset[index].save_pickle(self.datapath)
        self.dataset[index].file    = datafile
        datapath_parts = datafile.split('/')
        self.dataset[index].project = datapath_parts[-5]
        self.dataset[index].dss     = int(datapath_parts[-4][-2:])
        self.dataset[index].year    = int(datapath_parts[-3])
        self.dataset[index].doy     = int(datapath_parts[-2])
        self.dataset[index].source  = self.dataset[index].header['OBJECT'][0]
      return len(self.dataset)
    else:
      self.logger.error("load_HDF5_files: no files found in %s",
                        self.rawdatapath)
      raise RuntimeError("load_HDF5_files: no files found")
    
  def beam_switched_average(self, index=None, mode='LINE-PBSW'):
    """
    Returns pol average spectrum and total integration of given dataset(s)
    
    Also computes the average spectrum and total integration time for each pol
    and accumulates them in the SAOexaminer attributes avg_diff and tau.
    
    @param index : dataset index
    @type  index : int
    
    @param mode : observing mode (which aught to be in the header)
    @type  mode : str
    """
    if self.dataset == {}:
      self.get_datasets()
    if index == None:
      self.indices = list(self.dataset.keys())
    else:
      if type(index) == int:
        self.indices = [index]
      elif type(index) == list:
        self.indices = index
      else:
        self.logger.error("beam_switched_average: invalid indices: %s", index)
    first_scan = list(self.dataset[0].data.keys())[0]
    num_chan = self.dataset[0].data[first_scan].shape[0]
    for dindex in self.indices: # loop over datasets
      if mode == 'LINE-PBSW':
        spec, tau = self.dataset[dindex].average_calibrated_spectrum()
        if self.avg_diff == {}:
          # cumulative average
          # xan't do it in __init__ because num_chan is not known
          self.avg_diff = {0: zeros(num_chan), 1: zeros(num_chan)}
          self.tau = {0: 0.0, 1: 0.0}      # cumulative integration time
        for pol in [0,1]:
          self.avg_diff[pol] += spec[pol]
          self.tau[pol] += tau[pol]
      else:
        self.logger.error("get_datasets: mode %s not yet implemented", mode)
    for pol in [0,1]:
      self.avg_diff[pol] /= len(self.indices)
  
  def extract_window(self, xlimits, dsetidx=None, frame="RADI-LSR"):
    """
    Extracts a subset of a avg_diff spectrum
    """
    if self.avg_diff == {} or dsetidx:
      self.beam_switched_average(index=dsetidx, mode='LINE-PBSW')
    # we assume that all scans in all datasets have the same number of channels
    first_idx = self.indices[0]
    first_scan = list(self.dataset[first_idx].data.keys())[0]
    first_ds = self.dataset[first_idx]
    num_chan = first_ds.data[first_scan].shape[0]
    if frame:
      if frame == "RADI-LSR":
        x = ds.compute_X_axis("RADI-LSR",ds.dss)
      elif frame == "RADI-OBJ":
        if re.match('67P', first_ds.source):
          vspline = self.load_ephemeris('67P')
          x = first_ds.compute_X_axis("RADI-OBJ", first_ds.dss, vspline=vspline)
        else:
          self.logger.error("plot_spectra: no ephemeris for %s", self.object)
      else:
        x = first_ds.compute_X_axis(frame, first_ds.dss)
    else:
      frame = firs_ds.frame
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
    extract = {}
    for pol in list(self.avg_diff.keys()):
      extract[pol] = self.avg_diff[pol][ch1:ch2]
    d = Data(x[ch1:ch2], extract)
    d.frame = frame
    d.channels = (ch1,ch2)
    return d

#################################### Methods ##################################

def parse_filename(fname):
  """
  Parse an SAO .hf5 filename for date and source
  """
  parts = fname[7:-8].split('.')
  UNIXtime = float('.'.join(parts[:2]))
  sourcename = '.'.join(parts[2:])
  return UNIXtime, sourcename

