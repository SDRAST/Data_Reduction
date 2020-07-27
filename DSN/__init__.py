"""
Data_Reduction.DSN
==================

Subclasses for reducing data taken with DSN-like open loop recorders.

Open-loop recorders are raw IF voltage recorders that are not synchronized with
the communications between the spacecraft and the ground station.  As such, they
are the most basic kind of recorder possible in radio astronomy, equivalent to
VLBI recorders.  Indeed, an early implementation was known as the "VLBI Science
Recorder" (VSR), followed by later varieties of VSR and eventually, the OSR.

OLR recordings at different stations are indeed combined for VLBI measurements
of spacecraft with respect to distant radio sources, a powerful navigation tool.

Raw IF recordings can be computational converted into any of the standard
signal types used in radio astronomy -- square-law detected power, spectra, 
Stokes parameters, VLBI U-V maps, high time ans spectral resolution pulsar data,
*etc.*
"""
import datetime
import glob
import logging
import numpy as NP
import os.path
import time

import Data_Reduction as DR
import Data_Reduction.DSN.OLSR as OLSR
import DatesTimes as DT

logger =  logging.getLogger(__name__)

class Observation(DR.Observation):
  """
  Class for observations based on open-loop recordings
  
  The arguments for the superclass initialization are ``parent`` (typically 
  ``self`` or ``None``), ``name`` (will be set to a YEAR/DOY default if not 
  provided), ``dss`` (required), ``date`` (YEAR/DOY required), ``start``, 
  ``end``, and ``project`` (required).
    """
  def __init__(self, parent=None, name=None, dss=None,
                     date=None, start=None, end=None,
                     project=None):
    """
    The project, station, and date are needed to locate the directory for the
    working files.
  
    The header item ``STATION_ID`` is not correct in files from Malargue.
    The start time is in the header as ``TIME_TAG_YEAR``, ``TIME_TAG_DOY``,
    ``TIME_TAG_SECOND_OF_DAY`` and ``TIMETAG_PICOSECONDS_OF_THE_SECOND``.
  
    These channel metadata can be extracted from the header: ``SAMPLE_RATE`` for
    bandwidth, ``RF_TO_IF_DOWNCONV`` for the band and receiver center frequency,
    ``IF_TO_CHANNEL_DOWNCONV`` for the channel center frequency.
    
  Args:
    parent (Session): (optional) session to which this observation belongs
    name (str):       (optional) an identifier; default is station ID + "obs"
    dss (int):        (required) station number
    date (str):       (required) "YEAR/DOY"
    start (float):    (optional) UNIX time at the start
    end (float):      (optional) UNIX time at the end
    project (str):    (required) directory under /usr/local/projects
    
    """
    mylogger = logging.getLogger(logger.name+".Observation")
    DR.Observation.__init__(self, parent=parent, name=name, dss=dss, date=date,
                                  project=project)
    self.logger =  mylogger
    
class Map(Observation):
  """
  """
  def __init__(self):
    """
    """
    pass
  

class Recording(DR.Recording):
  """
  Superclass for DSN-based VSR, WVSR, OLR, etc. data files.
  
  This also covers data files made at DSA-3 TTCP recorder
  """
  headerFormats = {
    'RDEF': {'RECORD_LABEL':                     {'Size':4, 'Type':'CHARACTER'},
             'RECORD_LENGTH':             {'Size':4, 'Type':'UNSIGNED INTEGER'},
             'RECORD_VERSION_ID':         {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'STATION_ID':                {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'SPACECRAFT_ID':             {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'SAMPLE_SIZE':               {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'SAMPLE_RATE':               {'Size':4, 'Type':'UNSIGNED INTEGER'},
             'VALIDITY_FLAG':             {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'AGENCY_FLAG':               {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'RF_TO_IF_DOWNCONV':           {'Size':8, 'Type':'FLOATING POINT'},
             'IF_TO_CHANNEL_DOWNCONV':      {'Size':8, 'Type':'FLOATING POINT'},
             'TIME_TAG_YEAR':             {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'TIME_TAG_DOY':              {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'TIME_TAG_SECOND_OF_DAY':
                                          {'Size':4, 'Type':'UNSIGNED INTEGER'},
             'TIMETAG_PICOSECONDS_OF_THE_SECOND':
                                            {'Size':8, 'Type':'FLOATING POINT'},
             'CHANNEL_ACCUMULATED_PHASE':   {'Size':8, 'Type':'FLOATING POINT'},
             'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT0':{'Size':8,
                                                       'Type':'FLOATING POINT'},
             'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT1':{'Size':8,
                                                       'Type':'FLOATING POINT'},
             'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT2':{'Size':8,
                                                       'Type':'FLOATING POINT'},
             'CHANNEL_PHASE_POLYNOMIAL_COEFFICIENT3':{'Size':8,
                                                       'Type':'FLOATING POINT'},
             'EMPTY_FIELDS1':                       {'Size':36, 'Type':'EMPTY'},
             'PREDICT_PASS_NUMBER':       {'Size':2, 'Type':'UNSIGNED INTEGER'},
             'UPLINK_BAND':               {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'DOWNLINK_BAND':             {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'TRACK_MODE':                {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'UPLINK_DSS_ID':             {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'OLR_ID':                    {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'OLR_SOFTWARE_VERSION':      {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'CHANNEL_POWER_CALIBRATION_FACTOR':
                                            {'Size':4, 'Type':'FLOATING POINT'},
             'TOTAL_FREQUENCY_OFFSET':      {'Size':8, 'Type':'FLOATING POINT'},
             'CHANNEL_NUMBER':            {'Size':1, 'Type':'UNSIGNED INTEGER'},
             'EMPTY_FIELDS2':                       {'Size':19, 'Type':'EMPTY'},
             'END_LABEL':                 {'Size':4, 'Type':'UNSIGNED INTEGER'}}
             }
  headerSize = {"RDEF": 176}
            
  def __init__(self, session=None, name=None, filename=None,
                     dss=None, date=None, project=None):
    """
    initiate a RS data set
    
    Args
    ====
      channels: (dict) associates files with names
      parent:   (Session) optional
      name:     (str) optional; defaults to station name + "rec"
      dss:      (int) required
      date:     (str) required as "YEAR/DOY"
      project:  (str) required; must be known project name
    
    Notes
    =====
    The project, station and date are used to find the path to the data file.
    """
    mylogger = logging.getLogger(logger.name+".Recording")
    if filename:
      pass
    elif session:
      mylogger.warning("A filename will be constructed from the session")
    else:
      raise RuntimeError("A session or a filename must be given")
    DR.Recording.__init__(self, session=session, name=name)
    self.logger = logging.getLogger(logger.name+".Recording")
    
  def not_needed(self):
    try:
      self.format = self.checkFormat()
    except FileNotFoundError:
      self.logger.error("check file path: %s", self.file)
      raise FileNotFoundError
    except NameError:
      self.logger.error('%s is not a recognized format', self.format)
      raise NameError
    if self.format == "RDEF":
      pass
    else:
      self.logger.error('not currently coded for format %s', self.format)
      raise RuntimeError('invalid file format')
    self.header = {}

  def checkFormat(self):
    """
    check the file format
    """
    # Open file, read first 4 bytes as string
    try:
      self.f = open(self.session_dir+self.filname, 'rb')
    except FileNotFoundError:
      raise FileNotFoundError("could not find %s" % self.file)
      self.format = None
      return self.format
    dataBytes = np.fromfile(self.f, dtype=np.uint8, count=4)
    data = "".join(map(chr, dataBytes))
    self.f.close()
    # If string matches a particular format, return format
    if data == 'NJPL':
      return 'SFDU'
    elif data == 'VSRD':
      return 'VDR'
    elif data == 'RDEF':
      return data
    else:
      return None

  def readHeader(self, record=0):
    """
    get the header
    """
    self.f = open(self.sessionpath+self.file, 'rb')
    self.f.seek(0,0) # default file pointer position is start of file
    if record:
        self.f.seek(self.header['RECORD_LENGTH']*record, 0)
        
    headerBytes = np.fromfile(self.f, 
                              dtype=np.uint8, 
                              count=self.headerSize[self.format])
    self.logger.debug("readHeader: got %d bytes", len(headerBytes))
    self.f.close()

    # Check for errors
    if len(headerBytes) != self.headerSize[self.format]:
      self.header = {'bytes': self.headerSize[self.format]}
      self.header['EOF'] = True
      self.logger.warning("readHeader: end of file")
      return False
    else:
      self.header = {}
    
    startByte = 0
    for header in self.headerFormats[self.format]:
        self.logger.debug("readHeader: header is %s", header)
        # Extract data bytes
        endByte = startByte + self.headerFormats[self.format][header]['Size']
        dataBytes = headerBytes[startByte:endByte]
        self.logger.debug("readHeader: dataBytes: %s", dataBytes)
        # Convert
        self.logger.debug("readHeader: header type is %s",
                                self.headerFormats[self.format][header]['Type'])
        if self.headerFormats[self.format][header]['Type'] == 'CHARACTER':
            data = "".join(map(chr, dataBytes))
        elif self.headerFormats[self.format][header]['Type'] == 'UNSIGNED INTEGER':
            data = self.data2int(dataBytes)
        elif self.headerFormats[self.format][header]['Type'] == 'FLOATING POINT':
            data = self.data2float(dataBytes)

        # Save data
        if self.headerFormats[self.format][header]['Type'] != 'EMPTY':
            self.header[header] = data

        # Increment start byte
        startByte = endByte
    return self.header

  def loadOLSRfile(self, filename):
    """
    loads from a DSS open loop science recorder
    
    Currently type RDEF is supported.
    """    
    self.file = open(filename)
    

# --------------------------- module functions ---------------------------

def get_file_metadata(project, dss, year, doy, pattern):
  """
  """
  projdatapath, sessionpath, rawdatapath = \
                          DR.get_obs_dirs(project, dss, year, doy, datafmt=None)
  files = glob.glob(sessionpath+pattern)
  logger.info("get_file_metadata: files: %s", files)
  header = {}
  metadata = {}
  options = {'toPrintHeaders': False}
  for fname in files:
    fidx = files.index(fname)
    fmt = OLSR.checkFormat(fname)
    if fidx:
      pass
    else:
      options.update({'format': fmt})
    filename = os.path.basename(fname)
    # get the first header
    if fmt == "RDEF":
      header[filename] = OLSR.readHeaders(fname, options)
    elif fmt == "VDR":
      header.append(OLSR.VDR(fname, options))
    elif fmt == "SFDU":
      header.append(OLSR.SFDU(fname, options))
    else:
      logger.error("get_file_metadata: %s has unknown format %s", 
                   filename, fmt)
    # get the metadata
    metadata[fidx] = {"file":   fname,
                      "bpS":    header[filename]["SAMPLE_SIZE"],
                      "bw":     header[filename]["SAMPLE_RATE"]/1e6,
                      "f_ofst":(header[filename]['RF_TO_IF_DOWNCONV']-31950000000\
                       +header[filename]['IF_TO_CHANNEL_DOWNCONV'])/1e6,
          "freq":  (31950000000+header[filename]['IF_TO_CHANNEL_DOWNCONV'])/1e6,
                      "unixtime": DT.VSR_tuple_to_timestamp(
                                              header[filename]['TIME_TAG_YEAR'], 
                                               header[filename]['TIME_TAG_DOY'], 
                                     header[filename]['TIME_TAG_SECOND_OF_DAY']
                  +header[filename]['TIMETAG_PICOSECONDS_OF_THE_SECOND']//1e12),
                      "year":   header[filename]['TIME_TAG_YEAR'],
                      "DOY":    header[filename]['TIME_TAG_DOY']}
    metadata[fidx]['size'] = os.path.getsize(fname)
    Bps = header[filename]["SAMPLE_RATE"]*header[filename]["SAMPLE_SIZE"]*2/8
    metadata[fidx]['duration'] = metadata[fidx]['size']/Bps
  return header, metadata
  
def print_file_metadata(project, dss, year, doy, pattern):
  """
  """
  header, metadata = get_file_metadata(project, dss, year, doy, pattern)
  print(27*" "+" Freq."+4*" "+"BW   File"+32*" "+"duration")
  output = []
  for fidx in list(metadata.keys()):
    output.append("%24s %7.1f %5.1f %38s %4d"
                                     % (time.ctime(metadata[fidx]['unixtime']),
                                        metadata[fidx]['freq'],
                                        metadata[fidx]['bw'],
                                        os.path.basename(metadata[fidx]['file']),
                                        metadata[fidx]['duration']
                                       ))
  output.sort()
  return output
  
    
