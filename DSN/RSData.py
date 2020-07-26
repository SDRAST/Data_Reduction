"""
RSData - Radio Science Data

Data handler for such recorders as VSR, WVSR, OLR, etc. It recognizes the
'RDEF', 'VDR', and 'SFDU' formats.

This uses:
 
  fileObject.seek(offset[, whence])

  Parameters:
  
    offset − This is the position of the read/write pointer within the file.
    whence − This is optional and::
    
      * defaults to 0 which means absolute file positioning;
      * 1 means seek relative to the current position and
      * 2 means seek relative to the file's end.

"""
import logging
import numpy as np

import Data_Reduction.DSN as DSN
import DatesTimes as DT

logger = logging.getLogger(__name__)

class Observation(DSN.Observation):
  """
  class to read and convert data in radio science recorder data files
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
            
            
  def __init__(self, channels, parent=None, name=None, dss=None,
                     date=None, start=None, end=None,
                     project=None):
    """
    initiate a RS data set
    
    Args
    ====
      channels: (dict) associates files with names
      parent:   (Session) optional
      name:     (str) optional; defaults to station name + "obs"
      dss:      (int) required
      date:     (str) required as "YEAR/DOY"
      start:    (datetime.datetime) optional
      end:      (datetime.datetime) optional
      project:  (str) required; must be known project name
    
    Notes
    =====
    The project, station and date are used to find the path to the data file.
    """
    self.logger = logging.getLogger(logger.name+".RSData")
    if name:
      self.name = name
    else:
      self.name = None
    DSN.Observation.__init__(self, parent=parent, name=name, dss=dss,
                     date=date, start=start, end=end,
                     project=project)
    for key,value in channels:
      self.channel = self.Channel
    self.file = filename
    self.curRecord = 0
    # set the ``f`` and ``format`` attributes
    try:
      self.format = self.checkFormat()
    except FileNotFoundError:
      self.logger.error("check file path: %s", self.file)
    except NameError:
      self.logger.error('%s is not a recognized format:', self.format)
  
  def finish(self):
    self.f.close()
    self.logger.info("finish: done")
  exit = finish
  quit = finish
  
  def checkFormat(self):
    """
    check the file format
    """
    # Open file, read first 4 bytes as string
    try:
      self.f = open(self.sessionpath+self.file, 'rb')
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

  def data2int(self, data):
    """
    data sample to integer
    
    The digitizer only tests for > or < some level. This makes a 1s-complement
    notation preferred.  The -0 value is avoided. >0 is represented as 00000000
    and <0 by 11111111.  So then::

      In [11]: bin2sint('00000000') 
      Out[11]: 1
      In [12]: bin2sint('11111111')
      Out[12]: -1

      In [13]: bin2sint('00000001')
      Out[13]: 3
      In [14]: bin2sint('11111110') 
      Out[14]: -3

      In [15]: bin2sint('00000011')
      Out[15]: 7
      In [16]: bin2sint('11111100')
      Out[16]: -7
      
      In [17]: bin2sint('10000000') 
      Out[17]: -255
      In [18]: bin2sint('01111111')
      Out[18]: 255
    """
    dec = 0
    if self.format == 'RDEF':
        for power in range(len(data)-1, -1, -1):
            dec = dec + data[power]*(2**(8*power))
    elif self.format == 'VDR' or self.format == 'SFDU':
        for power in range(0, len(data)):
            dec = dec + data[power]*2**(8*(len(data)-power-1.0))

    dec = int(dec)
    return dec

  def data2float(self, data):
    """
    data sample to float
    """
    self.logger.debug("data2float: data = %s", data)
    # Create binary
    binVec = '';
    if self.format == 'RDEF':
        for index in range(len(data)-1, -1, -1):
            binStr = "{0:b}".format(data[index])
            binStr = binStr.zfill(8)
            binVec = binVec + binStr
        self.logger.debug("data2float: binVec = %s", binVec)
    elif self.format == 'VDR' or self.format == 'SFDU':
        for index in range(0, len(data), 1):
            #binStr = "{0:b}".format(data[index])
            binStr = np.binary_repr(data[index], width=8)
            binVec = binVec + binStr
            
    if len(data) == 8:
        # 64-bit Float
        sign = int(binVec[0], 2)
        expo = int(binVec[1:12], 2) - 2**10 + 1
        frac = 1 + int(binVec[12:], 2)/float(2**52)
        value = ((-1)**sign) * (2**expo) * frac
    elif len(data) == 4:
        # 32-bit Float
        sign = int(binVec[0], 2)
        expo = int(binVec[1:9], 2) - 2**7 + 1
        frac = 1 + int(binVec[9:], 2)/float(2**24)
        value = ((-1)**sign) * (2**expo) * frac
    else:
        raise NameError('Floating point number must contain 32 or 64 bits!')
    return value

  def timeOfRecord(self):
    """
    convert time in header to people format
    
    Note that the default way of displaying a datetime object is without the
    fractional seconds, and to the nearest integer microsecond if using
    ``strftime()``. So anything less than 500,000 ps will be seen as 0.
    """
    return DT.VSR_to_datetime((self.header['TIME_TAG_YEAR'],
                                 self.header['TIME_TAG_DOY'],
                                 self.header['TIME_TAG_SECOND_OF_DAY']
                       +self.header['TIMETAG_PICOSECONDS_OF_THE_SECOND']/1.e12))

  # records should be private classes
  def readRecordData(self, record=0):
    """
    read and convert data to array of complex samples
    
    The array returned will have records [startRecord:endRecord-1] in the usual
    Python way of defining ranges.
    """
    self.f = open(self.sessionpath+self.file, 'rb')
    # position to start of record
    self.f.seek(self.header['RECORD_LENGTH']*record, 0)
    # read the header
    headerBytes = np.fromfile(self.f, 
                              dtype=np.uint8, 
                              count=self.headerSize[self.format])
    # read the data
    dataBytes = self.self.headerSize[self.format]
    
    if len(dataBytes) == 0:
      self.logger.warning("readRecordData: no data read")
    else:
      # Get samples per word, etc
      samplesPerWord = 16/headers['SAMPLE_SIZE']
      if samplesPerWord in [2,4,8,16]:
        pass
      else:
        logger.error('RDEF: SAMPLE SIZE %d is not currently supported!',
                     headers['SAMPLE_SIZE'])
        return None
    
    
