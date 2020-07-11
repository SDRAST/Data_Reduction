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
    
    
  def load_file(self, num_recs=5, datafile=None, schdfile=None, catlfile=None):
    """
    loads data from an OLR file
    """
    if datafile and schdfile:  #  and catlfile:
      pass
    else:
      self.logger.error("load_file: missing file: data=%s, sched=%s catlog=%s")
      raise RuntimeError("need data, a recording schedule, and coordinates")
    # set the datafile reader options
    options = {'toPrintHeaders':False, 'toPrintData': False,
               'lookForDate': True, 'fixTime': True}
    fname = self.sessionpath+datafile
    fmt = OLSR.checkFormat(fname)
    options.update({'format':fmt})
    # extract the datasets
    self.schedule = NP.genfromtxt(schdfile, dtype=[('time',     'S24'),
                                                   ('duration', int),
                                                   ('name',     'S8')])
    self.times = self.schedule['time']
    self.durations = self.schedule['duration']
    num_samples = len(self.times)
    self.records = {}
    if True:
      for start in self.times[:num_recs]:
        idx = list(self.times).index(start)
        options.update({'startDate': datetime.datetime.strptime(
                                                          start.decode("utf-8"),
                                                          '%Y/%m/%d/%H:%M:%S')})
        options.update({'duration': float(self.durations[idx])})
        if fmt == "RDEF":
          self.records[idx] = OLSR.RDEF(fname, options)
        elif fmt == "VDR":
          self.records[idx] = OLSR.VDR(fname, options)
        elif fmt == "SFDU":
          self.records[idx] = OLSR.SFDU(fname, options)
        else:
          logger.error("get_file_metadata: %s has unknown format %s", 
                   filename, fmt)
    
class Map(Observation):
  """
  """
  def __init__(self):
    """
    """
    pass

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
  options = {'toPrintHeaders': False, 'format': fmt}
  for fname in files:
    fidx = files.index(fname)
    fmt = OLSR.checkFormat(fname)
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
                   +header[filename]['TIMETAG_PICOSECONDS_OF_THE_SECOND']/1e12),
                      "year":   header[filename]['TIME_TAG_YEAR'],
                      "DOY":    header[filename]['TIME_TAG_DOY']}
  return header, metadata
  
def print_file_metadata(project, dss, year, doy, pattern):
  """
  """
  header, metadata = get_file_metadata(project, dss, year, doy, pattern)
  print(27*" "+" Freq."+4*" "+"BW   File")
  output = []
  for fidx in list(metadata.keys()):
    output.append("%24s %7.1f %5.1f %38s"
                                     % (time.ctime(metadata[fidx]['unixtime']),
                                        metadata[fidx]['freq'],
                                        metadata[fidx]['bw'],
                                        os.path.basename(metadata[fidx]['file'])
                                       ))
  output.sort()
  return output
  
    
