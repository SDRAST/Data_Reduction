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
import glob
import logging
import os.path

import Data_Reduction as DR
import Data_Reduction.DSN.OLSR as OLSR

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
                                  start=start, end=end, project=project)
    self.logger =  mylogger
    

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
  for fname in files:
    fmt = OLSR.checkFormat(fname)
    filename = os.path.basename(fname)
    # get the first header
    if fmt == "RDEF":
      header[filename] = (OLSR.readHeaders(fname,
                                           {'toPrintHeaders': False, 
                                            'format': fmt}))
    elif fmt == "VDR":
      header.append(OLSR.VDR(fname, {'toPrintHeaders':True, 'format': fmt}))
    elif fmt == "SFDU":
      header.append(OLSR.SFDU(fname, {'toPrintHeaders': True, 'format': fmt}))
    else:
      logger.error("get_file_metadata: %s has unknown format %s", 
                   filename, fmt)
  return header
