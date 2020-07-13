"""
Data_Reduction.Malargue
=======================

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
import Data_Reduction.DSN as DSN
import Data_Reduction.DSN.OLSR as OLSR
import DatesTimes as DT

logger =  logging.getLogger(__name__)

class Observation(DSN.Observation):
  """
  Class for observations based on open-loop recording made as DSA-3 (DSS-84)
  
  The arguments for the superclass initialization are::
    ``parent`` (typically ``self`` or ``None``),
    ``name``   (will be set to a YEAR/DOY default if not provided),
    ``dss``    (required),
    ``date``   (YEAR/DOY required),
    ``start``  (optional, usually taken from file header), 
    ``end``    (otional, inferred from header and number of records), and
    ``project`` (required).
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
    DSN.Observation.__init__(self, parent=parent, name=name, dss=dss, date=date,
                                  project=project)
    self.logger =  mylogger
    
  def parse_obs_file(self, filename=None):
    """
    An ``.obs`` file has some preliminary observation followed by a section for
    each ``scan``, which might be a map or series of spectra::
    
      #=========================================================================================================
      # SCAN_NUM  SRC_ID            START_TIME         STOP_TIME                  RA         DEC      TFREQ
      # --------  ----------------  -----------------  -----------------  ----------  ----------  ----------
      S 001       ...               2020-163T11:44:46  2020-163T11:54:49  999.000000  999.000000  31950000000.000000
      # DATAFILE                               COH_FLAG  DOR_MULT                 FSUB    HARMONIC
      # -------------------------------------  --------  -----------  ----------------  ----------
      D NET4n001tSsMG12rOPc01-20163114446.prd  F         xxx              -3.00000e+06           0
      D NET4n001tSsMG12rOPc02-20163114446.prd  F         xxx              -1.00000e+06           0
      D NET4n001tSsMG12rOPc03-20163114446.prd  F         xxx              +1.00000e+06           0
      D NET4n001tSsMG12rOPc04-20163114446.prd  F         xxx              +3.00000e+06           0
      D NET4n001tSsMG12rOPc05-20163114446.prd  F         xxx              -3.00000e+06           0
      D NET4n001tSsMG12rOPc06-20163114446.prd  F         xxx              -1.00000e+06           0
      D NET4n001tSsMG12rOPc07-20163114446.prd  F         xxx              +1.00000e+06           0
      D NET4n001tSsMG12rOPc08-20163114446.prd  F         xxx              +3.00000e+06           0
      Z
      #=========================================================================================================
    """
    if filename == None:
      # search for it
      self.logger.debug("parse_obs_file: looking in %s", self.sessionpath)
      filenames = glob.glob(self.sessionpath+"*.obs")
      self.logger.debug("parse_obs_file: found %s", filenames)
      if len(filenames) == 1:
        pass
      elif len(filenames) == 0:
        self.logger.error("No observation file found in %s", self.sessionpath)
        raise RuntimeError("Data reduction requires an observation file")
      else:
        self.logger.warning("Found %d observation files; using %s",
                            len(filenames), filename)
      filename = filenames[0]
    fd = open(filename, "rt")
    lines = fd.readlines()
    scan = {}
    for line in lines:
      parts = line.strip().split()
      if line[0] == "#":
        pass
      elif line[0] == "S":
        skip_this_scan = False
        scan_ID = int(parts[1])
        if float(parts[7]) == 0.0:
          # no useful data
          skip_this_scan = True
          continue
        else:
          scan[scan_ID] = {"start":  DT.ISOtime2datetime(parts[3]),
                           "stop":   DT.ISOtime2datetime(parts[4]),
                           "ra2000": float(parts[5]),
                           "dec000": float(parts[6]),
                           "freq":   float(parts[7])/1e6}
        scan[scan_ID]['channel'] = {}
      elif line[0] == "D" and not skip_this_scan:
        filename = parts[1]
        ch_ptr1 = filename.index('c')
        ch_ptr2 = filename.index('-')
        chan = int(filename[ch_ptr1+1:ch_ptr2])
        scan[scan_ID]['channel'][chan] = {'file': filename}
    return scan
        
        
        
class Map(Observation):
  """
  """
  def __init__(self):
    """
    """
    pass

