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
import Data_Reduction.DSN.RSData as RSData
import DatesTimes as DT

logger =  logging.getLogger(__name__)

class Observation(DSN.RSData.Observation):
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
    
  def load_file(self, num_recs=5, datafile=None, schdfile=None, catlfile=None):
    """
    loads data from an OLR file
    
    This is Malargue soecific because of the catalog and schedule file formats.
    
    We load five records at a time or poor old Python really bogs down.
    
    I think we proceed as follows:
    
    #. Get a list of data files in the directory.
    #. We parse the ``*.obs`` file in the directory to get:
    
       #. a list of scan numbers,
       #. the base file name for each scan (without scan or channel number),
       #. a list of channel numbers from the first scan.
       
    #. For each scan (which is a map):
    
       #. Ignore if there is no corresponding data file (``NET*.prd``)
       #. Find the corresponding schedule (``sch*.txt``) and load the times and 
          position names.
       #. Open the matching catalog (``cat*.txt``) and for each position get the
          coordinates.
       #. for each channel:
       
          #. Form the file name: ``NET4n%03dtSsMG12rOPc%02d*.prd`` where the 
             first formatted item is scan number and the second is channel.
          #. Process the ordered list of data files found with the above mask.
          
             #. Read and parse the file header.
             #. For each record (in groups of five records at a time) for 
                efficiency):
          
                #. For the first record of the file only:
                   #. Put the time and coordinates in the first row of the 
                      structured numpy array.
                #. Read the record data from the datafile.
                #. Process the data (e.g. sum square average, power spectrum, 
                   spectrogram, etc.)
                #. Save the processed data for each record in the numpy dict 
                   keyed on channel number.
             
       #. save the numpy array of reduced data.
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
    fmt = self.checkFormat(fname)
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
          self.records[idx] = self.RDEF(fname, options)
        elif fmt == "VDR":
          self.records[idx] = self.VDR(fname, options)
        elif fmt == "SFDU":
          self.records[idx] = self.SFDU(fname, options)
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


class Session(DR.Session):
  """
  Class for a Malargue observing session on a given year and DOY
  
  Attributes Inherited from Superclass::
    doy (int)               - day of year for session
    logger (logging.Logger) - logging.Logger object
    parent (object)         - a data reduction session (mult. observ. sessions)
    year (int)              -
    doy (int)               -
    project (str)           - 
    session_dir (str)       - path to results from this session
  
  A session usually refers to a telescope, date and project.  This will
  normally define a path to the session directory.
  """
  def __init__(self, parent=None, date=None, project=None, dss=None,
                     path=None):
    """"
    initialize data reduction for one observing session
    
    Args
    ====
      parent:  (object) optional class for a data reduction tool
      date:    (int) required
      project: (str) required
      dss      (int) required
      path     (str) optional
      
    If `path` is given for a non-standard observing files location, and it does
    not exist, it will be created.  Then the Recording and Observation instances
    must be directed to where the files are.
    """
    mylogger = logging.getLogger(logger.name+".Session")
    DR.Session.__init__(self, parent=parent, date=date, project=project,
                              dss=dss, path=None)
    self.logger = mylogger
    metafiles = self.select_data_files(name_pattern="*.obs") # automatically
      
                     
class Recording(DSN.Recording):
  """
  Metadata and directory for raw data files from DSA-3.
  """
  def __init__(self, session=None, name=None, filename=None,
                     dss=84, date=None, project=None):
    """
    Initialize a Recording object
    
    Args
    ====
      session:  (Session) optional
      name:     (str) optional; will be constructed if not given
      filename: (str) optional: will glob for `*.obs` if not given
      dss:      (int)
      date:     (str) `"YEAR/DOY"`
      project:  (str) `"SolarPatrol"`
    """
    DSN.Recording.__init__(self, session=session, name=name, filename=filename,
                                 dss=dss, date=date, project=project)
    self.metafiles = self.session.select_data_files(name_pattern="*.obs")
    self.scans = {}
    for mfile in self.metafiles:
      self.scans[os.path.basename(mfile)] = self.parse_obs_file(filename=mfile)
    
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


