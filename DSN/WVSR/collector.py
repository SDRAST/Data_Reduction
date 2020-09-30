"""
provides class WVSRmetadataCollector to get data from logs
"""
import glob
import logging
import re

from datetime import datetime
from os.path import basename, splitext

import MonitorControl.Configurations.projects as projects
from MonitorControl.Configurations.projects.sources import get_all_source_data
from Data_Reduction.DSN.WVSR import parse_scan_files
from Data_Reduction.DSN.WVSR.SpecData import get_channel_IDs
from DatesTimes import VSR_script_time_to_timestamp
#from MonitorControl.Configurations.DSN_standard import standard_equipment
from support.lists import unique


logger = logging.getLogger(__name__)

class WVSRmetadataCollector:
  """
  Class to assemble data for a session data file
  
  Public attributes::
    chan_names   - list of channel IDs
    datadir      - session directory in 'project_data' directory
    date         - datetime.datetime for observation date
    doy          - day of year of observation
    dss          - station where observations were done
    equip        - telescope equipment description
    logger       - logging.Logger object
    project      - name of AUTO project (not real name)
    real_obs_dir - project's path to session directory
    scaninfo     - metadata for every scan
    scankeys     - sorted list of scan numbers
    sourcedata   - details of sources and calibrators
    wvsrdir      - location of WVSR files
    wvsrlogs     - list of log file names
    year         - year of observation
    
  'scaninfo' example::
    In [12]: collector.scaninfo[6]
    Out[12]: 
    {'subch 1': '16-237-001-s0006_d14_RCP_LCP.wvsr.fft',
     'subch 2': '16-237-002-s0006_d14_RCP_LCP.wvsr.fft',
     'end': datetime.datetime(2016, 8, 24, 9, 42),
     'source': 'w5w-fregg51-ref',
     'start': datetime.datetime(2016, 8, 24, 9, 40, 30)}
     
  Methods::
    get_WVSR_names(self):
    get_metadata_from_WVSR_logs(self):
    parse_WVSR_log(self, logname):
    get_metadata_from_FFT_logs(self):
    parse_fft_log_name(self, logname):
    parse_fft_log(self, lines):
    get_metadata_from_WVSR_scripts(self, band='X', pol="RCP"):
    
  """
  def __init__(self, activity, dss, year, doy, time):
    """
    initiate a WVSR configuration description collector
    
    Mainly this creates two dicts::
      wvsr_cfg - configuration for one WVSR
      fft_meta - details on the processing all the scans
      
    @param project : AUTO project name
    @type  project : str
    
    @param dss : DSN station
    @type  dss : int
    
    @type year : int
    @type doy  : int
    
    @param time : start time without colon: HHMM
    @type  time : str
    """
    self.activity = activity
    self.project = projects.activity_project(activity)
    self.dss = dss
    self.year = year
    self.doy = doy
    self.time = time
    self.logger = logging.getLogger(logger.name+".WVSRmetadataCollector")
    self.logger.info("__init__: WVSRmetadataCollector initializing")
    # get all the directories involved
    auto_obsdir, self.real_obsdir, self.activity_dir, self.project_dir, \
        self.datadir, self.wvsrdir, self.fftdir = \
        projects.get_session_dirs( activity, dss, year, doy)
    if self.wvsrdir:
      pass
    else:
      self.logger.error("__init__: no WVSR data directory")
      raise RuntimeError("no WVSR data directory")
    # get the high-level metadata
    self.get_WVSR_names()
    self.get_metadata_from_WVSR_logs()
    # get the parameters for each FFT configuration
    self.scaninfo = {} # defined for each scan
    self.get_metadata_from_FFT_logs()
    # specify station equipment for each IF
    #self.equip = {}
    for wvsr in self.wvsrnames:
      #self.equip[wvsr] = {}
      self.logger.debug('__init__: %s config: %s', wvsr, self.wvsr_cfg[wvsr])
      for IF in self.wvsr_cfg[wvsr]['channels']:
        self.logger.debug('__init__: %s IF %s config: %s',
                          wvsr, IF, self.wvsr_cfg[wvsr][IF])
        band = self.wvsr_cfg[wvsr][IF]['IF_source'].split('_')[1]
        #self.equip[wvsr][IF] = standard_equipment(dss, band)
        # add source name to 'scaninfo'
        self.get_metadata_from_WVSR_scripts(band=band,
                                            pol=self.wvsr_cfg[wvsr][IF]['pol'])
    
  def get_WVSR_names(self):
    """
    find the recorders used
    
    creates a dict of WVSR identifiers based on simple name like 'wvsr2'::
      In [17]: collector.wvsrnames
      Out[17]: {'wvsr2': 'WVSR2_10.16-237-084500'}
    """
    wvsrlogs = glob.glob(self.wvsrdir+"/WVSR*")
    if wvsrlogs:
      pass
    else:
      self.logger.error("get_WVSR_names: no WVSR logs found")
      raise RuntimeError("get_WVSR_names: no WVSR logs found")
    self.logger.debug("get_WVSR_names: logs: %s", wvsrlogs)
    self.wvsrnames = {}
    self.wvsrlogs = []
    for f in wvsrlogs:
      logname = basename(f)
      if re.search(self.time, logname):
        self.wvsrlogs.append(f)
        self.wvsrnames[logname[:5].lower()] = logname
    self.logger.debug("get_WVSR_names: %s", self.wvsrnames)
  
  def get_metadata_from_WVSR_logs(self):
    """
    parse the WVSR FFT logs for metadata
    
    Each WVSR typically has two digital IF channels with orthogonal pols.
    'wvsr_cfg' is a dict.  Example::
      In [3]: collector.wvsr_cfg
      Out[3]: 
      {'wvsr2': {1: {'DSP': 'dsp2',
                     'IF_source': '14_X_RCP',
                     'chan_id 1': {'bandwidth': 8000000,
                                   'bits': 8,
                                   'sfro': 209386000},
                     'chan_id 2': {'bandwidth': 8000000,
                                   'bits': 8,
                                   'sfro': 484824000},
                     'expid': '16_237_DSS-14_ARCP',
                     'pol': 'RCP',
                     'rf_to_if_lo': 8100,
                     'subchannels': ['chan_id 1', 'chan_id 2']},
                 2: {'DSP': 'dsp1',
                     'IF_source': '14_X_LCP',
                     'chan_id 1': {'bandwidth': 8000000, 'bits': 8,
                                   'sfro': 209386000},
                     'chan_id 2': {'bandwidth': 8000000, 'bits': 8,
                                   'sfro': 484824000},
                     'expid': '16_237_DSS-14_ALCP',
                     'pol': 'LCP',
                     'rf_to_if_lo': 8100,
                     'subchannels': ['chan_id 1', 'chan_id 2']},
                 'channels': [1, 2],
                 'from log': 'WVSR2_10.16-237-084500',
                 'user': 'cjn,'}}
    """
    self.wvsr_cfg = {}
    for logname in self.wvsrlogs:
      wvsr = basename(logname)[:5].lower()
      self.logger.debug("get_metadata_from_WVSR: WVSR ID: %s", wvsr)
      self.parse_WVSR_log(logname)
      self.logger.debug("get_metadata_from_WVSR: from log %s",
                                               self.wvsr_cfg[wvsr]['from log'])
      self.logger.debug("get_metadata_from_WVSR: user: %s",
                                                   self.wvsr_cfg[wvsr]['user'])
      self.wvsr_cfg[wvsr]["channels"] = []
      for key in list(self.wvsr_cfg[wvsr].keys()):
        if type(key) == int:
          self.wvsr_cfg[wvsr]["channels"].append(key)
          self.wvsr_cfg[wvsr][key]['subchannels'] = []
          for subkey in list(self.wvsr_cfg[wvsr][key].keys()):
            if subkey[:7] == 'chan_id':
              self.wvsr_cfg[wvsr][key]['subchannels'].append(subkey)
          self.wvsr_cfg[wvsr][key]['subchannels'].sort()
      self.logger.debug("get_metadata_from_WVSR: %s", self.wvsr_cfg)

  def parse_WVSR_log(self, logname):
    """
    Extracts metadata from the WVSR log
    
    Most metadata are not consider time sensitive so that the last value found
    is the value that is returned.  The exception is attenuator data which may
    change during a track.
    
    In quite a few cases we simply exec() a statement this is in the line or
    is constructed from parts of the line.  A typical line looks like this:
    16/237 08:45:15 wvsr2 ATT[2]: att = 15, des_amp = -10, cur_amp = -10.07, \
                                                     max_amp = 0, min_amp = -50
    """
    logfile = open(logname,'r')
    logtext = logfile.read()
    logfile.close()
    lines = logtext.split('\n')
    self.logger.debug("parse_WVSR_log: read %s", basename(logname))
    EXPID = {}
    RF_TO_IF_LO = {}
    IFS = {}    # for IFs from IFS[*] ... CRITICAL
    altIFS = {} # for IFs from IFS[*] ... PROGRESS
    CHAN = {}
    ATT = {}
    # parse the rest of log file
    for line in lines[1:]:
      if re.search("EXPID\\[[12]\\]", line): # experiment ID
        # parse a log file line like::
        #   16/237 08:45:13 wvsr2 EXPID[1]: EVT 301 PROGRESS: \
        #              exp 16_237_DSS-14_ARCP created by Client 3625 cjn 4-2053
        if line.find("Welcome") > -1:
          continue
        # take the variable name from the 4th item leaving off the colon
        # and add item to EXPID dict
        parts = line.split()
        self.logger.debug("parse_WVSR_log: EXPID line parts: %s", parts)
        wvsrID = parts[2] # third item
        user = parts[-2]  # second last item
        exec(parts[3][:-1]+"= "+"'"+str(parts[8])+"'")
        # this is the clumsier way to get the signal for the IF
        subparts = parts[8].split('_')
        altIFS[int(parts[3][6])] = subparts[0]+"_"+subparts[3][:-3]+"_"+subparts[3][-3:]
        self.logger.debug("parse_WVSR_log: IFS from EXPID is %s", IFS)
      if re.search('RF_TO_IF_LO\[',line): # receiver or first LO
        # parse a log file line like::
        #   16/237 08:45:15 wvsr2 RF_TO_IF_LO[2]: RF_TO_IF_LO: value = 8100
        # add to RF_TO_IF_LO for the IF channel to dict
        parts = line.split()
        exec(parts[3][:-1]+"= "+str(parts[-1]))
      if re.search("CHAN\\[.\\]: dsp", line): # digital signal processor module
        # parse a lof file line like::
        #   16/237 08:45:17 wvsr2 \
        #          CHAN[2]: dsp1:1 chan_id = 001, bandwidth = 8000000, bits = 8
        # add subchannel and its parameters to CHAN dict
        parts = line.split()
        chan_str = parts[3][:-1]
        subchan_str = "'"+parts[5]+" ' + str(int("+parts[7][:-1]+"))"
        self.logger.debug("parse_WVSR_log: for CHAN doing %s with %s",
                     chan_str, subchan_str)
        # create subchannel item in CHAN if needed
        if CHAN == {}:
          exec(chan_str+"={}") # creates a dict for each CHAN
        try:
          eval(chan_str+".has_key("+subchan_str+")")
        except:
          exec(chan_str+"={}")
          exec(chan_str+"["+subchan_str+"]={}")
        else:
          exec(chan_str+"["+subchan_str+"]={}")
        # add assigned DSP
        exec(chan_str+"['DSP']='"+parts[4].split(":")[0]+"'")
        # add DDCLO, add bandwidth and bits/sample
        exec(chan_str+"["+subchan_str+"]['"+parts[8]+"']="
             +str(int(parts[10][:-1])))
        exec(chan_str+"["+subchan_str+"]['"+parts[11]+"']="
             +str(int(parts[13])))
      if re.search("rsp_ddclo", line): # subchannel offset LO
        # create subchannel offset
        parts = line.split()
        self.logger.debug("parse_WVSR_log: line parts: %s", parts)
        dsp, subchan = parts[8].split(':')
        offset = int(float(parts[9])) # don't need fractional Hz
        # assign to correct channel
        for key in list(CHAN.keys()):
          if CHAN[key]['DSP'] == dsp:
            if "chan_id "+subchan in CHAN[key]:
              CHAN[key]["chan_id "+subchan]["sfro"] = offset
            else:
              self.logger.info(
                  "parse_WVSR_log: no subchannel ID %s for IF channel %s",
                     subchan,key)
      if re.search("IFS\\[.\\]:.*PROGRESS", line): # IF signal source for
        # IF switch inputs
        self.logger.debug("parse_WVSR_log: IFS from line: %s", line)
        parts = line.split()
        self.logger.debug("parse_WVSR_log: line parts: %s", parts)
        # exclude the colon from the signal source
        exec("alt"+parts[3][:-1]+"= '"+str(parts[10][:-1])+"'")
        self.logger.debug("parse_WVSR_log: IFS from PROGRESS is %s", IFS)
      if re.search("IFS\\[.\\]:.*CRITICAL", line): # better IF signal source
        parts = line.split()
        exec(parts[3][:-1]+"= '"+str(parts[10][:-1])+"'")
      if re.search("ATT\\[.\\]:.*cur_amp", line): # attenuator and power
        self.logger.debug("parse_WVSR_log: ATT line parts: %s", parts)
        parts = line.split()
        UT = WVSR_script_time_to_timestamp(*parts[:2])
        try:
          # this executes something like" "ATT[2]['att']=15"
          self.logger.debug("parse_WVSR_log: trying %s",
                            parts[3][:-1] + \
                            "["+str(UT) + "]['" + parts[4] + "']" + \
                            "="+parts[6].strip(','))
          exec(parts[3][:-1] + "["+str(UT) + "]['" + parts[4] + "']" + \
               "="+parts[6].strip(','))
        except KeyError as details:
          # ATT[2] is not defines as a dict
          self.logger.debug("parse_WVSR_log: failed because %s", details)
          # this executes something like:'ATT[2]={1472028315: {}}'
          self.logger.debug("parse_WVSR_log: trying %s",
                            parts[3][:-1]+"={"+str(UT)+": {}}")
          exec(parts[3][:-1]+"={"+str(UT)+": {}}")
          # this executes something like "ATT[2][1472028315]['att']=15"
          self.logger.debug("parse_WVSR_log: trying %s",
                            parts[3][:-1] + "["   +str(UT)+ "]['"+parts[4] + \
                            "']"+"="+parts[6].strip(','))
          exec(parts[3][:-1] + "[" + str(UT) + "]['" + parts[4] + "']" + \
               "="+parts[6].strip(','))
          # now that the dict exists add the power
          exec(parts[3][:-1] + "[" + str(UT) + "]['" + parts[10] + "']" + \
               "="+parts[12].strip(',')) 
    self.logger.info("parse_WVSR_log: EXPID = %s", EXPID)
    self.logger.info("parse_WVSR_log: RF_TO_IF_LO = %s", RF_TO_IF_LO)
    # fill in missing IF signals
    for key in altIFS:
      if key in IFS:
        continue
      else:
        IFS[key] = altIFS[key]
    self.logger.info("parse_WVSR_log: IFS = %s", IFS)
    self.logger.info("parse_WVSR_log: CHAN = %s", CHAN)
    self.logger.info("parse_WVSR_log: ATT = %s", ATT)
    # add to WVSR configuration
    self.wvsr_cfg[wvsrID] = {'user': user,
                             'from log': basename(logname)}
    # now build the dict
    #   CHAN keys will also apply to ATT
    for key in list(CHAN.keys()):
      self.wvsr_cfg[wvsrID][key] = {}
      self.wvsr_cfg[wvsrID][key]['pol'] = EXPID[key][-3:]
      self.wvsr_cfg[wvsrID][key]['expid'] = EXPID[key]
      self.wvsr_cfg[wvsrID][key]['rf_to_if_lo'] = RF_TO_IF_LO[key]
      if RF_TO_IF_LO[key] == 1325:
        band = 'L'
      elif RF_TO_IF_LO[key] == 2000:
        band = 'S'
      elif RF_TO_IF_LO[key] == 8100:
        band = 'X'
      else:
        raise RuntimeError('unknown LO frequency', RF_TO_IF_LO[key])
      try:
        self.wvsr_cfg[wvsrID][key]['IF_source'] = IFS[key]
      except:
        # hope the EXPID is sensible
        self.wvsr_cfg[wvsrID][key]['IF_source'] = EXPID[key][-7:].replace(
                                                       EXPID[key][-4],band+'_')          
      subchans = list(CHAN[key].keys())
      for subch in subchans:
        self.wvsr_cfg[wvsrID][key][subch] = {}
        if subch[:7] == 'chan_id':
          for param in list(CHAN[key][subch].keys()):
            self.wvsr_cfg[wvsrID][key][subch][param] = CHAN[key][subch][param]
        else:
          self.wvsr_cfg[wvsrID][key][subch] = CHAN[key][subch]
      self.wvsr_cfg[wvsrID][key]['attenuation'] = {}
      for time in list(ATT[key].keys()):
        self.wvsr_cfg[wvsrID][key]['attenuation'][time] = {}
        for param in list(ATT[key][time].keys()):
          self.wvsr_cfg[wvsrID][key]['attenuation'][time][param] = ATT[key][time][param]
    self.logger.debug("parse_WVSR_log: %s results: %s", wvsrID,
                       self.wvsr_cfg[wvsrID])

  def get_metadata_from_FFT_logs(self):
    """"
    Gets metadata about the FFT processing of WVSR files
    
    This must be invoked after 'get_metadata_from_WVSR_logs'.
    
    This assumes that two IFs are combined into one signal with full
    polarization information.  Then the subchannels are the same for both.
    
    The problem with this right now is that there is no clean way to handle
    simultaneous WVSRs so this is going to have to assume that 'wvsrnames' has
    only one item for now.
    
    Returns a dict like::
      In [2]: collector.fft_meta[3]
      Out[2]: 
      {1: {'chan_id 1': {'datafile':  '/data2/16g237/16-237-001-s0003_d14_RCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90},
           'chan_id 2': {'datafile':  '/data2/16g237/16-237-002-s0003_d14_RCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90}},
       2: {'chan_id 1': {'datafile':  '/data2/16g237/16-237-001-s0003_d14_LCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90},
           'chan_id 2': {'datafile':  '/data2/16g237/16-237-002-s0003_d14_LCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90},
           'n_freqs':   131072,
           'n_samples': 8000000}}
    
    Notes
    =====
    Log File Names
    --------------
    Somewhere along the way Chuck changed filenames from::
      16-364-001-s0001_d14_RCP_LCP.wvsr.log
    to::
      16-364-001-s0001_d14_EGG_RCP_LCP.wvsr.log
    """
    fftlogs = glob.glob(self.fftdir+"*.log")
    fftlogs.sort()
    if fftlogs == []:
      raise RuntimeError("There are no logs in %s", self.fftdir)
    scans = []
    subchls = []
    # get the scan numbers
    for log in fftlogs:
      # parse file name
      logname = splitext(basename(log))[0]
      self.logger.debug("get_metadata_from_FFT_logs: log basename: %s",
                        logname)
      YR, DOY, scan, subch, projcode, band, dss  = \
                                               self.parse_fft_log_name(logname)
      self.logger.debug("get_metadata_from_FFT_logs: parsed: %s",
                        (YR, DOY, scan, subch, projcode, band))
      scans.append(scan)
      subchls.append(subch)
    scannums = unique(scans)
    self.logger.debug("get_metadata_from_FFT_logs: scannums: %s", scannums)
    subchannels = unique(subchls)
    self.logger.debug("get_metadata_from_FFT_logs: subchannels: %s",
                      subchannels)
    # now parse the FFT log files
    self.fft_meta = {}
    for scan in scannums:
      self.fft_meta[scan] = {} # first level dict on scan
      self.scaninfo[scan] = {}
      # create a dict for each IF
      for channel in self.wvsr_cfg[list(self.wvsrnames.keys())[0]]["channels"]:
        self.fft_meta[scan][channel] = {} # second level dict on IF channel
      # each IF has the same number of subchannels
      for subch in subchannels:
        channelID = 'subch '+str(subch)
        # parse the FFT log
        if projcode and band:
          logname = "%2d-%03d-%03d-s%04d_d%2d_%3s_%s_RCP_LCP.wvsr.log" % \
              (self.year-2000, self.doy, subch, scan, dss, projcode, band)
        elif projcode:
          logname = "%2d-%03d-%03d-s%04d_d%2d_%3s_RCP_LCP.wvsr.log" % \
                  (self.year-2000, self.doy, subch, scan, dss, projcode)
        elif dss:
          logname = "%2d-%03d-%03d-s%04d_d%2d_RCP_LCP.wvsr.log" % \
                  (self.year-2000, self.doy, subch, scan, dss)
        else:
          logname = "%2d-%03d-%03d-s%04d_RCP_LCP.wvsr.log" % \
                  (self.year-2000, self.doy, subch, scan)
        try:
          logfile = open(self.fftdir+logname)
        except IOError as details:
          self.logger.debug("get_metadata_from_FFT_logs: no %s",
                            self.fftdir+logname)
          break
        fftfile = logname.replace('.log','.fft')
        self.scaninfo[scan][channelID] = fftfile
        lines = logfile.readlines()
        logfile.close()
        subch_key = "chan_id "+str(subch)
        self.logger.debug("get_metadata_from_FFT_logs: for %s", subch_key)
        extracted = self.parse_fft_log(lines)
        self.logger.debug("get_metadata_from_FFT_logs: extracted %s",
                          extracted)
        if extracted:
          # the start and end time are the same for both channels so use 1
          try:
            self.scaninfo[scan]['start'] = datetime.strptime(
                                                 extracted[1]['first rec'],
                                                 "%Y %j:%H:%M:%S.%f")
          except ValueError as details:
            self.logger.error(
                 "get_metadata_from_FFT_logs: scan %d: %s", scan, str(details))
          try:
            self.scaninfo[scan]['end'] = \
                                datetime.strptime(extracted[1]['last rec'],
                                                   "%Y %j:%H:%M:%S.%f")
          except ValueError as details:
            self.logger.error(
                 "get_metadata_from_FFT_logs: scan %d: %s", scan, str(details))
        if 'start' in self.scaninfo[scan] and \
                                            'end' in self.scaninfo[scan]:
          # move some data to fft_meta structure
          extracted_keys = list(extracted.keys())
          IFs = []
          for key in extracted_keys:
            if type(key) == int:
              # integers are IF number; collect for loop over IFs
              IFs.append(key)
            else:
              # str are other dict items
              self.fft_meta[scan][key] = extracted[key]
          # create a dict for each IF
          for ch in IFs:
            self.fft_meta[scan][ch][subch_key] = {} # third level dict on subch
            for key in list(extracted[ch].keys()): # these keys are variable names
              self.fft_meta[scan][ch][subch_key][key] = extracted[ch][key]
      # did we get good time data from either subchannel?
      if 'start' in self.scaninfo[scan] and \
                                            'end' in self.scaninfo[scan]:
        self.logger.debug("get_metadata_from_FFT_logs: scan %d info: %s",
                          scan, self.scaninfo[scan])
        self.logger.debug("get_metadata_from_FFT_logs: scan %d FFT metadata: %s",
                        scan, self.fft_meta[scan])
        pass
      else:
        # we can't use this scan
        self.logger.debug(
          "get_metadata_from_FFT_logs: could not get all metadata for scan %d", 
          scan)
        self.scaninfo.pop(scan)
    if self.scaninfo == {}:
      raise RuntimeError("There are no FFT log files with valid data")
    
  def parse_fft_log_name(self, logname):
    """
    Gets metadata from FFT log file name
    
    The earliest form of the log file name is::
      YY-DOY-sub-sSCAN_RCP_LCP.wvsr.fft
    where 'sub' is the 3-digit sub-channel number and 'SCAN' is the 4-digit
    scan number.
      
    A later form of the filename is::
      YY-DOY-sub-sSCAN_dSS_RCP_LCP.wvsr.log
    where 'SS' is the 2-digit station number
    
    Some of these were later changed to::
      YY-DOY-sub-sSCAN_dSS_PRJ_RCP_LCP.wvsr.log
    where PRJ is a 3-digit project code. This change did not affect the year,
    DOY, channel, or scan number.
    
    Later still the filenames became::
      YY-DOY-sub-sSCAN_dSS_PRJ_?band_RCP_LCP.wvsr.log
    
    Eventually the project code was changed to a four digit activity code.
    
    Returns year, DOY, scan, subch, activity, band, dss
    """
    self.logger.debug("parse_fft_log_name: %s", logname)
    # split on the hyphens to get year, DOY and subchannel
    name_parts = logname.split('-')
    # extract year
    YR = int(name_parts[0])
    if 2000+YR != self.year:
      self.logger.error("parse_fft_log_name: %s is for wrong year",
                        logname)
      return False
    # extract DOY
    DOY = int(name_parts[1])
    if DOY != self.doy:
      self.logger.error("parse_fft_log_name: %s is for wrong DOY",
                        logname)
      return False
    # extract subchannel
    subch = int(name_parts[2])
    # split the underscore delimited part (which starts with scan number)
    lastparts = name_parts[3].split("_")
    # extract scan, the first item in the last parts list
    scan = int(lastparts[0][1:]) # removes 's'
    dss = int(lastparts[1][1:])  # removes 'd'
    if dss != self.dss:
      self.logger.error("parse_fft_log_name: %s is for wrong DSS",
                          logname)
      return False
    elif len(lastparts) == 4:
      return YR, DOY, scan, subch, None, None, None
    # at this point the filename takes different forms
    if len(lastparts) == 5:
      # project code in log name
      return YR, DOY, scan, subch, lastparts[2], None, dss
    elif len(lastparts)  == 6:
      # band also in log name
      # ['s0001', 'd14', 'EGG', 'Xband', 'RCP', 'LCP.wvsr']
      return YR, DOY, scan, subch, lastparts[2], lastparts[3], dss
    else:
      self.logger.debug("parse_fft_log_name: last parts: %s", lastparts)
      self.logger.error("parse_fft_log_name: unmanaged completion")
      
      
  def parse_fft_log(self, lines):
    """
    Get metadata from one scan's log file.
        
    Each log has a preamble section with details about the processing followed
    by reports for each FFT.
    
    The output is a dict like this example for scan 1::
      In [2]: collector.fft_meta[1]
      Out[2]: 
      {1: {'datafile': '/data2/16g237/16-237-002-s0001_d14_RCP.wvsr',
           'first rec': '2016 237:09:30:31.000',
           'last_rec': '2016 237:09:32:00.000',
           'n_recs': 91,
           'n_secs': 91},
       2: {'datafile': '/data2/16g237/16-237-002-s0001_d14_LCP.wvsr',
           'first rec': '2016 237:09:30:31.000',
           'last_rec': '2016 237:09:32:00.000',
           'n_recs': 91,
           'n_secs': 91},
      'n_freqs': 131072,
      'n_samples': 8000000}
      where the first key is the IF channel.
            
    We don't know yet how to handle FFT logs if there are two or more WVSRs so
    'wvsr_cfg' for one and only one has a similar structure with WVSR ID in
    place of scan number.
    """
    in_preamble = True
    first_part = True
    found = {}
    for line in lines:
      parts = line.strip().split()
      if parts == []:
        continue
      if in_preamble:
        if line[:10] == "Input file": # this starts a channel section
          datafile = parts[2].lstrip('<').rstrip('>')
          VSRtype = line[-1]
        if parts[0] == "Channel" and len(parts) > 2:
          channel = int(parts[1][:-1]) # strip off colon
          found[channel] = {}
          found[channel]["datafile"] = datafile
        if re.search("First", line):
          found[channel]["first rec"] = ' '.join(parts[2:])
        if re.search("Last", line):
          found[channel]["last rec"] = ' '.join(parts[2:])
        if line[:4] == "nrec":
          if int(parts[2]):
            # non-zero
            found[channel]["n_recs"] = int(parts[2])
            found[channel]["n_secs"] = int(parts[5])
          else:
            # no data recorded
            return {}
        if parts[0] == 'Unequal':
          in_preamble = False
      else:
        if re.search("ns.*npts.*power_two", line):
          found["n_samples"] = int(parts[2][:-1])
          found["n_freqs"] = int(parts[5])
          break
    return found

  def get_metadata_from_WVSR_scripts(self, band='X', pol="RCP"):
    """
    @param band : S, X or K
    @type  band : str
      
    @param pol : RCP or LCP (or R or L)
    @type  pol : str
    """
    if len(pol) == 1:
      pol += "CP"
    template = self.activity_dir+"*"+band+pol+".scr"
    self.logger.debug("get_metadata_from_WVSR_scripts: template: %s",
                      template)
    scripts = glob.glob(template)
    if len(scripts) == 0:
      # try an older template
      template = self.activity_dir+"*"+pol+".scr"
      self.logger.debug("get_metadata_from_WVSR_scripts: trying %s",
                      template)
      scripts = glob.glob(template)
    if len(scripts) == 0:
      self.logger.error("get_metadata_from_WVSR_scripts: no scripts found")
      raise RuntimeError("no WVSR scripts found")
    scripts.sort()
    script = scripts[-1] # use the last version of the script
    self.logger.debug("get_metadata_from_WVSR_script: script is %s", script)
    scriptfile = open(script, 'r')
    lines = scriptfile.readlines()
    scriptfile.close()
    # get source for each scan
    for line in lines:
      line = line.strip()
      if line[:5] == 'srcid':
        self.logger.debug("get_metadata_from_WVSR_scripts: processing '%s'",
                          line)
        parts =  line.split()
        scanID, source = parts[1].split(':')
        scan = int(scanID)
        if scan in self.scaninfo:
          self.scaninfo[scan]['source'] = source
    
      

