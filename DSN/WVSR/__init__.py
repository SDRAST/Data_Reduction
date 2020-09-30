"""
Processes .scans files for metadata and get data from the .scr file

.scans files are used to write WVSR script files.  They coordinate antenna
motion with WVSR recordings and assign scan numbers to WVSR recordings

The relevant commands in the .scr file are::
  DDCLO       - selects the reatime signal processor (RSP) downconversion
                filter. Allowed values are AUTO, 80, 160, 240, 320, 400, 480,
                and 560. AUTO specifies that the predicts in current use are to
                be used to set the LO so that the channel filters are close to 
                centered in the 160 MHz pass band.
  FROV        - FRequency OVerride allows the user to override the sky
                frequency predicts with a fixed frequency. If the value is 0,
                then the predicts will be used instead.
  RF_TO_IF_LO - the RF to IF LO used to calculate baseband frequencies.
                If value is left at 0 the software sets it based on the 
                DOWNLINK_BAND from predicts and the following values::
                  L = 1200,
                  S = 2000,
                  X = 8100,
                  Ka = 31700 MHz 
  SFRO        - sets the specified WVSR channel frequency offset to the
                specified frequency value. The offset is relative to the
                predicts carrier plus any FRO frequency offset that has been
                specified, including any rate. The values  must be real with 
                optional K, M or G appended for kHz, MHz or GHz respectively.

"""
import logging
import re

from datetime import datetime
from os.path import basename, exists, splitext

from DatesTimes import VSR_script_time_to_timestamp

wvsr_dir = "/data/"   # was data2 on crab14
postproc_dir = wvsr_dir+"post_processing/"
wvsr_fft_dir = postproc_dir+"auto/" # was cjnaudet/auto on crab14

logger = logging.getLogger(__name__)

def get_last_scr_file(obsdir, pol):
  """
  """
  logger.debug("get_last_scr_file: looking in %s", obsdir)
  files = glob.glob(obsdir+"*"+pol+".scr")
  logger.debug("get_last_scr_file: found %s", files)
  if files == []:
    return None
  else:
    files.sort()
    return files[-1]  

def get_frequency(text):
  """
  """
  if text[-1].isdigit():
    freq = float(text[:-1])
  elif text[-1] == "M":
    freq = float(text[:-1])*1e6
  else:
    raise RuntimeError("Unknown FROV suffix %s", text[-1])
  return freq
  
def get_WVSR_parameters(scrdata):
  """
  get WVSR data from script file.
  
  Obsolete!  Should come from log files
  """  
  metadata = {}
  bandwidth = {}
  bits_per_sample = {}
  sfro = {}
  for line in scrdata:
    parts = line.strip().split()
    if parts == []:
      continue
    if parts[0].upper() == "RF_TO_IF_LO":
      # the RF to IF LO used to calculate baseband frequencies.
      rf_to_if_lo = parts[-1]
      metadata['rf_to_if_lo'] = rf_to_if_lo
      continue
    elif parts[0].upper() == "FROV":
      # sky frequency to be used instead of what is in the predicts
      freq_override = get_frequency(parts[-1])
      metada
      continue
    elif parts[0].upper() == "SFRO":
      # sets the specified WVSR channel frequency offset to the specified
      # frequency value.
      chan = int(parts[1])
      sfro[chan] = get_frequency(parts[-1])
      continue
    elif parts[0].upper() == "CHAN":
      # Allocates a channel to this client and configures its bandwidths and
      # bitrates 
      if parts[1].upper() == "ALL":
        continue
      chan = int(parts[1])
      bandwidth[chan] = get_frequency(parts[2])
      bits_per_sample[chan] = int(parts[3])
    elif parts[0].upper() == "DDCLO":
      chan = int(parts[1])
    elif parts[0].upper() == "EXPID":
      pass
  ch_freq = {}
  if int(rf_to_if_lo)*1e6 != freq_override:
    logger.warning("get_WVSR_parameters: RF_TO_IF_LO does not match FROV")
  for chan in list(sfro.keys()):
    ch_freq[chan] = freq_override + sfro[chan]
  return ch_freq, bandwidth, bits_per_sample

def get_obs_session(obsdir):
  """
  Makes a name of the form YYdDOY where::
    YY  - two-digit year
    g   - lowercase Complex code, e.g., g for GDSCC
    DOY - must be three digits
  """
  parts = obsdir.split('/')
  dss  = int(parts[-4][-2:])
  year = int(parts[-3])
  doy  = int(parts[-2])
  project = parts[4]
  return project, dss, year, doy

def make_datadir_name(obsdir):
  """
  """
  project, dss, year, doy = get_obs_session(obsdir)
  src = "%02d" % (year - 2000)
  if dss == 14:
    src += "g"
  elif dss == 43:
    src += "c"
  elif dss == 63:
    src += "m"
  else:
    raise Exception("%s is not a valid 70-m antenna", dss)
  src += "%03d/" % doy
  return src
  
  
def get_Naudet_FFT_dir(obsdir):
  """
  Parses the observation directory name to find the FFT files
  
  The FFT files are in /data2/cjnaudet/auto/. The name convention for the
  directory with the FFT files is YYsDOY, where YY is the two-digit year, s is
  the lowercase first letter of the Complex name followed by the DOY.
  
  @param obsdir : directory for the observation session
  @type  obsdir : str
  
  @return: str (full path to FFT directory)
  """
  # get the observing session details
  project, dss, year, doy = get_obs_session(obsdir)
  
  # relative links are compatible with sshfs
  return wvsr_fft_dir + make_datadir_name(obsdir)

def get_latest_scan_files(obsdir):
  """
  returns the latest versions of the .scans files
  
  This gets the .scans file(s) from the observation directory. It then breaks
  down the file names to get the year, DOY and pol.  It constructs a template
  without the version number.  It then globs with this template to get all the
  versions and takes the latest version. In then reads the latest version.
  
  @param obsdir : directory for the observing session
  @type  obsdir : str
  """
  scanfiles = glob.glob(obsdir+"*.scans")
  logger.debug("get_latest_scan_files: scanfiles: %s", scanfiles)
  latest_files = []
  for filename in scanfiles:
    # get date and pol information from the file name
    parts = basename(filename).split("_")
    filenamebase = parts[0]
    year = 2000+int(filenamebase[:2])
    doy = int(filenamebase[2:])
    pol = parts[-1].split('.')[0]
    pol_template = obsdir + "_".join(parts[:-2]) + "*_"+pol+".scans"
    logger.debug("get_latest_scan_files: looking for %s", pol_template)
    files = glob.glob(pol_template)
    # take the highest version number
    files.sort()
    lastfile = files[-1]
    logger.debug("parse_latest_scan_files: getting %s", lastfile)
    latest_files.append(lastfile)
  return latest_files
  
def parse_scan_files(obsdir):
  """
  Get the metadata from all the scans from the *.scans files
  
  The names of the .wvsr files are used to determine the channel and
  polarization for each 
  
  The datafiles are in a large data area which has been linked to subdir FFT in
  the observation directory. Association between scans, the sources observed
  and the datafiles is put into the dict returned by this function.  The dict
  keys are::
    fftdir - where the FFT files are
    scanno - scan number (int) which is a dict. Each scan dict has::
             start  - datetime.datetime object
             end    - datetime.datetime object
             source - str
             files  - the names of the FFT files
  
  NOTES
  =====
  Redundancy
  ----------
  Note that while there are two scans files for each observation, only one
  needs to be parsed because the contents are the same for each polarization.
  
  Data file path
  --------------
  The filenames in the .scan files may have a path prepended if the file is not
  initially in the usual place. However, this is corrected after the session
  and so just strip off the path.
  
  @param obsdir :
  @type  obsdir : str
  
  @return: dict with item for each IF channel
  """
  # get observation details
  logger.debug("parse_scan_files: checking %s", obsdir)
  parts = obsdir.split('/')
  year = int(parts[-3])
  DOY = int(parts[-2])
  project = parts[4]  
  
  FFTdir = get_Naudet_FFT_dir(obsdir) # obsdir+"FFT/"
  logger.debug("parse_scan_files: FFT files are in %s", FFTdir)

  scans = {}
  scanfiles = get_latest_scan_files(obsdir)
  if scanfiles:
    logger.debug("parse_scan_files: found %s", scanfiles)
    # scan file contents are the same except for pol
    scanfile = scanfiles[0]
    logger.debug("parse_scan_files: examining %s", scanfile)
    fd = open(scanfile,'r')
    scan_meta = fd.readlines()
    logger.debug("parse_scan_files: %d lines of metadata", len(scan_meta))
    fd.close()
    # Extract metadata from the scan  and fft file names.
    #   Example base filename: '16208_1045_014000_AUTO_V00_LCP.scans'
    #                           YYDDD HHMM DSS    PROJ VER POL
    # Polarization for this file
    #pol = splitext(basename(scanfile))[0].split('_')[-1]
    #logger.debug("parse_scan_files: processing %s polarization", pol)
    for line in scan_meta:
      logger.debug("parse_scan_files: processing %s", line.strip())
      lineparts = line.split()
      logger.debug("parse_scan_files: parts: %s", lineparts)
      # extract the scan number
      scannum = int(lineparts[0])
      scans[scannum] = {}
      # extract the times and source
      startUXtime = VSR_script_time_to_timestamp(year,lineparts[1])
      endUXtime = VSR_script_time_to_timestamp(year,lineparts[2])
      logger.debug("parse_scan_files: scan %d from %s to %s", scannum,
                        startUXtime, endUXtime)
      scans[scannum]['start'] = datetime.utcfromtimestamp(startUXtime)
      scans[scannum]['end'] = datetime.utcfromtimestamp(endUXtime)
      scans[scannum]['source'] = lineparts[-1]
      logger.debug("parse_scan_files: scan %d from %s to %s", scannum,
                        scans[scannum]['start'], scans[scannum]['end'])
      # the raw file names in the scans file are used to get the FFT file names
      rawfiles = lineparts[3:5]
      logger.debug("parse_scan_files: looking for %s", rawfiles)
      # get the channel IDs from the fft file names
      # filename examples::
      #   16-167-001-LCP-s0001.wvsr
      #   16-167-002-LCP-s0001.wvsr
      #   16-167-001-RCP-s0001.wvsr
      #   16-167-002-RCP-s0001.wvsr
      #   or
      #   16-208-001-s0012_LCP.wvsr
      #   16-208-001-s0012_RCP.wvsr
      #   16-208-002-s0012_LCP.wvsr
      #   16-208-002-s0012_RCP.wvsr
      #   or
      #   16-238-001-s0043_d14_LCP.wvsr
      #   16-238-001-s0043_d14_RCP.wvsr
      #   16-238-002-s0043_d14_LCP.wvsr
      #   16-238-002-s0043_d14_RCP.wvsr
      for filename in rawfiles:
        rawfilename = basename(filename)
        logger.debug("parse_scan_files: using %s", rawfilename)
        rawparts = rawfilename.split('-')
        # for example: ['16', '196', '001', 's0020_RCP.wvsr']
        # or
        #              ['16', '238', '001', 's0043_d14_LCP.wvsr']
        logger.debug("parse_scan_files: parts: %s", rawparts)
        chan_str = rawparts[2]
        chanID = int(chan_str)
        logger.debug("parse_scan_files: channel ID is %s", chanID)
        if len(rawparts) == 5:
          scan_str = rawparts[-1].split('.')[0]
        else:
          lastparts = rawparts[3].split('_')
          logger.debug("parse_scan_files: last parts: %s", lastparts)
          scan_str = lastparts[0]
        scanID = int(scan_str[1:])
        logger.debug("parse_scan_files: scan ID = %s", scanID)
        # no need for station or pol; we already know that         
        namebase = '-'.join(rawparts[:2]+[chan_str]+[scan_str])
        # e.g. 16-196-001-s0020
        logger.debug("parse_scan_files: FFT file name base is %s", namebase)
        fft_filenames = glob.glob(FFTdir+namebase+"*.fft")
        logger.debug("parse_scan_files: found %s", fft_filenames)
        for index in range(len(fft_filenames)):
          namestr = basename(fft_filenames[index])
          name,ext = splitext(namestr)
          scans[scannum]["subch %d" % chanID] = namestr
  else:
    return None
  return scans

def parse_WVSR_logs(wvsrlogs):
  """
  Get metadata from WVSR FFT log files.
  
  This function gets the variable names and values from the logs. Example::
  
  
  A WVSR can be configured for two IF channels, each 160 MHz wide centered at
  80, 160, ..., 560 MHz.  Each IF channel configuration is called an 
  'experiment'. The experiment name is defined by the EXPID command.
  
  Each channel can be configured with multiple
  sub-channels of varying widths.
  
  Note
  ====
  This function has 'exec' and therefore may not have subfunctions. Also,
  the 'exec's cannot be done in a global function because they modify local
  variables
  """
  # parse the WVSR log files
  for log in wvsrlogs:
    logfile = open(log,'r')
    logtext = logfile.read()
    logfile.close()
    lines = logtext.split('\n')
    # WVSR and user
    line = lines[0].split()
    wvsrID, user = line[2], line[-3]
    EXPID = {}
    RF_TO_IF_LO = {}
    IFS = {}
    CHAN = {}
    WVSR_meta = {wvsrID: {}}
    
    for line in lines[1:]:
      if re.search("EXPID\\[[12]\\]", line):
        # parse a log file line like::
        #   16/237 08:45:13 wvsr2 EXPID[1]: EVT 301 PROGRESS: \
        #              exp 16_237_DSS-14_ARCP created by Client 3625 cjn 4-2053
        # take the variable name from the 4th item leaving off the colon
        # and add item to EXPID dict
        parts = line.split()
        logger.debug("parse_WVSR_FFT_logs: EXPID line parts: %s", parts)
        exec(parts[3][:-1]+"= "+"'"+str(parts[8])+"'")
      if re.search('RF_TO_IF_LO\[',line):
        # add to RF_TO_IF_LO for the IF channel to dict
        parts = line.split()
        exec(parts[3][:-1]+"= "+str(parts[-1]))
      if re.search("CHAN\\[.\\]: dsp", line):
        # add subchannel and its parameters to CHAN dict
        parts = line.split()
        chan_str = parts[3][:-1]
        subchan_str = "'"+parts[5]+" ' + str(int("+parts[7][:-1]+"))"
        logger.debug("parse_WVSR_FFT_logs: for CHAN doing %s with %s",
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
      if re.search("rsp_ddclo", line):
        # create subchannel offset
        parts = line.split()
        logger.debug("parse_WVSR_FFT_logs: line parts: %s", parts)
        dsp, subchan = parts[8].split(':')
        offset = int(float(parts[9])) # don't need fractional Hz
        # assign to correct channel
        for key in CHAN.keys():
          if CHAN[key]['DSP'] == dsp:
            if CHAN[key].has_key("chan_id "+subchan):
              CHAN[key]["chan_id "+subchan]["sfro"] = offset
            else:
              logger.warning(
                     "parse_WVSR_FFT_logs: no subchannel ID %s for channel %s",
                     subchan,key)
      if re.search("IFS\\[.\\]:.*PROGRESS", line):
        # IF switch inputs
        logger.debug("parse_WVSR_FFT_logs: IFS from line: %s", line)
        parts = line.split()
        logger.debug("parse_WVSR_FFT_logs: line parts: %s", parts)
        # exclude the colon from the signal source
        exec(parts[3][:-1]+"= '"+str(parts[10][:-1])+"'")
    logger.debug("parse_WVSR_FFT_logs: EXPID = %s", EXPID)
    logger.debug("parse_WVSR_FFT_logs: RF_TO_IF_LO = %s", RF_TO_IF_LO)
    logger.debug("parse_WVSR_FFT_logs: IFS = %s", IFS)
    logger.debug("parse_WVSR_FFT_logs: CHAN = %s", CHAN)
    # reformat into a single dict
    WVSR_meta[wvsrID]['user'] = user
    WVSR_meta[wvsrID]['from log'] = basename(log)
    for key in CHAN.keys():
      WVSR_meta[wvsrID][key] = {}
      WVSR_meta[wvsrID][key]['pol'] = EXPID[key][-3:]
      WVSR_meta[wvsrID][key]['expid'] = EXPID[key]
      WVSR_meta[wvsrID][key]['rf_to_if_lo'] = RF_TO_IF_LO[key]
      WVSR_meta[wvsrID][key]['IF_source'] = IFS[key]
      subchans = CHAN[key].keys()
      for subch in subchans:
        WVSR_meta[wvsrID][key][subch] = {}
        if subch[:7] == 'chan_id':
          for param in CHAN[key][subch].keys():
            WVSR_meta[wvsrID][key][subch][param] = CHAN[key][subch][param]
        else:
          WVSR_meta[wvsrID][key][subch] = CHAN[key][subch]
  return WVSR_meta

