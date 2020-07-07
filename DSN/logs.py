"""
Functions for extracting data from the various DSN logs.
"""
import glob
import os
import re
import shutil
import time

import Astronomy as A
import DatesTimes as DT
import Data_Reduction.DSN as DRDSN

diag = False

def eac_macro_file_name(datadir,year,DOY):
  """
  Returns the name of the EAC macro log file name.

  There should be only
  one in a observation day directory.  If there are more it returns the
  last created one.

  @param datadir : string
    Observation date subdirectory full path

  @param year : int

  @param DOY : int

  @return: string
    EAC macro log file name
  """
  pattern = "PESD-pm-%4d-%03d*.log" % (year,DOY)
  if diag:
    print("Looking for file pattern",pattern)
  files = glob.glob(datadir+pattern)
  if len(files) == 0:
    return []
  elif len(files) > 1:
    print("Warning! more than one eac macro log. Using the last one.")
    return files[-1]
  else:
    return files[0]

def find_EAC_logs(logdir,year,doy):
  """
  Adds EAC status, Tsys and five-point logs to 'log_dict'.
  
  This puts data from the pfYYYYDDD.txt, tYYYYDDD.A/B and 5pYYYYDDD.pri/sec logs
  on the EAC into a dictionary of log fil names.

  @param logdir : string
    Full path to the directory where the logs are.

  @param year : int
    Year of observation

  @param doy : int
    Day of observation

  @return:
    A dictionary with the keys:
      - track:      antenna status log names
      - Ts1:        primary channel system temperature log names
      - Ts2:        secondary channel system temperature log names
      - fivept.pri: primary channel five-point log names
      - fivept.sec: secondary channel five-point log names
  """
  global log_dict
  if logdir[-1] != '/':
    logdir += '/'
  log_dict = {}
  log_dict['track']      = find_log_file(logdir,'pf',year,'',doy,'.txt',7)
  log_dict['Ts1']        = find_log_file(logdir,'t',year,'',doy,'.A',7)
  log_dict['Ts2']        = find_log_file(logdir,'t',year,'',doy,'.B',7)
  log_dict['fivept.pri'] = find_log_file(logdir,'5p',year,'',doy,'.pri',7)
  log_dict['fivept.sec'] = find_log_file(logdir,'5p',year,'',doy,'.sec',7)
  find_minical_logs(logdir,year,doy)
  return log_dict

def find_EAC_macro_log(year,DOY):
  """
  Get the EAC macro processor log file names.
  
  This finds the macro processor log(s) which are created by the ``at`` script
  which starts the macro processor.  There maybe several log names which include
  a time.

  Notes
  =====
  This first looks for a new-style filename pattern::
  
    "PESD-pm-%4d-%03d-*.log" % (year,DOY)
    
  and then for the old style format::
  
    "PESD-pm-%4d%03d.log" % (year,DOY)
    
  The new style allows for a time code HHMM

  @param year : Year of observation.
  @type  year : int

  @param DOY : Day of observation
  @type  DOY : int

  @return: (list)Macro processor log file names
  """
  # Try the new name format first
  # new_style_macro_pattern = ("%4d-%03d-*.log" % (year,DOY))
  prefix = 'PESD-pm-'
  midfix = '-'
  suffix = '-*.log'
  # There could be more than one for this date
  files = find_log_file(DRDSN.eaclogdir,prefix,year,midfix,DOY,suffix,1)
  if files != []:
    pm_logs = []
    for f in files:
      pm_logs.append(f)
  else:
    # Try the old style format
    # old_style_macro_logname = 'PESD-pm-'+("%4d%03d.log" % (year,DOY))
    prefix = 'PESD-pm-'
    suffix = '.log'
    pm_logs = find_log_file(DRDSN.eaclogdir,prefix,year,'',DOY,suffix,1)
  return pm_logs
  
def find_log_file(logdir,prefix,year,midfix,doy,suffix,max_look_back):
  """Find a generic log file.

  Notes
  =====
  Tries the current date and then work
  backwards as much as 'max_look_back' days until one is found.
  The pseudo-code for this is::
    
    initialize output list to empty
    count = 0 # look bad only max_look_back days
    while output list is empty
      form name pattern
      if name pattern has *
        result = glob(pattern)
        if result is not empty:
          return result
        elif max_look_back exceeded:
          return empty result
        else:
          try an earlier date
      else:
        while log name == None:
          see if file exists
          if yes:
            return log name as a list
          elif max_look_back exceeded:
            return empty list
          else:
            try an earlier date
    return output list

  @param logdir : string
    Full path to the directory with the log(s)

  @param prefix : string
    Whatever precides the date, such as 'PESD-pm-'. Could include '*'

  @param year : int
    Four digit year

  @param midfix : string
    What, if anything, that goes between the year andd the day of year.
    Could be ''.

  @param doy : int
    Day of year, no leading zero(s)

  @param suffix : string
    The rest of the pattern.  Could include '*'.

  @param max_look_back : int
    Number of days before the specified DOY to look for the file

  @return: list
    File names of the latest file fitting the pattern
  """
  result = []
  count = 0
  while result == []:
    test_name = logdir+prefix+YYYYDDD_datecode(year,midfix,doy)+suffix
    if test_name != "" and re.search('\*',test_name):
     # search for possibly multiple files
      files = glob.glob(test_name)
      if files != []:
        return files
      elif count > max_look_back:
        return []
      else:
        # Try a prior day in case the log stayed open
        doy -= 1
        if doy == 0:
          year -= 1
          doy = 366
        count += 1
    else:
      # See if the file exists
      if os.path.exists(test_name):
        file_yr,file_mn,file_dy = time.gmtime(os.stat(test_name).st_mtime)[0:3]
        file_doy = DT.day_of_year(file_yr,file_mn,file_dy)
        if file_doy >= doy:
          # could be 'tomorrow' if there was a midnight transition but we
          # don't want files closed earlier than today
          return [test_name]
      elif count > max_look_back:
        return []
      else:
        # Try a prior day in case the log stayed open
        doy -= 1
        if doy == 0:
          year -= 1
          doy = 366
        count += 1
  
def find_minical_logs(logdir,year,doy):
  """Finds the EAC minical log files.

  Notes
  =====
  A log may have been started prior to the given doy so we look
  back as many as seven days.
  
  @param logdir : string
    Full path of the directory to search

  @param year : int

  @param doy : int

  @return: Null
    The results are put in global 'log_dict'
  """
  global log_dict
  pattern = logdir+"m"+YYYYDDD_datecode(year,'',doy)+'*'
  print("Initial pattern:",pattern)
  logs = glob.glob(logdir+"m"+YYYYDDD_datecode(year,'',doy)+'*')
  requested_doy = doy
  while logs == []:
    # Try the previous day
    doy -= 1
    print("Trying DOY",doy)
    pattern = logdir+"m"+YYYYDDD_datecode(year,'',doy)+'*'
    logs = glob.glob(pattern)
    if requested_doy-doy > 7:
      return []
  for f in logs:
    b = os.path.basename(f)
    if b[0] == 'm':
      # minicals
      freq = f.split('.')[1]
      key = 'cal'+str(freq)
      log_dict[key] = f
    else:
      pass
  
def get_boresights(bore_files):
  """
  Extract data from five-point boresight files.

  @param bore_files : Full paths to the boresight files.
  @type  bore_files : list

  @return: (dict) The keys are::
      - X_BW:       measured beamwidth in cross-elevation or cross-declination
      - DOY:        day of year
      - UTC:        coordinated universal time
      - SOURCE:     source name
      - HA:         hour angle
      - DEC:        declination
      - AZ:         azimuth
      - EL:         elevation
      - AX/HXPOS:   azimuth/hour angle map position
      - EL/DEPOS:   elevation/declination map position
      - FREQ:       frequency
      - POL:        polarization
      - ELOFF:      elevation offset
      - XELOFF:     cross-elevation offset
      - DECOFF:     declination offset
      - XDECOFF:    cross-declination offset
      - XPOS:       subreflector X position
      - YPOS:       subreflector Y position
      - ZPOS:       subreflector Z position
      - YOFF:       subreflector Y offset
      - ZOFF:       subreflector Z offset
      - BASE+:      baseline in the scan + direction, e.g.: bore_data['fivept.pri']['BASE+'] = ['40.7885']
      - BL+:        bore_data['fivept.pri']['BL+'] = ['40.5802']
      - BL-:        bore_data['fivept.pri']['BL-'] = ['40.6755']
      - BASE-:      baseline in the scan - direction, e.g.: bore_data['fivept.pri']['BASE-'] = ['40.6401']
      - DELTA+:     bore_data['fivept.pri']['DELTA+']   = ['0.0304588',   0.124607']
      - DELTA(0):   bore_data['fivept.pri']['DELTA(0)'] = ['0.0985753',   0.0985793']
      - DELTA-:     bore_data['fivept.pri']['DELTA-']   = ['0.0359992',   0.0969167']
      - RMS:        bore_data['fivept.pri']['RMS'] = ['0.0000','0.0000']
      - ERR:        bore_data['fivept.pri']['ERR'] = ['0.0000']
      - POS:        bore_data['fivept.pri']['POS'] = ['0.0000']
      - X_ERR:      bore_data['fivept.pri']['X_ERR'] = ['0.0000']
      - NUMPTS:     bore_data['fivept.pri']['NUMPTS'] = ['0', '0']
      - BW:         bore_data['fivept.pri']['BW'] = ['0.0000']
      - ON(1):      bore_data['fivept.pri']['ON(1)'] = ['0', '0']
      - ON(2):      bore_data['fivept.pri']['ON(2)'] = ['0', '0']
  """
  print("Bore files:",bore_files)
  data = {}
  for bore_file in bore_files:
    print("Bore file:",bore_file)
    fd = open(bore_file,'r')
    # The first line has the keys for the columns
    keys = fd.readline().strip().split()
    for key in keys:
      if key not in data:
        data[key] = []
    # The remaining lines are data
    lines = fd.readlines()
    fd.close()
    for line in lines:
      l = line.strip().split()
      if len(l) == 0:
        pass
      for item in range(len(l)):
        data[keys[item]].append(l[item])
  return data

def get_cals(cal_file):
  """
  Extract the data from a minical file.
  
  This gets the data from the minical log files and puts them in a
  dictionary keyed with 'calXXXX' where XXXX is the frequency in MHz.

  @param cal_file : File with minical data
  @type  cal_file : string

  @return: (dict of dicts) A dictionary keyed with the DOY/UT of a record.
  
  The values are also dictionaries whose keys are::
      - year:       year
      - dss:        station number
      - freq:       frequency in MHz
      - bw:         bandwidth in MHz
      - t_samp:     seconds per sample
      - n_samp:     number of samples per measurement
      - UT HH.HHH:  minical data for time HH.HHH in UT hours, consisting of
      - DOY:        day of year
      - El:         telescope elevation
      - Freq:       receiver frequency
      - Gain:       linear gain
      - Hum:        humidity
      - Lin:        linearity
      - Pol:        polarization
      - Press:      atmospheric pressure
      - R1:         reduced minical measurements
      - R2:
      - R3:
      - R4:
      - R5:
      - R6:
      - T(follow):  post-LNA follow-on temperature
      - T(lna):     LNA effective noise temperature
      - T(load):    load effective noise temperature
      - Tdiode:     effective noise diode temperature
      - Tphys:      load physical temperature
      - Tsys:       system noise temperature
      - UT:         UT in decimal hours
      - Wind_dir:   wind direction in degrees
      - Wind_vel:   wind velocity in ?
      - WxTemp:     air temperature
  """
  fd = open(cal_file,'r')
  lines = fd.readlines()
  result = {}
  data = []
  for line in lines:
    l = line.strip().split()
    if len(l) == 0:
      pass
    # following process file header data
    elif l[0] == 'Year':
      result['year'] = l[1]
    elif l[0] == 'Dss':
      result['dss'] = l[1]
    elif l[0] == 'Frequency':
      result['freq'] = l[1]
    elif l[0] == 'Bandwidth':
      result['bw'] = l[1]
    elif l[0] == 'Sample_time(sec)':
      result['t_samp'] = l[1]
    elif l[0] == 'Samples/reading':
      result['n_samp'] = l[1]
    elif l[0] == 'Reads/point':
      pass
    # This creates the keys for the row data
    elif l[0] == 'DOY':
      keys = l
    # This processes the rows
    else:
      linedata = {}
      for item in range(len(l)):
        linedata[keys[item]] = l[item]
      result["UT "+linedata['DOY']+'/'+linedata['UT']] = linedata
  fd.close()
  return result
  
def get_EAC_macro_log(year,DOY,dest_path):
  """
  Copy the EAC macro processor log
  
  This gets the macro processor log which is created by the 'at' script
  which starts the macro processor.

  Notes
  =====
  This uses find_EAC_macro_log() to get the log names.

  @param year : Year of observation

  @param DOY : Day of observation

  @param dest_path : Full path to the destination directory.

  @return: list
    EAC macro processor logs copied.
  """
  print("Entered get_EAC_macro_log for",year,DOY,dest_path)
  pm_logs = find_EAC_macro_log(year,DOY)
  if pm_logs!= None:
    # We found one or more logs
    for f in pm_logs:
      try:
        shutil.copy(f,dest_path)
        print(os.path.basename(f),"copied to",dest_path)
      except:
        print("Could not copy",os.path.basename(f),'because', sys.exc_info()[0])
  return pm_logs

def get_EAC_logs(year,DOY,dest_path):
  """
  Copy the EAC logs and scripts.to a local directory.
   
  @param year : int
    Year of observation

  @param DOY : int

  @param dest_path : string

  @return: dictionary
    EAC logs copied
  """
  eac_logs = find_EAC_logs(DRDSN.eaclogdir,year,DOY)
  for key in list(eac_logs.keys()):
    print("In get_EAC_logs, processing",key)
    # an item in the logs list could itself be a list
    if type(eac_logs[key]) == list:
      # process the item as a list
      print("EAC log files:",eac_logs[key])
      for logfile in eac_logs[key]:
        # This overwrites anything of the same name
        shutil.copy(logfile,dest_path)
        print(os.path.basename(logfile),"copied")
    else:
      # process the item as a file
      print("EAC log file:",eac_logs[key])
      shutil.copy(eac_logs[key],dest_path)
      print(os.path.basename(eac_logs[key]),"copied")
  # also get the 'process_macro' input file
  eac_logs['macro'] = find_log_file(DRDSN.eaclogdir,'PESD-pm-',year,'-',DOY,'*.log',1)
  if eac_logs['macro'] != None:
    if type(eac_logs['macro']) == list:
      for macrofile in eac_logs['macro']:
        shutil.copy(macrofile,dest_path)
        print(macrofile,"copied")
    else:
      shutil.copy(eac_logs['macro'],dest_path)
      print(eac_logs['macro'],"copied")
  return eac_logs

def get_log_data(logs):
  """
  Gets the data from EAC logs of each type.

  Extracts the minical log data, the system temperature data,
  the boresight data and the antenna tracking status data.

  @param logs : dictionary
    Dictionary of various log types.  The values are lists
    of log file names.

  @return: tuple
    Dictionaries with the actual data, 
    (cal_data,tsys_data,bore_data,track_data)
  """
  cal_data = {}
  tsys_data = {}
  bore_data = {}
  track_data = {}
  for key in list(logs.keys()):
    if key[0:3] == 'cal':
      cal_data[key] = get_cals(logs[key])
    elif key[0] == 'T':
      tsys_data[key] = get_tsys(logs[key])
    elif key[0:4] == 'five':
      print("In get_log_data:", logs[key])
      bore_data[key] = get_boresights(logs[key])
    elif key == 'track':
      track_data[key] = get_track_status(logs[key])
  return cal_data,tsys_data,bore_data,track_data

def get_on_off_times(eaclog):
  """
  Get the times when the antenna is on-source or off-source.

  This parses an EAC macro log for the times when the antenna is on
  source and when it is off.  It also gets the offset amount and direction
  and the source name. It returns a list of
  [on time, off time, source, offset, deg] list.

  @param eaclog : The full path to an EAC macro response file.
  @type  eaclog : str
    
  @return: (list of lists)  The sublist consist of::
      - on_time:   UNIX timestamp for when the antenna is on-source
      - offtime:   UNIX timestamp for when the antenna is off-source
      - source:    name of the source (string)
      - direction: direction of the offset (string)
      - amount:    the amount of the offset (degrees, float)
  """
  # Get the year and day-of-year from the filename
  parts = os.path.basename(eaclog).split('-')
  if diag:
    print(parts)
  if len(parts) == 3:
    # old style file name
    year = int(parts[2][:4])
    DOY = int(parts[2][4:])
  else:
    year = int(parts[2])
    DOY = int(parts[3])
  # get the log data
  fd = open(eaclog,'r')
  loglines = fd.readlines()
  fd.close()
  # process the data
  times = []
  on_time = None
  offtime = None
  # initially no source is defined and we are not on-source
  source_requested = False
  on_source = False
  source = None
  for line in loglines:
    parts = line.strip().split()
    if parts != []:
      if parts[1] == "source":
        source = parts[2].upper()
        source_requested = True
        on_source = False
      elif source_requested == True and re.search("Tracking",line):
        # now on source
        on_source = True
        on_time = VSR_script_time_to_timestamp(year,parts[0].replace('_','/'))
        # parse_eac_time(year,parts[0])
      elif on_source == True and parts[1] == "poffset":
        # going off source
        on_source = False
        thistime = VSR_script_time_to_timestamp(year,parts[0].replace('_','/'))
        direction = parts[2]
        if len(parts) == 4:
          amount = parts[3]
        else:
          amount = parts[3:4]
        if (on_time != None) and (thistime != on_time):
          offtime = thistime
          times.append([on_time,offtime,source,direction,amount])
      elif source == None:
        on_source = False
      elif not on_source and re.search("clr PO",line):
        on_source = True
        on_time = VSR_script_time_to_timestamp(year,parts[0].replace('_','/'))
      elif parts[1] == "source":
        source = parts[2]
        on_source = False
        thistime = VSR_script_time_to_timestamp(year,parts[0].replace('_','/'))
        if (on_time != None) and (thistime != on_time):
          offtime = thistime
          times.append([on_time, offtime, source, direction, amount])
  return times

def get_RAC_logs(year,DOY,dest_path):
  """
  Get the RAC logs

  @param year : int
    Year of observation

  @param DOY : int
    Day of observation

  @param dest_path : string
    Full path to the directory where the files should go.
  """
  rac_logs = glob.glob(DRDSN.raclogdir+'rac_'+("%4d.%03d" % (year,DOY)))
  while rac_logs == []:
    DOY -= 1
    rac_logs = glob.glob(DRDSN.raclogdir+'rac_'+("%4d.%03d" % (year,DOY)))
  for log in rac_logs:
    print("Copying",os.path.basename(log),"to",dest_path)
    # This overwrites anything of the same name
    shutil.copy(log,dest_path)
    print(os.path.basename(log),"copied")
  return rac_logs

def get_RAVI_logs(year,DOY,dest_path):
  """
  Get the RAVI logs and scripts and data

  @param year : int
    Year of observation

  @param DOY : int
    Day of observation

  @param dest_path : string
    Full path to the directory where the logs should go.

  @return: tuple
    Names of the RAVI logs, RAVI 'at' jobs, and means data files.
  """
  ravi_logs = glob.glob(DRDSN.raviLogDir+"vsr1*"+str(year)[2:]+("-%03d" % DOY)+"*.log")
  for log in ravi_logs:
    # This overwrites anything of the same name
    shutil.copy(log,dest_path)
    print(os.path.basename(log),"copied")
  ravi_at_jobs = glob.glob(DRDSN.ravi_atJobsDir
                      +"PESD-at-"+str(year)+("%03d" % DOY)+"*.ravi")
  for job in ravi_at_jobs:
    # This overwrites anything of the same name
    shutil.copy(job,dest_path)
    print(os.path.basename(job),"copied")
  means = glob.glob(DRDSN.ravi_data_dir
               +"STATS_NP1000_vsr1*"+str(year)[2:]+("-%03d" % DOY)+"*-qlook")
  for datafile in means:
    # This overwrites anything of the same name
    destination = dest_path + "/means_" + os.path.basename(datafile)[12:34]
    try:
      shutil.copy(datafile,destination)
      print(os.path.basename(datafile),"copied to",destination)
    except:
      print("Could not copy",os.path.basename(datafile),"to",destination)
  return ravi_logs, ravi_at_jobs, means
  
def get_tsys(tsys_files):
  """
  Get Tsys data.
   
  Gets system temperature data from the log files and puts them in a
  dictionary keyed with 'TsX' where X is the channel.  The keys are::
  
    - UTC:          hh:mm:ss string
    - RAC_chan:     int string
    - DOY:          int string
    - Tsys:         float string
    - Polarization: string
    - Frequency:    float string
    - Gain:         float string
    - Az:           float string
    - El:           float string
    
  Notes
  =====
  ``int string`` means a string amenable to conversion
  to a number with ``int(str)``.  Similarly for float.

  @param tsys_files : Full paths to the Tsys files.
  @type  tsys_files : list

  @return: (dict) System temperature data as lists associated with each key.
  """
  tsys_files.sort()
  data = {}
  for tsys_file in tsys_files:
    fd = open(tsys_file,'r')
    # This makes an ordered list of keys
    keys = fd.readline().strip().split()
    lines = fd.readlines()
    fd.close()
    # This initializes the data lists
    for key in keys:
      if key not in data:
        data[key] = []
    for line in lines:
      l = line.strip().split()
      if len(l) == 0:
        # ignore blank lines
        pass
      for item in range(len(l)):
        data[keys[item]].append(l[item])
  return data

def get_VSR_logs(year,DOY,dest_path):
  """
  Get logs from the VSR.
  
  This gets the VSR logs, consisting of VSR* and vrt* logs in the VSR
  log directory and PESD logs in the ops home directory.

  @param year : int
    Year of observation

  @param DOY : int
    Day of observation

  @param dest_path : string
    Full path of directory where log copies should go.

  @return: tuple
    Names of VSR logs, VSR scripts, and VSR 'at' jobs.
  """
  
  vsr_HW_logs = glob.glob(DRDSN.vsrLogDir+'VSR*log.'+("%03d-" % DOY)+'*')
  for log in vsr_HW_logs:
    # This overwrites anything of the same name
    shutil.copy(log,dest_path)
    print(os.path.basename(log),"copied")
  ses_date = dest_path.split('/')[-2]
  print("Session",ses_date)
  vsr_SW_logs = glob.glob(DRDSN.obs_dir+ses_date+"/vrt*log."+("%03d" % DOY)+"*")
  for log in vsr_SW_logs:
    # This overwrites anything of the same name
    shutil.copy(log,dest_path)
    print(os.path.basename(log),"copied")
  
  vsr_scripts = glob.glob(DRDSN.vsrScriptDir+'PESD-'+str(year)+'-'+("%03d" % DOY)+'*.vsr')
  for script in vsr_scripts:
    # This overwrites anything of the same name
    shutil.copy(script,dest_path)
    print(os.path.basename(script),"copied")
    
  vsr_at_jobs = glob.glob(DRDSN.vsr_atJobsDir+"PESD-at-"+str(year)+'-'+("%03d" % DOY)+"*.vsr")
  for job in vsr_at_jobs:
    # This overwrites anything of the same name
    shutil.copy(job,dest_path)
    print(os.path.basename(job),"copied")
  vsr_logs = vsr_HW_logs + vsr_SW_logs
  return vsr_logs, vsr_scripts, vsr_at_jobs

# These functions are unfinished
    
def report_logs(cal_data,tsys_data,bore_data,track_data):
  """
  Prints out the logs of each type.
  
  Notes
  =====
  This function needs more work to be useful.

  @param cal_data : dictionary
    Data from minicals

  @param tsys_data : dictionary
    Data from system temperature logs

  @param bore_data : dictionary
    Data from boresights

  @param track_data : dictionary
    Data from antenna status (pf*) files.

  @return: None
  """
  print("cal_data:")
  for key in list(cal_data.keys()):
    print(key)
  print("\ntsys_data:")
  for key in list(tsys_data.keys()):
    print(key,"for",tsys_data[key]['Frequency'][0])
  print("\nbore_data:")
  for key in list(bore_data.keys()):
    print(key,"for",bore_data[key]['FREQ'][0])
  print("\n track status data")
  for key in list(track_data.keys()):
    print(key)

def parse_boresight_line(line):
  """Parse a line of boresight data.
   
  Proxesses data from an EAC ``5pYYYYDDD.xxx`` file, where
  xxx is ``pri`` (primary, the preferred frequency for the boresight) or ``sec``
  (secondary). ``line`` is a line from the log file.  The columns in the table
  are::
  
    - DOY:      day of year (3 digit integer)
    - UTC:      Coordinated Universal Time (HH:MM:SS)
    - SOURCE:   name of the source used for boresight
    - FREQ:     receiver frequency used (MHz)
    - POL:      feed polarization (LCP/RCP)
    - HA:       approximate hour angle (deg) at center position
    - DEC:      approximate declination (deg) at center position
    - AZ:       approximate azimuth (deg) at center position
    - EL:       approximate elevation (deg) at center position
    - AX_HXPOS: total initial cross-elevation or cross-declination offset (deg)
    - EL_DEPOS: total initial elevation or declination offset (deg)
    - DELTAP_X: delta Ta (K) at X-direction plus offset (nominallY +HPBW/2)
    - DELTAZ_X: delta Ta (K) at X-direction zero offset
    - DELTAM_X: delta Ta (K) at X-direction minus offset (nominally -HPBW/2)
    - NUMPTS_X: always 0 (number of samples at a position?)
    - RMS_X:    always 0.0 (r.m.s. of sample mean)
    - BASEP_X:  Ta (K) of baseline in plus X offset direction
    - ON1_X:    always 0
    - ON2_X:    always 0
    - BASEM_X:  Ta (K) of baseline in minus X offset direction
    - BW_X:     measured beamwidth (deg) in X direction
    - ERR_X:    estimated error (deg) of beamwidth in X direction
    - POS_X:    total boresight offset (deg) in X not including manual offset
    - DELTAP_Y: delta Ta (K) at Y-direction minus offset (nominallY -HPBW/2)
    - DELTAZ_Y: delta Ta (K) at Y-direction zero offset position
    - DELTAM_Y: delta Ta (K) at Y-direction plus offset (nominallY +HPBW/2)
    - NUMPTS_Y: always 0
    - RMS_Y:    always 0.0
    - BASEP_Y:  Ta (K) of baseline in plus Y offset direction
    - ON1_Y:    always 0
    - ON2_Y:    always 0
    - BASEM_Y:  Ta (K) of baseline in minus Y offset direction
    - BW_Y:     measured beamwidth (DEG) in Y direction
    - ERR_Y:    estimated error (deg) of beamwidth in Y direction
    - POS_Y:    total boresight offset (deg) in Y not including manual offset
    
  The following summarize the boresight results.  One pair of columns should
  be the same as POS_X and POS_Y and the other pair calculated from these,
  depending on whether the boresight was done in XDEC,DEC or XEL,EL::
  
    - XELOFF:   total boresight offet (deg) in cross-el, not incl. manual offs.
    - ELOFF:    total boresight offet (deg) in elevation, not incl. manual offs.
    - XDECOFF:  total boresight offet (deg) in cross-dec, not incl. manual offs.
    - DECOFF:   total boresight offet (deg) in declination, not incl. manual offs.
    - The following may vary according to a subreflector position model
    - YPOS:     subreflector Y axis position (inches)
    - ZPOS:     subreflector Z axis position (inches)
    
  The following are offsets from the model::
  
    - YOFF:     subreflector Y axis offset (inches)
    - ZOFF:     subreflector Z axis offset (inches)
    
  The information extracted from this and returned in a dictionary is::
  
    - time:     UNIX timestamp in seconds since the epoch 1970.0
    - SOURCE:
    - FREQ:
    - POL:
    - HA:
    - DEC
    - AZ:
    - EL:

  Notes
  =====
  This function is unfinished.
  """
  T = logtime_to_timetuple(line[1])
  sec_since_O_UT = ((T[0]*60)+T[1]*60)+T[2]
  sec_since_Jan0 = (line[0]-1)*24*60*60
  return sec_since_Jan0

def get_track_status(status_files):
  """I've forgotten what is is supposed to do."""
  for status_file in status_files:
    fd = open(status_file,'r')
    fd.close()
    return None

def parse_tsys_log_line(line):
  """
  Parse a line in a Tsys log.
   
  A Tsys log has the name tYYYYDDD.x (where x is A or B).
  It has the following columns::
  
    - DOY:          (3-digit int)
    - UTC:          (HH:MM:SS) Coordinated Universal Time
    - Tsys:         (K) - system temperature
    - RAC_chan:     integer
    - Frequency:    (MHz)
    - Polarization: (string) LCP or RCP
    - Gain:         (K/W)
    - Az:           (deg)
    - El:           (deg)

  Notes
  =====
  This function is unfinished
  """
  T = logtime_to_timetuple(line[1])
  sec_since_O_UT = ((T[0]*60)+T[1]*60)+T[2]
  sec_since_Jan0 = (line[0]-1)*24*60*60
  return sec_since_Jan0
