# -*- coding: utf-8 -*-
"""
Functions for scripts, logs and data files on the EAC, RAC, RAVI and VSR.

A lot of stuff here is either obsolete or no longer accurate

This package has functions to generate scripts, download and parse DSN logs,
and manage data files. The module also initializes path names relative
to a root directory whose default is::

  root_dir = "/usr/local/projects/PESD/"

The observations are kept in subdirectories with names like YYYY-MM-DD in::

  obs_dir  = root_dir + "observations/"

Notes
=====
This module and the others in this package assume that sshfs has been used to
mount the relevant Flight Ops radio astronomy hosts as sub-directories of
localhost:/usr/local/projects/PESD/rahost.  The logs and some post-processed
data files will be or have been copied into
localhost:/usr/local/projects/PESD/observations/YYYY-MM-DD.  This directory
will be created when 'make_scripts' is run.
"""

diag_read = False
diag = False
STATS_binary_record_size = 76

root_dir = "/usr/local/projects/PESD/"
obs_dir  = root_dir + "observations/"

import logging
import os
import re
import shutil
import socket
import struct
import subprocess
import sys
import time

import DatesTimes as DT

mylogger = logging.getLogger("__main__."+__name__)

def send_ssh_command(local_cmd,remote_cmd):
  """
  Sends an ssh comand to a remote host.
  
  @param local_cmd : string
    Something like 'ssh -l napoleon venus-eac1'
  @param remote_cmd : string
    What is to be executed on the remote host, like 'date'
  """
  # proc_in,proc_out,proc_err = os.popen3(local_cmd+' '+remote_cmd)
  p = subprocess.Popen(local_cmd+' '+remote_cmd,
            shell=True,
          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (proc_in,
   proc_out,
   proc_err) = (p.stdin, p.stdout, p.stderr)
  STDOUT = proc_out.readlines()
  proc_out.close()
  STDERR = proc_err.readlines()
  proc_err.close()
  proc_in.close()
  return STDOUT,STDERR

def verify_ssh_mount(server):
  """
  Verify that a server is mounted.

  @param server : string
    One of 'EAC', 'VSR', 'RAVI', 'RAC'.  These names are 
    used for the scripts in /usr/local/projects/PESD/bin.
  """
  servername = {}
  servername['EAC'] = 'venus-eac'
  servername['VSR'] = 'venus-vsr'
  servername['RAVI'] = 'ravi2'
  servername['RAC'] = 'venus-rac'
  svr = server.upper()
  # The following tests for something unique in the server file system
  if svr == 'EAC':
    test = "scripts"
  elif svr == 'RAC':
    test = "nvidia"
  elif svr == 'VSR':
    test = 'var/home/ftp/vsr/scripts'
  elif svr == 'RAVI':
    test = 'media'
  target = root_dir + ""+svr+"/"+test
  if diag:
    print(("Looking for",target))
  if not os.path.exists(target):
    try:
      os.system(root_dir + "bin/"+servername[svr]+"-mount")
      print((svr,"has been mounted"))
      return True
    except:
      print(("sshfs mount",svr+":/home/ops with"))
      print(("   /usr/local/projects/PESD/bin/"+servername[svr]+"-tunnel"))
      print(("   /usr/local/projects/PESD/bin/"+servername[svr]+"-mount"))
      return False
  else:
    return True

def check_vsr():
  """
  Check the usage of the VSR BLS disk area.

  Reports names of files and amount of disk used.
  """
  fd = os.popen('/usr/local/projects/PESD/bin/venus-vsr-ssh ls -l /bls/meta')
  response = fd.readlines()
  fd.close()
  timecode = {}
  for line in response:
    items = line.strip().split()
    if len(items) == 9:
      timedata = ' '.join(items[5:8])
      name = items[8]
      if re.search(":",items[7]):
        timestr = time.ctime().split()[-1]+' '+ timedata
        timecode[name] = time.strptime(timestr,"%Y %b %d %H:%M")
      else:
        timecode[name] = time.strptime(timedata,"%b %d %Y")
  fd = os.popen('/usr/local/projects/PESD/bin/venus-vsr-ssh bls_ls')
  response = fd.readlines()
  print(response)
  fd.close()
  size = {}
  for line in response:
    print((">",line))
    items = line.strip().split()
    if len(items) < 4:
      name = items[0]
      size[name] = int(items[1])/1024./1024/1024
    else:
      total = int(items[0])/1024./1024
      free = int(items[3])/1024./1024
  return total,free,timecode,size

def clear_vsr(year,DOY):
  """
  Delete the data files on venus-vsr1 for a given date.

  @param year : int
    Year of observation
  @param DOY : int
    Day of observation
  """
  ssh_command = root_dir + "bin/venus-vsr-ssh "
  std_in, std_out, std_err = \
    os.popen3(ssh_command+"bls_ls")
  status = std_err.readline().strip()
  std_in.close()
  std_err.close()
  if status == '':
    response = std_out.readlines()
    std_out.close()
    sessionID = "%2s-%03d-" % (str(year)[2:],DOY)
    for line in response:
      if re.search(sessionID, line.strip()):
        data_file = line.split()[0]
        command = "bls_rm "+data_file
        print(command)
        std_in, std_out, std_err = \
          os.popen3(ssh_command+command)
        status = std_err.readline().strip()
        std_in.close()
        std_out.close()
        std_err.close()
        print(status)
  else:
    print(("Error:",status))


def backup_one_line(fd):
  """
  Moves the file pointer to right after the previous newline.

  @param fd : file descriptor
  Notes
  =====
  It ignores the current character in case it is a newline.
  """
  index = -2
  c = ""
  while c != "\n":
    fd.seek(index,2)
    c = fd.read(1)
    index -= 1

def get_last_line(fd):
  """
  Gets the last line in an ASCII file.
  
  This is useful for getting the last line in a very large ASCII file.
  It also reports of the last line was completed with a newline or interrupted.

  @param fd : file descriptor
  """
  fd.seek(-1,2)
  clast = fd.read(1)
  if clast == "\n":
    complete_line = True
  else:
    complete_line = False
  backup_one_line(fd)
  line = fd.readline()
  return (complete_line, line)

def header_report(year,doy,start_sec,freq,spc,vsr,nchan,bw,bps,nsamps):
  """
  Format a report of the observing parameters

  This reports parameters such as those found in a RAVI or VSR data file,
  such as those in a STATS file header.

  @param year : Year the file was made
  @type  year : int
    
  @param doy : Day of year that the file was made
  @type  doy : int
    
  @param start_sec : Seconds since midnight for the record
  @type  start_sec : int
    
  @param freq : Frequency in MHz
  @type  freq : float
    
  @param spc : Number identifying the SPC where the VSR resides
  @type  spc : int
    
  @param vsr: Number identifying the VSR
  @type  vsr: int
    
  @param nchan: Number of channels recorded
  @type  nchan: int
    
  @param bw : Recording bandwidth
  @type  bw : float
    
  @param bps : Bits per sample
  @type  bps : int
    
  @param nsamps : Samp/s (not the inverse of the bandwidth if data were averaged.
  @type  nsamps : int
  """
  report = []
  report.append('send_host   = '+os.popen('hostname').readline().strip())
  report.append('year        = '+str(year))
  report.append('doy         = '+str(doy))
  report.append('UTstart     = '+str(start_sec)+' // seconds')
  report.append('Freq        = '+str(freq)+' // MHz')
  report.append('VSR         = '+str(vsr))
  report.append('# channels  = '+str(nchan))
  report.append('bandwidth   = '+str(bw)+' // MHz')
  report.append('bits/sample = '+str(bps))
  report.append('stats/sec   = '+str(nsamps))
  return report

def get_obs_date(filename):
  """
  Parse a filename for the date of observation.

  Given a VSR or RAVI data file, return the year, month, day of
  observation as a tuple.  Generally, the data file has the year and DOY
  in its header but, if not, they are encoded in the file name.

  Notes
  =====
  
  A raw VSR UNIX file has the form::
  
    vsr1X.CCC.YY-DDD-HHMM.raw
    
  where X = a|b, CCC = 1w1 | 2w1 (and possibly 1N1 or 2N1) and the rest are
  the date info. The 1 in front of the X could be 2 or 3 but not at DSS-13.

  A signal statistics file has the form::
  
    STATS_NP1000_vsr1X.CCC.YY-DDD-HHMM
    
  encoded as above. It may also have the suffix -qlook and -bin, and HHMM
  may be 'mars' for some older files.

  FFT file names are::
  
    FFT_vsr1X.CCC.YY-DDD-HHMM.raw_DDD-HH:MM:SS_HH:MM:SS
    
  where the part after 'raw' gives the start and end of the data segment
  processed. The HHMM may be 'mars' for some older files, in which case the
  file name form is FFT_vsr1X.CCC.YY-DDD-mars_097-04:19:45_04:20:00
  
  @param filename : string
  """
  # This gets rid of "STATS_NP1000_", if any
  filename_parts = filename.split('_')
  if filename_parts[0] == 'STATS':
    filename = filename_parts[2]
  elif filename_parts[0] == 'FFT':
    filename = filename_parts[1]
  parts = filename.split('.')
  items = parts[2].split('-')
  yr = 2000 + int(items[0])
  doy = int(items[1])
  return calendar_date(yr,doy)

def set_file_pointer(fd,new_pos):
  """
  Position the file pointer in a binary file.

  Notes
  =====
  This overcomes a bug in Python 2.5.2 which does not allow Python
  long int offsets.

  @param fd : file descriptor

  @param new_pos : long
  """
  pos = fd.tell()
  if diag:
    print(("Current file position =",pos))
  if new_pos != pos:
    whence = 0
    if new_pos > 2000000000: # cover for a bug in Python 2.5.2
      fd.seek(2000000000,whence)
      whence = 1
      new_pos = int(new_pos - 2000000000)
      while new_pos > 2000000000:
        fd.seek(2000000000,whence)
        new_pos = int(new_pos - 2000000000)
    fd.seek(new_pos,whence)
    if diag:
      print(("New file position =",fd.tell()))
  
def set_binary_record_pointer(fd,file_header_size,record_size,index):
  """
  Position the file pointer to the record at the specified record index.

  @param fd : file descriptor

  @param file_header_size : int
    The header size may be zero.

  @param record_size : int

  @param index : int
    Record index starting from 0.  The Python -1 for the last record is
    not allowed.
  """
  pos = fd.tell()
  current_index = (pos - file_header_size)/record_size
  if diag:
    print(("Current record index =",current_index))
    print(("Moving to",index))
  if index != current_index:
    new_pos = pos + record_size*(index-current_index)
    set_file_pointer(fd,new_pos)

def get_binary_record(fd,file_header_size,record_size,index):
  """
  Get the binary record at the specified index.

  @param fd : file descriptor

  @param file_header_size : int
    The file may have a file header, or 'file_header_size' may be 0.

  @param record_size : int
    This includes a record header if there is one.

  @param index : int
    Starting from zero.
  """
  set_binary_record_pointer(fd,file_header_size,record_size,index)
  buf = fd.read(record_size)
  return buf

def print_header(header):
  """
  Prints selected contents from a VSR data record header.

  @param header : tuple
    (spcid, vsrid, chanid, bps, srate, errflg, year, doy, sec, freq, orate, 
    nsubchan)
  """
  spcid, vsrid, chanid, bps, srate, errflg, year, doy, sec, freq, orate, \
    nsubchan = header
  print(("SPC %2d, VSR %1d, channel %1d" % (spcid,vsrid,chanid)))
  print(("%5f ksamp/sec, %2d b/samp" % (srate/1000,bps)))
  print(("error flag =",errflg))
  TS = DT.VSR_tuple_to_timestamp(year,doy,sec)
  print(("date/time at start: ",DT.timestamp_to_str_with_ms(TS)))
  print(("Frequency = %9.3f MHz" % (freq/1e6)))
  print(("Averages/sec:",orate))
  print((nsubchan,"subchannels"))

def rx_band(freq):
  """
  Associate a waveguide band with a frequency
  
  """
  # Which front end channel is this for?
  if freq < 2000.:
    band = 'L'
  elif freq < 5000.:
    band = 'S'
  elif freq < 10000.:
    band = 'X'
  elif freq < 40000.:
    band = 'Ka'
  else:
    band = None
  return band

def path_to_remote(workstation, local_path):
  """
  """
  thishost = socket.gethostname()
  if thishost == workstation:
    # program runs on this host
    datapath = local_path
  elif thishost == 'gpu1' or thishost == 'gpu2':
    return local_path
  else:
    workstationpath = "/home/kuiper/mnt/"+workstation
    if thishost == 'dto':
      datapath = workstationpath+local_path
    else:
      dtopath = "/home/kuiper/mnt/dto"
      if thishost == 'kuiper':
        datapath = dtopath +workstationpath+local_path
      else:
        return None
  return datapath

