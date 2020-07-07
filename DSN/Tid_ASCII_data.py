# -*- coding: utf-8 -*-
"""
Tid_ASCII_data
==============

This module handles Tidbinbilla spectroscopy data files written in
ASCII format by the HP 9825 controller.  The data are written in
scientific notation.  The precision is such that the frequencies
are given to the nearest 10 kHz.  As a result, computed velocities
may be as much as 0.3 km/s off.

NEEDS MAJOR UPDATES FROM MODULE Observatory TO MonitorControl
"""
import logging
import os.path

import Data_Reduction.DSN.Tid_data as DRtid

logger = logging.getLogger(__name__)

def get_scan_data(data,index,block_size):
  """Get a scan from an ASCII file

  A data block is everything between two lines which begin with
  the word 'Scan'.  The scan number is also encoded as the first
  data number.  The source name, however, did not begin to
  appear in the header until later.
  
  @param data : list of strings
    All the lines from a file

  @param index : int
    a line with the word 'Scan'

  @param block_size : int
    number of lines in the block, including the line with 'Scan'
  """
  global line
  global char_index
  numbers = []
  if index+block_size > len(data):
    # Incomplete file
    return None,""
  if data[index][12:].isspace() == False:
    name = data[index][12:].strip()
    # Fix things LaTeX doesn't like
    name = name.replace("_"," ")
  else:
    name = ""
  for i in range(index+1,index+block_size):
    line = data[i].rstrip()
    for i in range(5):
      # There are five (or less) numbers of the form  3.744546E 00
      # preceded by one blank
      char_index = 1+i*13
      try:
        # There may not be five numbers on a line
        number_str = line[char_index:char_index+13]
        if len(number_str) == 13:
          new_str = number_str.replace('E ','E+')
          try:
            numbers.append(float(new_str))
          except ValueError:
            print("Could not convert to float")
            print(number_str,'->',new_str)
            print("in:",line)
            exit()
      except:
        # No more numbers on this line
        # You'd think that this would be the end of the scan but in
        # at least one case, 1989/150 scan 229, there are a bunch of
        # blank lines in the middle of the scan.  So, we keep going.
        pass
  return numbers,name

def group_files(files):
  """
  Group files by DOY

  Data for one day may be in a single file::

    95_088.tid

  or in several files named sequentially like this::

    90_342.tid, 90_342a.tid, 90_342b.tid, 90_342c.tid

  or like this::

    95_091a.tid, 95_091b.tid, 95_091b.tid

  @param files : list
    the files to be grouped

  @return: list of lists of strings
  """
  groups = []
  group = []
  for f in files:
    name = os.path.basename(f)
    logger.debug("group_files: Examining %s",name)
    if name[6] == '.':
      # This may be a standalone file or the start of a new group
      # remember its name
      first_file = name
      if group != []:
        # If the previous group was not empty, append it to groups
        groups.append(group)
        logger.debug("group_files: Completed group: %s",group)
      # Initialize a new group
      group = [name]
      if f == files[-1]:
        # If it is the last file, append it to groups because we are done
        groups.append(group)
    elif name[6] == 'a':
      # This may continue a new group or start one
      if name[:6] != first_file[:6]:
        # If the first six characters of this file are not the same as
        # those of the first file in the group, then this starts a new
        # group.  Remember the file's name
        first_file = name
        if group != []:
          groups.append(group)
          logger.debug("group_files: Completed group: %s",group)
        group = [name]
      else:
        # This continues a group
        group.append(name)
    else:
      # This handles 'b', 'c', etc.
      group.append(name)
      if f == files[-1]:
        groups.append(group)
        logger.debug("Group_files: Last file: %s",name)
        logger.debug("Completed group: %s",group)
  return groups

def num_scans_in_group(path,groups):
  """
  Number of scans in each group

  A group is a set of data files for one day.

  @param path : string
    Path to the location of the data files

  @param groups : list of lists of strings
    Each sub-list has the names of the files in the group.
  """
  scan_count = []
  for group in groups:
    # Let's find out how many scans there are in a day's observing
    num_scans = 0
    for f in group:
      fd = open(path+'/'+f,'r')
      first_line = fd.readline().split()
      # At least one file does not have the number of scans
      if first_line[1] == 'scans':
        num_scans += int(first_line[0])
      else:
        # use a default of 50
        num_scans += 50
      fd.close()
    scan_count.append(num_scans)
  return scan_count

def get_scan_indices(data):
  """Get the line numbers where scans start

  This finds all the lines starting with the word 'Scan'. Although
  some of these lines have a source name, it appears to be correct
  only for the first scans in the file.

  Notes
  =====
  This is specific to Tid ASCII files

  @param data : list of strings

  @return: list of integers
    Line numbers starting with 'Scan'
  """
  scan_indices = []
  for i in range(len(data)):
    if data[i][1:5] == 'Scan':
      # This starts a new scan.  Save a sequential index into the
      # data file for this scan.
      scan_indices.append(i)
  # Assume that all the data blocks in a file are the same size
  block_size = scan_indices[1]-scan_indices[0]
  return scan_indices, block_size

def sort_Tid_ASCII_files(files,dss):
  """
  Sort ASCII data files from Tidbinbilla, written with HP9825 controller.

  The HP9825 controller optionally wrote ASCII data files named for the
  year and day-of-year of observation.  This takes a list of filenames,
  typically generated by globbing a CD-ROM, and sorts them by time.

  @param files : Full path and name of each file to be included in the sorting.
  @type  files : list of str

  @param dss : Nominal dss number
  @type  dss : int
    
  @return: (tuple of lists of lists)
  
  For example::
  
    ([['1989-05-28T09:47:30', '1989-05-28T16:19:19', '89_148.tid']
      ['1989-05-28T16:26:07', '1989-05-30T08:27:51', '89_148a.tid']
      ['1989-05-30T08:34:13', '1989-06-03T09:01:10', '89_150.tid']
      ...
      ['1995-05-05T02:43:45', '1995-05-05T09:58:51', '95_125.tid']
      ['1995-05-06T01:01:44', '1995-05-06T09:46:46', '95_126.tid']
      ['1995-05-08T02:43:08', '1995-05-08T09:48:05', '95_128.tid']],
     [])
     The first list has the sorted file names.  The second list has
     files with bad scans
  """
  datafiles = []
  bad_scans = []
  path = os.path.dirname(files[0])
  for f in files:
    # This is a dictionary of scans indexed by scan number
    # grab all the data from a file. Mainly, this consists of rows of
    # five numbers of the form n.nnnnnnE+ee, 68 header values followed
    # by data. There is an initial line in the file saying
    # '      50 lines stored' (or some other number, up to 50), and
    # a line preceding each data block with
    # ' Scan  nn. source_name'
    fd = open(f,'r')
    data = fd.readlines()
    fd.close()

    # find the beginning of each scan.  The
    scan_indices, block_size = get_scan_indices(data)
    filedata = []
    bn = os.path.basename(f)
    for scan_index in [0,-1]:
      numbers = get_scan_data(data,
                              scan_indices[scan_index],
                              block_size)
      if numbers == None:
        # premature end of file
        bad_scans.append([bn,
                          DRtid.Tid_scan.header['SCAN'],
                          "returned no data"])
        break
      tid_scan = DRtid.Tid_scan(dss,numbers)
      yr_short = int(numbers[2])-1900
      yr_code = int(bn.split('_')[0])
      if scan_index == 0 and yr_short != yr_code:
        bad_scans.append([bn,
                          tid_scan.header['SCAN'],
                          "file name error: year="+str(int(numbers[2]))
                                                  +", name="+str(yr_code)])
      if len(numbers) == 64+tid_scan.header['MAXIS1'] or \
         len(numbers) == 65+tid_scan.header['MAXIS1']:
        # SpectraData files have 64+320 numbers
        # DAVOS files have 65+16384 numbers
        ISOtime = tid_scan.header['DATE-OBS']
        filedata.append(ISOtime)
      else:
        bad_scans.append([bn,
                          tid_scan.header['SCAN'],
                          str(tid_scan.header['MAXIS1'])
                            + "channels, "
                            + str(len(numbers))
                            + " numbers"])
        filedata.append(ISOtime)
    filedata.append(os.path.basename(f))
    datafiles.append(filedata)
  datafiles.sort()
  return datafiles, bad_scans

def get_scans(path,groups,group_index):
  """
  Build up a dictionary of Tid scans for an observing session

  This processes all the files for an observing session to
  create a Python dictionary of Tid scans, indexed by scan
  number.  Sometimes, scan numbers were re-used if scans had
  invalid data.  In case a scan already exists in a dictionary,
  the one with the latest date/time is retained.

  Notes
  =====

  Some of the scans in the ASCII data files are incomplete.
  This can be detected by the number of numbers in the scan.
  SpectraData files have 64+320 numbers. DAVOS files have
  65+16384 numbers.  Anything else is invalid and rejected.

  scantable() does not reject bad scans (e.g. all data values the same)
  but makes them invisible to get_scan(), summary() and getscannos().
  The bad scans will show up, for example, in get_sourcename() as "" and
  get_direction() as (0.0,0.0). This is also causes scans to be lumped
  together in mini-scantables. This function rejects bad scans by
  checking that the data values are not all the same.

  @param path : str::
    Path to the original data files

  @param groups : list of lists of strings::
    files grouped by DOY (observing session)

  @param group_index : int::
    index of the group to process

  @return: dictionary of Tid_scan instances::
    indexed by scan number
  """
  # This is a dictionary of row number indexed by scan number
  tid_scans = {}
  logger.debug("get_scans: doing group %d in %s", group_index, groups)
  for f in groups[group_index]:
    fd = open(path+'/'+f,'r')
    data = fd.readlines()
    fd.close()
    # find the beginning of each scan.  The
    scan_indices, block_size = Tid.get_scan_indices(data)
    logger.debug("get_scans: File %d scan indices: %s", f, scan_indices)

    # Now process all the scans in the file.
    for scan_index in scan_indices:
      # Get the numbers for the scan
      numbers,name = Tid.get_scan_data(data,
                                       scan_index,
                                       block_size)
      if numbers == None:
        print("Premature end-of-file")
        break
      tid_scan = Tid.Tid_scan(dss,numbers)
      if len(numbers) == 64+tid_scan.header['MAXIS1'] or \
         len(numbers) == 65+tid_scan.header['MAXIS1']:
        # The ASCII data have the right length

        if min(tid_scan.data) < max(tid_scan.data):
          # All data the same is nonsense
          if name and tid_scan.header['OBJECT'][0] == 'U':
            # If the scan has no name but the scan first line does
            tid_scan.header['OBJECT'] = name
          if (tid_scan.header['SCAN'] in tid_scans) == False:
            # If this scan number is not yet a key
            tid_scans[tid_scan.header['SCAN']] = tid_scan
          else:
            # Find the more recent
            if tid_scan.header['DATE-OBS'] > \
              tid_scans[tid_scan.header['SCAN']].header['DATE-OBS']:
              # replace the scan
              tid_scans[tid_scan.header['SCAN']] = tid_scan
              # Otherwise leave it be
        else:
          print("Scan",tid_scan.header['SCAN'],"has bad data")
      else:
        print("Incomplete scan",str(numbers[0])+":",len(numbers),"numbers")
        break
  return tid_scans

def report_count(filename,scan_indices):
  """
  Print the number of scans in a file

  @param filename : str::
    name of the file

  @param scan_indices:: list
    indices to lines beginning each scan

  @return: None
  """
  print(filename,'has',len(scan_indices),'scans')

def check_backend(path,fname,backends):
  """
  Determine the type of back end and the number channels.

  This updates the backends dictionary with the polarizations
  supported and the number of channels.

  @param path : str::
    path to the directory where the Tid data files are

  @param fname : str:::
    the name of the Tid data file to be processed

  @param backends : dictionary
  """
  # We need a few details to create the FITS file.  Check the first
  # file in the first group
  fd = open(path+'/'+fname,'r')
  data = fd.readlines()
  fd.close()
  scan_indices, block_size = get_scan_indices(data)
  numbers,name = get_scan_data(data,scan_indices[0],block_size)
  backends = Tid_scan.update_backend(numbers,backends)
  return backends
