# -*- coding: utf-8 -*-
"""
Extension to ASAP to provide access to data not read from SDFITS files by ASAP.

class supertable
================
The class supertable is an ASAP scantable with the SDFITS binary table
attached in the form of a pyFITS header and data.

The method 'get_superscan' works like the scantable 'get_scan' with scan
numbers and returns a superscan.  (A version which also retrieves scans by
source name can easily be added when needed.)

The method 'get_frequency' extracts from a supertable all scans
with a rest frequency within a specified range of a specified frequency.

'tableheader' returns header cards of possible interest to a user.

'column_names' returns the names of the table columns.

'scan_header' returns everything from a row except DATA.

'get_source_blocks' groups adjacent scans which have the same source name,
rest frequency, LSR velocity and observing mode. 'report_blocks' returns a
nicely formatted summary of all the blocks.

class finder
============
This object is associated with designated sub-directories which it will
search for SDFITS files meeting certain criteria.

The method 'find' will return file names which have a specified source
and (optionally, or) rest frequency.

module functions
================

supermerge
----------
Merge supertables into a new supertable

clobber
-------
Replace extremely by the (mean of) the adjacent value(s).

is_skinny
---------
Identifies spikes, as opposed to wider strong features

trim_extremes
-------------
Applies 'clobber' to an entire spectrum.

report_Tid_SDFITS
-----------------
Finds the blocks of identical scans, averages the scans in each
block and then plots the averaged blocks.

set_plot_layout
---------------
Arranges a given number of subplots optimally. The subplots are
individually scalled.

plot_report_scantable
---------------------
Does the reporting and plotting for each block of scans.

Things to Do
============
 # save supertables to disk
 # merge supertables
 # Frequency-switching folding
 # Make sure that scans to be averaged are compatible 
"""
from asap import average_time, \
                 scantable, \
                 plotter, \
                 merge, \
                 asaplotbase, \
                 unique
import pyfits
import glob
import os
import re
import matplotlib as mpl
import numpy as np
import logging
from os.path import isdir, basename, splitext
from Data_Reduction.FITS import fits_support
from Data_Reduction.FITS import SDFITS as S
from math import sqrt
from scipy import array, argmin, argmax, mean, sqrt, var
from matplotlib import pyplot as plt

print 'If you got this message:'
print 'QLayout: Attempting to add QLayout "" to CustomToolbarQT4Agg "", which already has a layout'
print 'run ipython without the -pylab switch'

mylogger = logging.getLogger(__name__)
mylogger.setLevel(logging.DEBUG)
mylogger.info("superasap started")

diag = False

class supertable(scantable):
  """
  A supertable is a scantable with its pyFITS HDU (header and data
  unit) included,  These can be accessed as
  supertable.header
  supertable.data
  """
  def __init__(self, antecedent, average=False, unit='K', HDU=None):
    """
    Create a scantable from a saved one or make a reference.

    Notes
    =====
    I don't know what scantable() does with the original Tid scan numbers
    but you can't count on the new ASAP scans numbers, starting from zero,
    to be contiguous. There may be missing numbers.

    Scans by ASAP Scan Number
    -------------------------
    A list of Tid scan numbers is obtained by::
      list(getscannos())

    Scans by Tid Scan Number
    ------------------------
    To get a list of Tid scan numbers::
      list(hdu[1].data.field('SCAN'))

    Selecting a Scan by Index
    -------------------------
    To select a scan by index (position in the table), you need to do::
      get_scan(list(getscannos())[index]).

    Selecting a Tid Scan by ASAP Scan Number
    ----------------------------------------
    To get information about RVSYS for ASAP 'scan' from the
    pyfits bintable, you would do::
      hdu[1].data(list(getscannos()).index(scan)).field('RVSYS').
    
    @param antecedent : string or scantable instance::
      If string, it is the name of an asap table on disk or
      the name of a rpfits/sdfits/ms file (integrations within
      scans are auto averaged and the whole file is read)
      or
      [advanced] a reference to an existing scantable

    @param average : boolean::
      average all integrations withinb a scan on read.
      The default (True) is taken from .asaprc.

    @param unit : string::
      brightness unit; must be consistent with K or Jy.
      Over-rides the default selected by the reader
      (input rpfits/sdfits/ms) or replaces the value
      in existing scantables

    @param HDU : SDFITS bintable HDI::
      Must be supplied if 'antecedent' is a scantable instead
      of a file.
    """
    if type(antecedent) == str:
      # Create a new scantable
      scantable.__init__(self, antecedent, average, unit)
      self.filename = antecedent
      self.set_freqframe('LSRK')
      # Attach the file HDU
      hdulist = pyfits.open(antecedent)
      self.header,self.data,self.cols = fits_support.get_ext_parts(hdulist[1])
    elif isinstance(antecedent,scantable):
      if HDU != None:
        # Make a new scantable from it
        scantable.__init__(self, antecedent)
        # Create an HDU to go with it
        self.header,self.data,self.cols = HDU.header, HDU.data, HDU.columns
      else:
        mylogger.warning("Making a supertable from a scantable requires an HDU")

  def new_HDU(self,scan_list):
    """
    Make a new HDU from the ASAP scans in 'scan_list'

    @param scan_list : list of ints::
      the scans numbers used to select scans from the scantable

    @return:
      an SDFITS bintable HDU
    """
    # The header and column definitions stay the same
    new_head = self.header
    new_cols = self.cols
    num_rows = len(scan_list)
    new_hdu = pyfits.new_table(new_cols, header=new_head, nrows=num_rows)
    # Make sure that the new scan list is ordered
    scan_list.sort()
    old_scan_list = list(self.getscannos())
    new_row_index = 0
    for scan in scan_list:
      old_scan_index = old_scan_list.index(scan)
      for i in range(len(new_cols)):
        new_hdu.data.field(i)[new_row_index] = \
          self.data.field(i)[old_scan_index]
      new_row_index += 1
    return new_hdu
      
  def get_superscan(self,scan_list):
    """
    Make a supertable out of the ASAP scans in the list

    @param scan_list : list of ints::
      the ASAP scan numbers used to select scans from the scantable

    @return: supertable of selected scans
    """
    # Get a new ASAP scantable:
    sct = self.get_scan(scan_list)
    # Get the associated pyFITS HDU
    new_hdu = self.new_HDU(scan_list)
    return supertable(sct,HDU=new_hdu)

  def get_frequency(self,frequency,tolerance=1):
    """
    Select all the scans within a given rest frequency range.

    @param frequency : float::
      MHz

    @param tolerance : float::
      MHz, default +/- 1 MHz
    """
    scan_list = [] # list of matching scan indices
    scannos = self.getscannos()
    for index in range(len(self)):
      if  frequency - tolerance \
        <= self.data[index].field("RESTFREQ")/1.e6 \
        <= frequency + tolerance:
        scan_list.append(scannos[index])
    num_rows = len(scan_list)
    if num_rows:
      # Get a new ASAP scantable:
      sct = self.get_scan(scan_list)
      # Get the associated pyFITS HDU
      new_hdu = self.new_HDU(scan_list)
      return supertable(sct,HDU=new_hdu)
    else:
      return None
  
  def tableheader(self):
    """
    Get the observing information from the extension header.

    Makes a dictionary of the extension header cards with
    observation information, ignoring structure cards.

    @return: dictionary
    """
    self.cards = {}
    for card in self.header.ascardlist():
      if card.key == "XTENSION":
        continue
      elif card.key == "BITPIX":
        continue
      elif card.key[:5] == "NAXIS" or card.key[:5] == "MAXIS":
        continue
      elif card.key[1:] == "COUNT":
        continue
      elif card.key == "TFIELDS":
        continue
      elif card.key[:5] == "TTYPE" or card.key[:5] == "TFORM":
        continue
      elif card.key[:5] == "TUNIT" or card.key[:4] == "TDIM":
        continue
      elif card.key == "NMATRIX":
        continue
      else:
        self.cards[card.key] = card.value
    return self.cards

  def column_names(self):
    """
    Get column names

    @return: list of strings
      Names of the columns
    """
    return self.cols.names

  def scanheader(self,row_num):
    """
    Scan header data from a row.

    Returns the data from a row, except for the spectrum.

    Notes
    =====

    This can be improved by using TTYPE to give more meaningful names
    to CRVALx and CDELTx.

    @param row_num : int::
      row number
    
    @return: dictionary
      All the row data except the spectrum.
    """
    rowdata = {}
    for name in self.column_names():
      if name != "DATA":
        rowdata[name] = self.data[row_num].field(name)
    return rowdata

  def get_source_blocks(self):
    """
    Get scan ranges for similar spectra

    This returns ranges of scans for which the source name, the
    rest frequency and the LSR velocity are the same.

    Notes
    =====
    The scan numbers designated by 'row' in the code are sequence
    numbers starting with 0.  The actual ASAP scan numbers are
    obtained from self.getscannos()[row]

    @return: list of lists::
      sequential (from 0) scan numbers
    """
    source = ""
    restfreq = 0.0
    Vlsr = -999.0
    obsmode = ""
    blocks = []
    for row in range(len(self.data)):
      if self.data[row].field('OBJECT')   != source   or \
         self.data[row].field('RESTFREQ') != restfreq or \
         self.data[row].field('VELOCITY') != Vlsr     or \
         self.data[row].field('OBSMODE')  != obsmode:
        # This row starts a new block
        block = [row]
        if blocks == []:
          # This is the first block
          blocks = [block]
        else:
          blocks.append(block)
      else:
        try:
          block.append(row)
        except UnboundLocalError,details:
          mylogger.warning("Problem appending row %d to block",row,exc_info=True)
          block = [row]
      source =   self.data[row].field('OBJECT')
      restfreq = self.data[row].field('RESTFREQ')
      Vlsr =     self.data[row].field('VELOCITY')
      obsmode =  self.data[row].field('OBSMODE')
    return blocks

  def report_blocks(self,blocks,selection):
    """
    Generates a summary table.

    Generates a table of scan blocks along with selected data about
    the scans, according to selection.  Selections are given by a
    string consisting of these options::
      S - observation scan numbers   (13 chars)
      T - date and time of observ'n
      O - object name                (17 chars)
      F - rest frequency in MHz      ( 8 chars)
      V - source velocity w.r.t. LSR ( 7 chars)
      B - beam offset (HA and dec)   (10 chars)
      R - reference position offset  (10 chars)
    Space for 68 chars, less the column width for each column, is
    available on an 80 char line
    """
    result = []
    selected = list(selection)
    if len(blocks) < 1:
      mylogger.warning("No data blocks")
    else:
      nchars = 0
      # First make the report header"
      line = "ASAP scan #s"
      for code in selected:
        if code.upper() == "S":
          line += "  Tid scan #s"
        elif code.upper() == "O":
          line += "   Object name   "
        elif code.upper() == "F":
          line += " Rest freq"
        elif code.upper() == "V":
          line += "    Vlsr"
        elif code.upper() == "B":
          line += " Beam offst"
        elif code.upper() == "R":
          line += " Ref offset"
        elif code.upper() == "T":
          line += "   Date and time "
        else:
          line += " ? "
      result.append(line)
      # Put the report for each block on a line
      scannos = list(self.getscannos())
      mylogger.debug("supertable.report_blocks: got scannos %s",str(scannos))
      for block in blocks: # 12 chars
        line = "%4d -> %4d" % (scannos[block[0]],scannos[block[-1]])
        for code in selected:
          if code.upper() == "S": # 13 chars
            first_scan = self.data[block[ 0]].field('SCAN')
            last_scan  = self.data[block[-1]].field('SCAN')
            line += " %4d -> %4d" % (first_scan,last_scan)
          elif code.upper() == "O": # 17 chars
            line += " %16s" % self.data[block[0]].field('OBJECT')
          elif code.upper() == "F": # 9 chars
            line += " %8.3f" % (self.data[block[0]].field('RESTFREQ')/1.e6)
          elif code.upper() == "V": # 8 chars
            line += " %7.1f" % self.data[block[0]].field('VELOCITY')
          elif code.upper() == "B": # 10 chars
            line += " %4.1f %4.1f" % (self.data[block[0]].field('BEAMHOFF'),
                                      self.data[block[0]].field('BEAMDOFF'))
          elif code.upper() == "R": # 10 chars
            line += " %4.1f %4.1f" % (self.data[block[0]].field('REF_HOFF'),
                                      self.data[block[0]].field('REF_DOFF'))
          elif code.upper() == "T": # 17 chars
            line += " %16s" % self.data[block[0]].field('DATE-OBS')
          else:
            line += " ?"
        result.append(line)
    return result

  def superaverage(self, weight='tintsys'):
    """
    Average together all the scans in this supertable.

    Note
    ====
    There is also a module function with this name which takes as its
    arguments a sequence of supertables

    @param weight : str::
      same as the ASAP argument for 'average_time'
    """
    mylogger.debug("Calling average_time for scan table")
    new_scantable = self.average_time(weight=weight)
    mylogger.debug("New scantable is %s", str(new_scantable))
    if new_scantable:
      new_head = self.header
      new_cols = self.cols
      new_hdu = pyfits.new_table(new_cols, header=new_head, nrows=1)
      # assume that the scan header info for the average is mostly the
      # same as for the first scan being averaged
      ASAP_scan_numbers = list(self.getscannos())
      mylogger.debug("ASAP scans: %s", str(ASAP_scan_numbers))
      for column_index in range(len(new_cols)):
        mylogger.debug("Processing column %d", column_index)
        new_hdu.data.field(column_index)[0] = \
           self.data.field(column_index)[0]
      # some data have changed in the averaging
      new_hdu.data.field('EXPOSURE')[0] = new_scantable.get_inttime(row=0)
      new_hdu.data.field('TSYS')[0]  = new_scantable.get_tsys()[0]
      new_hdu.data.field('DATA')[0] = array(new_scantable[0])
      return supertable(new_scantable, HDU=new_hdu)
    else:
      return None

  def supersave(self, filename=None, format="SDFITS", overwrite=False):
    """
    Save the supertable to disk

    The supertable has more information than a scantable.save writes.
    This method adds to the SDFITS file on disk anything that the
    scantable.save() left out.  It does this by writing a temporary
    scantable and reading it back in.  Using the saved scantable as
    a basis, it adds to its header and columns anything in the supertable
    that it does not have, and then writing the permanent file.

    For now, only format "SDFITS" is supported.
    """
    if format.upper() != "SDFITS":
      format = "SDFITS"
    tempfile = "/home/kuiper/Spec_Data/temp.fits"
    self.save(name = tempfile,
              format = "SDFITS",
              overwrite = True)
    temp_hdulist = pyfits.open(tempfile)
    temp_header,temp_data,temp_cols = \
        fits_support.get_ext_parts(temp_hdulist[1])

    # add missing header data
    temp_keys = temp_header.ascardlist().keys()
    my_keys = self.header.ascardlist().keys()
    for key in my_keys:
      if not in_list(key, temp_keys):
        if key[:5] != "TTYPE" and key[:5] != "TFORM" and key[:4] != "TDIM":
          mylogger.warning("header keyword %s is not in the saved header", key)
          temp_header.update(key, self.header[key])
    # add missing columns
    my_names = self.cols.names
    #---------
    for col in self.cols:
      my_name = col.name
      if not in_list(my_name, temp_cols.names):
        mylogger.warning("Column %s is not in the saved table", my_name)
        # But has the column been promoted to a header?
        if not in_list(my_name, temp_keys):
          #No
          temp_cols.add_col(col)
          mylogger.warning("%s is not a saved header keyword", my_name)
    # make a new table
    hdulist = S.make_SDFITS(temp_hdulist[0],temp_header.ascardlist(),temp_cols)
    if filename:
      hdulist.writeto(filename,clobber=True)
    else:
      mylogger.warning("No file name given; not saved")


class finder:
  """
  Finds the files which meet certain criteria
  """
  def __init__(self, dss=[], years=[]):
    """
    Creates a finder instance to search RA_data/SDFITS for files
    with scans meeting specific criteria in the specified stations
    and years

    @param dss : list of str::
      stations, e.g. "dss43", to be searched.  An empty list means
      all stations are included.

    @param years : list of ints::
      years to be searched.  An empty list means all years are
      included.
    """
    self.stations = dss     # list of stations to match
    self.years = {}         # dictionary of directories keyed by dss
    self.frequencies = []   # frequencies to be matched
    self.sources = []       # sources to be matched
    self.matched_files = [] # result from 'find'
    
    # If there is no list of stations then check all stations
    if self.stations == []:
      dss_list = glob.glob("/usr/local/RA_data/SDFITS/*")
      for path in dss_list:
        if isdir(path):
          self.stations.append(basename(path))
    for station in self.stations:
      self.years[station] = []
      years_list = \
        glob.glob("/usr/local/RA_data/SDFITS/"+station+"/*")
      for path in years_list:
        mylogger.debug("Checking %s for %s", path, station)
        if isdir(path):
          mylogger.warning("% is a directory",path)
          year_str = basename(path)
          year = int(year_str)
          if years == []:
            self.years[station].append(int(year_str))
          else:
            try:
              index = years.index(year)
              self.years[station].append(int(year_str))
            except ValueError:
              pass
    mylogger.debug("finder instance created for %s in %s",
                   str(self.stations), str(self.years))

  def find(self,frequencies=[], sources=[], mode="and"):
    """
    Search for the specified frequencies and sources using AND
    (default) or OR (deprecated) logic.
    
    @param frequencies : list of floats::
      frequencies to be included, to within +/- one MHz
      An empty list means all frequencies are included.

    @param sources : list of strings::
      A match means that the string provided is includes in the
      source name.  Both the test string and the source name are
      converted to uppercase. An empty list means all sources
      are included

    @param mode : string::
      The default "and" means that one item from each list must be
      satisfied for at least one scan for a file to be included.
      "or" means that any item in any list matching at least one
      scan causes the file to be included.  This mode is not
      recommended

    @return: list of strings::
      the paths from RA_data to the files found
"""
    basepath = "/usr/local/RA_data/SDFITS/"
    files_found = []
    
    # Make a list of all the keys files
    keyfiles = []
    for station in self.stations:
      for year in self.years[station]:
        # Get keyword files
        path = basepath+station+"/"+str(year)+"/"
        mylogger.debug("Checking %s", path)
        keyfiles += glob.glob(path+"*.keys")
    mylogger.debug("Frequency check files: %s", str(unique(keyfiles)))

    # Process the key files
    if keyfiles:
      # Check for frequency match
      for keyfile in keyfiles:
        freqline = os.popen("grep Frequencies "+keyfile).readline().strip()
        if freqline:
          freq_text = freqline.split(":")[1]
          freqs = freq_text.split(',')
          mylogger.debug("Frequencies:", str(freqs))
          for freq in freqs:
            for test_freq in frequencies:
              if test_freq - 1 <= float(freq) <= test_freq + 1:
                files_found.append(keyfile)
                mylogger.debug("%s matches", keyfile)
      mylogger.debug("Frequency matches found: %s", str(files_found))

      # check search mode
      if mode.lower() == 'and':
        # if AND, process the files which matched frequency
        files_to_check = files_found
        files_found = []
      else:
        # if OR, re-process all the keys files
        files_to_check = keyfiles
      mylogger.debug("Source check files: %s", str(files_to_check))

      # Check for source name match
      if files_to_check:
        for keyfile in files_to_check:
          src_line = os.popen("grep Sources " +
                      keyfile).readline().strip()
          if src_line:
            srcs = src_line.split(":")[1].split(",")
            mylogger.debug("Sources: %s", str(srcs))
            for source in sources:
              for src in srcs:
                if re.search(source.upper(),src.upper()):
                  files_found.append(keyfile)
                  mylogger.debug("%s matches", keyfile)
        result = []
        for keyfile in files_found:
          fits_file = splitext(keyfile)[0]+".fits"
          result.append(fits_file)
        mylogger.debug("Result: %s", str(result))
        result = unique(result)
        result.sort()
        self.matched_files = result
        return result
      else:
        return []
    else:
      return []
    
# ===================== generale functions of the module ====================

def in_list(item,List):
  """
  Report True if item is in list
  """
  try:
    i = List.index(item)
    return True
  except ValueError:
    return False
 
def superaverage(*args):
  """
  Average the scans in the designated tables

  The arguments may be a comma-separated supertable sequence or a list
  of supertables
  """
  supertables = get_supertables(args)
  if supertables:
    mylogger.debug("Calling asap.average_time for %s", str(supertables))
    new_scantable = average_time(supertables)
    mylogger.debug("New scantable is %s", str(new_scantable))
    if new_scantable:
      supertab = attach_hdu(supertables, new_scantable)
      mylogger.debug("New supertable is %s", str(supertab))
      return supertab
  return None

def supermerge(*args):
  """
  Merge supertables

  The argument(s) may be a comma-separated sequence of supertables or a list
  of supertables.

  @return: supertable instance
  """
  supertables = get_supertables(args)
  new_scantable = merge(supertables)
  if new_scantable:
    supertab = attach_hdu(supertables, new_scantable)
    return supertab
  else:
    return None
  
def get_supertables(*args):
  """
  Returns a supertable sequence

  The argument(s) may be a comma-separated sequence of supertables or a list
  of supertables.
  """
  if type(args[0]) == list:
    mylogger.debug("get_supertables argument is a list")
    for arg in args[0]:
      mylogger.debug("%s %s", str(type(arg)), str(arg))
    supertables = args[0]
  elif type(args[0]) == supertable:
    mylogger.debug("get_supertables argument is a supertable sequence")
    supertables = args
  elif type(args[0]) == tuple:
    mylogger.debug("get_supertables argument is a tuple")
    for arg in args[0]:
      mylogger.debug("%s %s", str(type(arg)), str(arg))
    supertables = args[0][0]
  else:
    mylogger.warning("%d improper argument type(s)", len(args))
    for arg in args[0]:
      mylogger.warning("%s %s", str(type(arg)), str(arg))
    return []
  mylogger.debug("%d supertables found for averaging", len(supertables))
  for s in supertables:
    mylogger.debug("'Supertable' is type %s", str(type(s)))
    if type(s) != supertable:
      mylogger.warning("Arguments must be supertables or a supertable list")
      return []
    else:
      mylogger.info("'Supertable' has header of length %d",
                    len(s.header.ascardlist()))
      mylogger.info("'Supertable' has %d rows", len(s.data))
  mylogger.debug("Returning: %s", str(type(supertables)))
  return supertables
  
def attach_hdu(supertables, new_scantable):
  """
  Attach an HDU to a new scantable to make a supertable.

  Notes
  =====
  There are two possible situations to handle.

  Merging
  -------
  In this case, every scan in every input supertable has a corresponding
  scan in the new scantable.  The number of scans in the new scantable
  is > 1.

  Averaging
  ---------
  All the scans in the input supertables have been collapsed into one
  scan in the new scantable.  The number of scans in the new scantable = 1.

  @param supertables : list of supertables

  @param new_scantable :
  @type  new_scantable : scantable instance

  @return: supertable instance
  """
  mylogger.debug("attach_hdu called for %s and %s",
                 str(supertables), str(new_scantable))
  # Assume that the header
  new_head = supertables[0].header
  new_cols = supertables[0].columns
  num_rows = len(new_scantable)
  mylogger.debug("The new table has %d rows", num_rows)
  new_hdu = pyfits.new_table(new_cols, header=new_head, nrows=num_rows)
  mylogger.info("%s", str(new_hdu.data))
  if num_rows == 1:
    # Take the data from the first scan in the average
    for keys in new_cols.names:
      mylogger.debug("Processing column %s", key)
      #new_hdu.data.field(key) = supertables[0].data.field(key)
    # the data have changed
    #new_hdu.data[0].field('EXPOSURE') = new_scantable.average_time()[0]
    #new_hdu.data[0].field('TSYS')  = new_scantable.get_tsys()[0]
    #new_hdu.data[0].field('DATA') = array(new_scantable[0])
  else:
    # This assumes that the new table's scans are in the same order
    # as in the ordered list of original tables
    new_row_index = 0
    for sctable in supertables:
      mylogger.info("%d scans in this supertable", len(sctable))
      for scan_index in range(len(sctable)):
        mylogger.debug("Processing scan %d", scan_index)
        mylogger.debug("%d columns in this scantable", len(new_cols))
        for i in range(len(new_cols)):
          mylogger.debug("Processing column %d row %d", i,new_row_index)
          new_hdu.data.field(i)[new_row_index] = \
            sctable.data.field(i)[scan_index]
        new_row_index += 1
  return supertable(new_scantable, HDU=new_hdu)

def clobber(data_array,index):
  """
  Replace a data array value with the adjacent value(s)

  @param data_array : numpy array

  @param index : int
  """
  if index == 0:
    data_array[index] = data_array[index+1]
  elif index == len(data_array)-1:
    data_array[index] = data_array[index-1]
  else:
    data_array[index] = (data_array[index-1] + data_array[index+1])/2.
  return data_array

def is_skinny(data_array,index):
  """
  Test whether a data value is an isolated outlier

  Returns True if the data values adjacent to the test value are
  less that 1/10 of the test value, i.e., the data point is a spike

  @param data_array : numpy array

  @param index : int

  @return: boolean
  """
  amean   = mean(data_array)
  test_value = abs(data_array[index]-amean)
  if index == 0:
    ref_value = abs(data_array[index+1] - amean)
  elif index == len(data_array)-1:
    ref_value = abs(data_array[index-1] - amean)
  else:
    ref_value = (data_array[index-1] + data_array[index+1])/2. - amean
  if test_value > 10 * ref_value:
    return True
  else:
    return False
  
def trim_extremes(data):
  """
  Remove extreme values from a data array.

  Extreme values are those greater than 10x the standard deviation
  and are 'skinny'.: numpy array

  @param data : numpy array
  """
  data_array = array(data)
  amean   = mean(data_array)
  avar    = var(data_array)
  astdev  = sqrt(avar)
  amax = data_array.max()
  amin = data_array.min()
  # Check the maximum
  if abs(amax-amean) > 10*astdev:
    index = argmax(data_array)
    #mylogger.debug("trim_extremes: testing max index")
    if is_skinny(data_array,index):
      data_array = clobber(data_array,index)
      #mylogger.info("Clobbering index %d", index)
  # check the minimum
  if abs(amin-amean) > 10*astdev:
    index = argmin(data_array)
    #mylogger.debug("trim_extremess: testing min index")
    if is_skinny(data_array,index):
      data_array = clobber(data_array,index)
      #mylogger.info("Clobbering index %d", index)
  return data_array

def set_plot_layout(nplots, rows=1, cols=1):
  """
  Find the optimal subplot layout.

  This will ensure that there are enough rows and columns for all
  the plots.  If there are not, or rows and columns are not specified,
  then the layout will be square, or cols-rows = 1.

  @param nplots : int::
    number of plots to accomodate

  @param rows : int::
    number of rows preferred

  @param cols : int::
    number of columns preferred

  @return: tuple::
    (rows,columns)
  """
  if rows*cols < nplots:
    # This estimate of the number of rows and column will be too small
    # unless nplots is an integer squared
    sidelength = int(sqrt(nplots))
    if sidelength*sidelength < nplots:
      # Let's have more columns
      sidelength += 1
    cols = sidelength
    # For this test remember that sidelength may now be one greater than
    # initially estimated
    if sidelength*(sidelength-1) >= nplots:
      sidelength -= 1
    rows = sidelength
  return rows, cols
  
def plot_report_scantable(sct,scan_ranges,restfreqs, rows=0, cols=0):
  """
  Plot spectra with independently scaled axes

  When the number of rows is not specified, it is taken to be 1.

  Note
  ====
  If we only ever have supertables then this can be made into a
  supertable method.

  @param sct : scantable or supertable instance

  @param rows : number of rows (optional)

  @param cols : number of columns (optional)

  @return: tuple::
    textx,texty,text_xstep,text_ystep
    Needed for writing the next line of text in the plot
  """
  mpl.rcParams['xtick.labelsize'] ='small'
  mpl.rcParams['ytick.labelsize'] ='small'
  mpl.rcParams['axes.formatter.limits'] = [-3,3]
  mpl.rcParams['axes.labelsize'] = 'small'
  # mpl.rcParams['font.size'] = 12.0 # 6.0
  nplots = len(sct)
  rows,cols = set_plot_layout(nplots)
  scannos = list(sct.getscannos())
  numscans = len(scannos)
  plt.clf()
  if numscans < 5:
    mpl.rcParams['font.size'] = 12.0
  elif numscans < 10:
    mpl.rcParams['font.size'] = 10.0
  elif numscans < 17:
    mpl.rcParams['font.size'] =  8.0
  else:
    mpl.rcParams['font.size'] =  6.0
  for scan in scannos:
    index = scannos.index(scan)
    x = sct.get_abcissa(scan)
    y = sct.get_spectrum(scan)
    plt.subplot(rows,cols,scan+1)
    # x is a sequence: (data,unit)
    plt.plot(x[0],trim_extremes(y))
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    text_xstep = (xmax-xmin)/10
    textx = xmin + text_xstep
    text_ystep = (ymax-ymin)/10
    texty = ymax - 2*text_ystep
    # Put the scan ranges on the plots
    SC = "%4d - %4d" % (scan_ranges[index][0],scan_ranges[index][1])
    plt.text(textx,texty,SC)
    if numscans < 37:
      # Put the source name on the plot
      source = sct.get_sourcename(scan)
      texty -= text_ystep
      plt.text(textx,texty,source)
      # Put frequencies on the plots
      RF = "%10.3f MHz" % (restfreqs[index]/1e6)
      texty -= text_ystep
      plt.text(textx,texty,RF)
    if index == (rows-1)*cols:
      plt.xlabel(r"$V_{LSR}$")
  return textx,texty,text_xstep,text_ystep
  
def report_Tid_SDFITS(filename,plot=False,report=False,outfile=''):
  """
  Reports on the scan blocks in a scantable

  This finds the blocks of identical scans, averages the scans in each
  block and then plots the averaged blocks.

  Notes
  =====
  plotter.plot(scantable) plots all the spectra over a range of
  the minimum of all the spectra and the maximum of all the spectra in
  whatever the current abscissa unit.  This only makes sense if the
  spectra are more or less overlapping.  This needs the plots to be
  independent.

  @param filename : string::
    path and name of SDFITS file

  @param plot : boolean::
    True if plots are to be generated

  @param report : boolean::
    True if summary reports are to be printed

  @return: tuple::
    super scantable, list of scan_groups, standard scantable of
    merged averages
  """
  # First make a supertable
  s_table = supertable(filename)
  scan_blocks = s_table.get_source_blocks()
  #mylogger.debug("report_Tid_SDFITS: found %d scan blocks", len(scan_blocks))
  if scan_blocks:
    # Index into scannos to get the scan number asap wants
    scannos = s_table.getscannos()
    #mylogger.debug("report_Tid_SDFITS: got scannos: %s", str(scannos))
    if scannos:
      if report:
        reports = s_table.report_blocks(scan_blocks,'SOFVBR')
        for line in reports:
          print line
      # all these lists should have the same length when done
      scan_groups = []
      scan_aves = []
      scan_ranges = []
      group_objects =[]
      group_dates = []
      restfreqs = []
      for blk in scan_blocks:
        if blk:
          # Sometimes a scantable can end with an empty record.  Filter
          # it out here:
          if float(s_table.data[blk[0]].field('RESTFREQ')) > 0.0:
            # scantable with the scans in one block
            asap_scans = []
            for scan in blk:
              try:
                asap_scans.append(scannos[scan])
              except IndexError:
                #mylogger.warning("Scan %d has no matching scan number", scan)
                break
            blk_table = s_table.get_scan(asap_scans)
            try:
              # list of average scans, one per block
              averaged_block = blk_table.average_time(align=True)
            except AttributeError, details:
              #mylogger.warning("scantable for %s cannot be averaged.",
              #                 str(blk), exc_info=True)
              #mylogger.warning("These scans are probably missing from the file")
              pass
            else:
              scan_aves.append(averaged_block)
              # list of scantables
              scan_groups.append(blk_table)
              # list of ranges of scans in each block
              scan_ranges.append([s_table.data[blk[ 0]].field('SCAN'),
                                  s_table.data[blk[-1]].field('SCAN')])
              # list of rest frequencies for each block
              restfreqs.append(s_table.data[blk[ 0]].field('RESTFREQ'))
              # source names
              group_objects.append(s_table.data[blk[ 0]].field('OBJECT'))
              group_dates.append(s_table.data[blk[ 0]].field('OBJECT'))
        else:
          #mylogger.warning("No scans in block")
          pass
      if len(scan_aves) == 0:
        return s_table,[],[]
      elif len(scan_aves) == 1:
        merged_aves = scan_aves[0]
      else:
        merged_aves = (merge(scan_aves))
      if plot:
        # merged_aves can have, at most, 16 scans
        merged_aves.set_unit("channel")
        # The merged average scantable has only one row
        spectrum = array(merged_aves.get_spectrum(0))
        if len(spectrum) == 16384:
          # DAVOS tends to have bad channels.  Remove the two extremes
          max_chan = argmax(spectrum)
          min_chan = argmin(spectrum)
          merged_aves.create_mask([min_chan,min_chan],
                                  [max_chan,max_chan], invert=True)
        merged_aves.set_doppler('RADIO')
        merged_aves.set_freqframe('LSRK')
        merged_aves.set_unit('km/s')
        text_x,text_y,text_xstep,text_ystep = \
          plot_report_scantable(merged_aves,scan_ranges,restfreqs)
        if outfile:
          #mylogger.info("Saving to %s", outfile)
          plt.savefig(outfile, dpi = 100)
      return s_table,scan_groups,merged_aves
    else:
      #mylogger.warn("No scan numbers found")
      pass
  else:
    #mylogger.warn("No scan blocks found")
    return None,[],[]