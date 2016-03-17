# -*- coding: utf-8 -*-
"""
Routines for managing spectroscopic FITS files
"""

import glob
import os
import os.path as P
import pyfits
import Astronomy as A
import logging

module_logger = logging.getLogger(__name__)

debug = 0
default_scan = 0

def list_FITS_files(dataDir):
  """
  Build up a list of available scans in the selected directory.
  
  It creates a list of standard FITS data files (one spectrum per file
  in image format), each element itself a list containing:
    - scan       - scan number taken from FITS header
    - date       - date of observation
    - time       - time of observation
    - source     - source name
    - restfreq   - spectral line rest frequency (MHz)
    - vlsr       - velocity of source w.r.t. Local Standard of Rest (km/s)
    - integ      - integration time
    - instrument - instrument usedd to obtain the data
    - f          - file name

  Notes
  =====
  This is for old FITS files written by exp_control, etc.
  
  @param dataDir : string -
    Where the datafiles are
  """
  module_logger.debug("\nChecking %s",dataDir)
  filenames = glob.glob(dataDir+'/*.fits')
  module_logger.debug("\nFiles: %s",filenames)
  scanlist = []
  default_scan = 0
  path,dummy = P.split(filenames[0])
  for f in filenames:
    source, restfreq, vlsr, scan, integ, date, time, instrument \
                = get_FITS_header(f)
    scanlist.append([scan,date,time,source,restfreq,vlsr,integ,instrument,f])
  scanlist.sort()
  return scanlist

def get_FITS_header(f):
  """
  Gets summary information from a standard FITS image file.

  @param f : string
    Name of an image FITS file
  """
  hdulist = pyfits.open(f)
  h = hdulist[0].header
  try:
    if h['extend'] == True:
      # This is an extended header; try the next block
      h = hdulist[1].header
      print "Extended header:\n",h
    else:
      # This is a simple FITS file
      pass
  except:
    # This must be a simple FITS file
    pass
  try:
    source = "%12s" % (h['object'].strip())
  except KeyError:
    source = "Unknown source"
    print "Header:\n",h
  try:
    restfreq = "%7.3f" % (float(h['restfreq'])/1e9)
  except KeyError:
    print "Keyword RESTFREQ not found."
  try:
    vlsr = "%5.1f" % (float(h['velo-lsr'])/1e3)
  except KeyError:
    print "Keyword VELO-LSR not found"
    vlsr = "%5.1f" % 0.0
  try:
    scan = int(h['scan-num'])
  except KeyError:
    print "Keyword SCAN-NUM not found"
    default_scan += 1
    scan = "%4d" % default_scan
  except ValueError:
    scan = h['scan-num']
  # integration time
  try:
    integ = "%7.1f" % float(h['obstime'])
  except KeyError:
    print "Keyword OBSTIME not found"
    integ = "%7.1f" % 0.0
  # time of observation
  try:
    [date_obs, fraction] = h['date-obs'].split('.')
  except ValueError:
    date_obs = h['date-obs']
  try:
    module_logger.debug(" Date of observation:",date_obs)
    T = A.isotime_to_datetime(date_obs)
    module_logger.debug("%s ", T)
    date = T.date()
    time = T.time()
  except :
    print "Could not resolve DATE-OBS"
  # backend
  try:
    # The following is probably site-specific.
    if h['origin'] == 'CLASS-Grenoble':
      # This applies at CSO
      instrument = h['telescop'][4:]
    else:
      instrument = h['instrume']
  except KeyError:
    # Note a spectral line file
    print "Don't know data source in", f
  return source, restfreq, vlsr, scan, integ, date, time, instrument

def find_FITS_dirs():
  """
  Report all the directories with files which have the extent .fits
  """
  fd = os.popen("locate *.fits")
  response = fd.readlines()
  fd.close()
  result = []
  for item in response:
    directory = P.dirname(item.strip())
    if result.count(directory) == 0:
      result.append(directory)
  result.sort()
  return result

def read_sdfits(file):
  """
  Report on an SDFITS file

  Prints the FITS file structure, the primary header, and
  the first table header, as well as the first row.

  Notes
  =====
  This differs slightly from get_ext_parts(tbhdu) in that it
  handles a table HDU already in memory.

  @param file : string -
    Name of an SDFITS file
  """
  print "\nProcessing",file
  hdus=pyfits.open(file)
  print "File contents:\n",hdus.info()

  print "Primary header:\n",hdus[0].header

  tbhead = hdus[1].header
  print "Table header:\n",tbhead
  tbdata = hdus[1].data
  ncols = tbhead['tfields']
  print ncols,"columns"
  lendata = len(tbdata)
  print lendata,"rows"

  colheads = []
  for i in range(1,ncols+1):
    key = "ttype"+str(i)
    colheads.append(tbhead[key])

  hdus.close()
  return tbhead,tbdata,colheads

def select_FITS_file(pattern):
  """
  Select a FITS file.

  This examines the files in the current directory tp see if any
  start with the string SIMPLE which marks the start of a FITS header.
  All the fits files are then listed to allow one to be selected.

  @param pattern : string -
    A file name pattern suitable for globbing.
  """
  files = glob.glob(pattern)
  for i in range(len(files)):
    if P.isfile(files[i]):
      fd = open(files[i])
      card1 = fd.read(80)
      fd.close()
      if card1[:6] == "SIMPLE":
        print i,")",files[i],"is FITS"
  f = int(raw_input("Select file by number: "))
  return files[f]

def get_ext_parts(tbhdu):
  """
  Report on a FITS HDU
  
  This reports on the FITS extension structure, the header, and
  the first and last rows.

  @param tbhdu : extension HDU object
  """

  tbhead  = tbhdu.header
  tbdata  = tbhdu.data
  tbcols  = tbhdu.columns
  
  nrows = len(tbdata)
  ncols = len(tbcols)
  module_logger.debug("Table header: %s\n",tbhead.ascardlist())
  module_logger.debug(" %d columns", ncols)
  module_logger.debug(" %d rows", nrows)

  colheads = []
  for i in range(1,ncols+1):
    key = "ttype"+str(i)
    colheads.append(tbhead[key])
    module_logger.debug("%s",format_sdfits_row(tbdata,0,colheads))
    last_row = nrows-1
    module_logger.debug("%s",format_sdfits_row(tbdata,last_row,colheads))
  return tbhead,tbdata,tbcols

def format_sdfits_row(data,row,colheads):
  """
  Prints the contents of a row in an SDFITS file.

  @param data : pyFITS HDU data object

  @param row : pyFITS row object

  @param colheads : pyFITS BINTABLE column names

  @return: str
  """
  rowdata = data[row]
  text = "\nRow %d\n" % row
  for i in range(len(colheads)):
    text += "%d: %s = %s" %  (i+1,colheads[i],rowdata[i])

