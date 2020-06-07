# -*- coding: utf-8 -*-
"""
A package for managing FITS files.

The core module creates and operates on SDFITS-format files.

SDFITS files are composed mainly of FITS binary tables.  The basic method for
creating a simple 1-table file is this::

  1. Create a primary HDU::
  
       hdu = pyfits.PrimaryHDU()
     There are optional named arguments.
      - 'data=' allows a data array to be included in the primary HDU.
      - 'header=' allows additional header data to be included, beyond
      - the bare minimum.
      
  2. Create the columns::
  
      scan_array = numpy.array([0]*tablesize)
      scan_col   = pyfits.Column(name='SCAN',format='1I',array=scan_array)
      date_array = strings.array(['']*tablesize)
      date_col   = pyfits.Column(name='DATE-OBS',format='A16',array=date_array)
      
  3. Create a column definitions object::
      cols  = pyfits.ColDefs([scan_col,date_col])
      
  4. Create a binary table::
  
      tbhdu = pyfits.new_table(cols)
     This also initializes a table header with required information.
     There are optional named arguments:
       - 'header=' allows additional header parameters to be included;
         -this is in the form of a CardList() instance
         -default None
       - 'nrows=' specifies the number of rows in the table; default 0
       - 'fill=' specifies the initial value to put in the rows; default 0
       - 'tbtype=' specifies the table type; default 'BinTableHDU'
       
  5. Create an HDU list with the primary HDU and the table::
  
       hdulist = pyfits.HDUList([hdu, tbhdu])
     Additional tables may be appended with::
       hdulist.append(tb2hdu)
"""
import astropy.io.fits as pyfits
import numpy
import logging

mylogger = logging.getLogger(__name__)

##################################### Classes #################################



##################################### Functions ################################

def save_scan_data(row, scan, dateobs, source, interval, antenna, be,
                   spectrum):
  """
  Save a scan to the FITS table
  
  This collect relevant backend data to write to FITS data table. It is for a
  "standard" (e.g. ASAP-compatible) SDFITS file initialized with the steps in
  this module's docstring.

  @param scan : Scan number
  @type  scan : int

  @param dateobs : ISO format time
  @type  dateobs : str

  @param source : Source name
  @type  source : str

  @param interval : Integration time
  @type  interval : float

  @param antenna : Antenna parameters like RA and decl.
  @type  antenna : dict

  @param be : Backend parameters
  @type  be : dict

  @param spectrum : dictionary
    Spectra indexed by polarization

  @param row : int
    Row into which to put the data
  """
  tbdata[row].setfield('SCAN',scan)
  tbdata[row].setfield('DATE-OBS',dateobs)
  tbdata[row].setfield('OBJECT',source)
  tbdata[row].setfield('EXPOSURE',interval)
  tbdata[row].setfield('CRVAL2',antenna['RA'])
  tbdata[row].setfield('CRPIX2',0)
  tbdata[row].setfield('CRVAL3',antenna['dec'])
  tbdata[row].setfield('CRPIX3',0)
  for be in backends:
    tbdata[row].setfield('BANDWIDT',be["bandwidth"])
    tbdata[row].setfield('MAXIS1',  be["num_chan"])
    tbdata[row].setfield('FREQRES', be["freqres"])
    for pol in be['polarizations']:
      tbdata[row].setfield('DATA',spectrum[pol])
      tbdata[row].setfield('STOKES',pol)
      row += 1
  return row
