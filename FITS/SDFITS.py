"""
Support for SDFITS specifically.

This module is used with Data_Reduction/FITS specifically to make
SDFITS extensions.  It also depends on methods of the Observatory
classes to configure the extension.

There is considerable variation among SDFITS-like extensions.  This
package aims for compatibility with ASAP.  ASAP requires a data cube
with at least four axes but not more than 5.  ASAP uses the
PKSIO::SDFITSreader (/asap-3.0.0/external/atnf/PKSIO/SDFITSreader.cc)
which expects FREQ, STOKES, RA, and DEC in any order. TIME or BEAM
are optional. (As per Mark Calabretta email of Wed 2010/11/10
15:31:50 [-0800].)  BEAM can be an axis type or a column name.

A typical sequence for making an SDFITS data file is::
  
  #import Observatory as O
  import Data_Reduction.FITS as F
  import Data_Reduction.FITS.SDFITS as SDF
  prihdu, hdr, cols = F.init_SDFITS(tel,num_scans)
  
  tel.FITS_init_frontends(hdr,cols,num_scans)
  tel.FITS_init_receivers(hdr,cols,num_scans)
  hdr, cols, tdim_name, tdim_str = \
      tid.FITS_init_backends( hdr,cols,num_scans,stokes_column=True)
      
  hdulist = SDF.make_SDFITS(prihdu,hdr,cols)

followed by additional header card definitions and a loop to fill all
the rows and then writing the file to disk.
"""
import astropy.io.fits as pyfits
import logging
import numpy
import Data_Reduction.FITS as F

logger = logging.getLogger(__name__)

def init_SDFITS(DSS,tablesize,time_column=False):
  """
  Initialize an SDFITS file and extension.

  Call this first to initiate an SDFITS file.  It creates the primary HDU
  suitable for a binary table extension.  It then creates a binary table with
  columns that one always needs. The data array has at least four axes for
  compatibility with ASAP.

  Notes
  =====

  The default length of axis 4 (STOKES) is 1 for compatibility with a
  "one spectrum per row" convention.  If this is not to be followed then
  subsequently called methods (probably Observatory.FITS_init_backends)
  must change MAXIS4.

  TDIMxx is calculated from MAXIS1 through MAXISn at the end.

  @param DSS : instance of the Telescope class::
    A front end with no telescope makes no sense.

  @param tablesize : int::
    Number of rows in the table

  @param time_column : boolean::
    If True, there will be a fourth data axis for time.
    Otherwise, the data axes are frequency, RA and decl.

  @return: tuple::
    The primary HDU instance, the table CardList (header)
    instance with initial header data, and the initial column set.
  """
  # create the primary HDU and extension headers
  prihdu = pyfits.PrimaryHDU()
  hdr  = pyfits.CardList()
  cols = make_basic_columns(tablesize,time_column)
  
  # add telescope location data to the table header
  logger.debug("DSS: %s", DSS)
  if type(DSS) == list:
    # This may seem odd but in the most general case there could be two or
    # more antennas, like in an interferometer.  In that case, however,
    # "single dish" FITS format doesn't apply.  We'll just assume a list of
    # length 1.
    dss = DSS[0]
  else:
    dss = DSS
  if dss !=0 :
    hdr.append(pyfits.Card('telescop', dss.name))
    hdr.append(pyfits.Card('sitelong', dss['longitude']))
    hdr.append(pyfits.Card('sitelat',  dss['latitude']))
    hdr.append(pyfits.Card('siteelev', dss['elevation']))
    hdr.append(pyfits.Card('obsgeo-x', dss['geo-x']))
    hdr.append(pyfits.Card('obsgeo-y', dss['geo-y']))
    hdr.append(pyfits.Card('obsgeo-z', dss['geo-z']))
    hdr.append(pyfits.Card('TIMESYS',  'UTC'))
    
  # there will always be four axes in the data array
  hdr.append(pyfits.Card('MAXIS',4))
  # we will always have the first data axis with frequency in the
  # from of the observatory, or time-delay for correlation functions
  # (cannot set MAXIS1 until we know the size of the spectrum)
  # hdr.append(pyfits.Card('MAXIS1',?))
  hdr.append(pyfits.Card('CTYPE1','FREQ-OBS'))
  
  # the second and third axes will be right ascension and declination
  hdr.append(pyfits.Card('MAXIS2',1))
  hdr.append(pyfits.Card('CTYPE2','RA---GLS'))
  
  hdr.append(pyfits.Card('MAXIS3',1))
  hdr.append(pyfits.Card('CTYPE3','DEC--GLS'))

  # the fourth axis is polarization.  As a default
  hdr.append(pyfits.Card('MAXIS4',1))
  hdr.append(pyfits.Card('CTYPE4','STOKES'))

  if time_column:
    # the optional fifth data axis will be time
    # (cannot set MAXIS5 until we know the number of spectra)
    # hdr.append(pyfits.Card('MAXIS4',?))
    hdr.append(pyfits.Card('CTYPE5','TIME'))
    
  return prihdu, hdr, cols

def init_bintable(obs,tel,num_scans):
  """
  Initialize an SDFITS file
  """
  # Create a blank SDFITS object
  prihdu, hdr, cols = init_SDFITS(tel,num_scans)
  prhead = prihdu.header
  logger.debug("init_bintable: Primary header:\n%s", prhead.ascardlist())
  # initialize the SDFITS file with fixed data
  tel.FITS_init_frontends(hdr,cols,num_scans)
  tel.FITS_init_receivers(hdr,cols,num_scans)
  hdr, cols, tdim_name, tdim_str = \
        obs.FITS_init_backends( hdr,cols,num_scans)
  hdulist = make_SDFITS(prihdu,hdr,cols)

  # We now have an initialized table that we can start filling
  # with generic stuff
  tbhdu = hdulist[1]
  tbhead, tbdata, tbcolumns = F.fits_support.get_ext_parts(tbhdu)
  tbhead.update(tdim_name,tdim_str)
  tbhead.update('extname', 'SINGLE DISH')
  tbhead.update('nmatrix', 1)
  tbhead.update('veldef' , 'FREQ-OBS')
  tbhead.update('projid' , 'MWAVE SPECT')
  tbhead.update('observer','UNKNOWN')
  tbhead.update('ctype1' , 'FREQ-OBS')
  tbhead.update('ctype2' , 'RA---GLS')
  tbhead.update('ctype3' , 'DEC--GLS')
  return hdulist,tbhead,tbdata,tbcolumns

def default_is_valid(format, default_value):
  """
  Checks whether the format code is consistent with the default value type

  @param format : str::
    A TFORM keyword code

  @param default_value : any::
    a value consistent with 'format'

  @return: boolean::
    True -  if the format matches the default value type
    False - if no TFORM key matched the default value type
  """
  format_types = {"I":int, "A":str, "D":float, "E":float}
  for key in format_types.keys():
    try:
      index = format.index(key)
      if format_types[key] == type(default_value):
        return True
    except ValueError:
      valid = False
  return valid
    
def make_column(name,format,default_value,tablesize,unit=''):
  """
  This creates one column in a binary table.

  @param name : str::
    The name of the column

  @param format : str::
    A TFORM keyword code

  @param default_value : any::
    a value consistent with 'format'

  @param tablesize : int::
    number of rows in the table

  @param unit : str::
    units for a physical value, e.g. deg, K, Jy

  @return: pyFITS Column instance
  """
  if default_is_valid(format, default_value):
    default_array = numpy.array([default_value]*tablesize)
    col   = pyfits.Column(name=name,format=format,array=default_array)
    return col
  else:
    return None

def make_basic_columns(tablesize,time_column=False):
  """
  Make the minimum set of columns needed by SDFITS

  This make the REQUIRED columns for an SDFITS binary table:
  * SCAN - scan number
  * CYCLE - subscan number
  * DATE-OBS - ISO format date and time
  * OBJECT - source name
  * EXPOSURE - integration time, sec
  * OBSMODE - observing mode
  * BANDWIDT - spectrometer bandwidth, Hz
  * TSYS - system temperature
  * CRVAL1 - data axis 1 reference value
  * CRPIX1 - data axis 1 reference pixel
  * CDELT1 - data axis 1 value step
  * CRVAL2 - data axis 1 reference value
  * CRPIX2 - data axis 1 reference pixel
  * CDELT2 - data axis 1 value step
  * CRVAL3 - data axis 1 reference value
  * CRPIX3 - data axis 1 reference pixel
  * CDELT3 - data axis 1 value step
  * CRVAL4 - data axis 1 reference value
  * CRPIX4 - data axis 1 reference pixel
  * CDELT4 - data axis 1 value step
  * EQUINOX - 
  * BEAMAOFF - azimuth offset
  * BEAMXOFF - cross-elevation offset
  * BEAMEOFF - elevation offset
  * BEAMHOFF - hour angle offset
  * BEAMCOFF - cross-declination offset
  * BEAMDOFF - declination offset
  * AZIMUTH -
  * ELEVATIO -

  @param tablesize : int
    Number of records (or rows) in the table.
  """
  # create empty column data arrays
  # create required columns.
  scan_col = make_column('SCAN',    '1I',  0,tablesize)
  cycl_col = make_column('CYCLE',   '1I',  0,tablesize)
  date_col = make_column('DATE-OBS','A16','',tablesize)
  obj_col  = make_column('OBJECT',  'A16','',tablesize)
  mode_col = make_column('OBSMODE', 'A8', '',tablesize)
  intg_col = make_column('EXPOSURE','1E',0.0,tablesize)
  tsys_col = make_column('TSYS',    '1E',0.0,tablesize,unit='K')
  freq_col = make_column('RESTFREQ','1D',0.0,tablesize,unit='Hz')
  obs_col  = make_column('OBSFREQ', '1D',0.0,tablesize,unit='Hz')
  cols = pyfits.ColDefs([scan_col,
                         cycl_col,
                         date_col,
                         obj_col,
                         mode_col,
                         intg_col,
                         tsys_col,
                         freq_col,
                         obs_col])
  # Add columns for the data
  
  # create columns for the first data axis: frequency
  ax1_col   = make_column('CRVAL1','1D',0.0,tablesize,unit='Hz')
  rfp1_col  = make_column('CRPIX1','1I',  0,tablesize)
  dfreq_col = make_column('CDELT1','1D',0.0,tablesize,unit='Hz')
  cols += pyfits.ColDefs([ax1_col,rfp1_col,dfreq_col])
  
  # second data axis: right ascension
  ra_col   = make_column('CRVAL2','1D',0.0,tablesize,unit='degree')
  rfp2_col = make_column('CRPIX2','1I',  0,tablesize)
  dra_col  = make_column('CDELT2','1D',0.0,tablesize,unit='degree')
  cols += pyfits.ColDefs([ra_col,rfp2_col,dra_col])
  
  # third data axis: declination
  dec_col  = make_column('CRVAL3','1D',0.0,tablesize,unit='degree')
  rfp3_col = make_column('CRPIX3','1I',  0,tablesize)
  ddec_col = make_column('CDELT3','1D',0.0,tablesize,unit='degree')
  cols += pyfits.ColDefs([dec_col,rfp3_col,ddec_col])

  # the fourth axis -- STOKES -- is made during backend initialization

  if time_column:
    # fourth data axis: Elapsed seconds since the start
    time_col  = make_column('CRVAL4','1E',0.0,tablesize,unit='sec')
    rfp4_col  = make_column('CRPIX4','1I',  0,tablesize)
    dtime_col = make_column('CDELT4','1E',0.0,tablesize,unit='sec')
    cols += pyfits.ColDefs([time_col, rfp4_col, dtime_col])
    
  # Velocity data
  rvsys_col  = make_column('RVSYS',   '1E',0.0,tablesize)
  vframe_col = make_column('VFRAME',  '1E',0.0,tablesize)
  vel_col    = make_column('VELOCITY','1E',0.0,tablesize)
  cols += pyfits.ColDefs([rvsys_col, vframe_col, vel_col])
  
  # ASAP expects a 1-relative beam number
  beam_col = make_column('BEAM',   '1I',  1,tablesize)
  eq_col   = make_column('EQUINOX','1E',0.0,tablesize)
  cols += pyfits.ColDefs([beam_col, eq_col])
  
  # beam offsets
  off_az_col = make_column('BEAMAOFF','1E',0.0,tablesize)
  off_ax_col = make_column('BEAMXOFF','1E',0.0,tablesize)
  off_el_col = make_column('BEAMEOFF','1E',0.0,tablesize)
  cols += pyfits.ColDefs([off_az_col, off_ax_col, off_el_col])
  off_ha_col = make_column('BEAMHOFF','1E',0.0,tablesize)
  off_cd_col = make_column('BEAMCOFF','1E',0.0,tablesize)
  off_de_col = make_column('BEAMDOFF','1E',0.0,tablesize)
  cols += pyfits.ColDefs([off_ha_col,off_cd_col,off_de_col])
  ref_ha_col = make_column('REF_HOFF','1E',0.0,tablesize)
  ref_de_col = make_column('REF_DOFF','1E',0.0,tablesize)
  cols += pyfits.ColDefs([ref_ha_col, ref_de_col])

  # frequency switching offset
  ref_fr_col = make_column('FOFFREF1','1E',0.0,tablesize)
  cols += pyfits.ColDefs([ref_fr_col])
  
  # azimuth and elevation
  az_col = make_column('AZIMUTH', '1E',0.0,tablesize)
  el_col = make_column('ELEVATIO','1E',0.0,tablesize)
  cols += pyfits.ColDefs([az_col, el_col])
  return cols

def make_SDFITS(prihdu,hdr,cols):
  """
  Assembles a FITS file in meomry with an SDFITS extension.

  @param prihdu : pyFITS primary HDU object

  @param hdr : a set of card images for the extension header

  @param cols : a set of (possibly empty) columns for the extension

  @return: pyFITS HDU list
  """
  header = pyfits.Header(cards=hdr)
  tbhdu = pyfits.new_table(cols, header=header)
  hdulist = pyfits.HDUList([prihdu,tbhdu])
  return hdulist


def pol_code(pol_str):
  """
  Returns numeric polarization code.

  Following the Green Bank and RPFITS conventions, the polarizations
  have been assigned numeric values.

  @param pol_str : string
    Valid polarization symbols:
      - I,Q,U,V - Stokes parameters
      - RCP or RR, LCP or LL, and cross-products RL and LR
      - HPOL or XX, VPOL or YY, and cross-products HV or XY and VH or YX

  @return: int
    in the range -6 to 4.
  """
  pstr = pol_str.upper()
  if pstr == "V":
    return 4
  elif pstr == "U":
    return 3
  elif pstr == "Q":
    return 2
  elif pstr == "I":
    return 1
  elif pstr == "RCP" or pstr == "RR":
    return -1
  elif pstr == "LCP" or pstr == "LL":
    return -2
  elif pstr == "RL":
    return -3
  elif pstr == "LR":
    return -4
  elif pstr == "HPOL" or pstr == "XX":
    return -5
  elif pstr == "VPOL" or pstr == "YY":
    return -6
  elif pstr == "XY" or pstr == "HV":
    return -7
  elif pstr == "YX" or pstr == "VH":
    return -6
  else:
    return 0

def set_offset_values(row,offset_dir,offset_amt,tbdata):
  """
  Fill beam offset columns with values
  
  This sets field values for the position offset based on what is
  in the EAC macro log.  It is a DSN extension of a GBT convention.

  @param row : int
    Table row to be filled

  @param offset_dir : string
    String as used in EAC logs: AZ, AZEL, etc.

  @param offset_amt : float or list of two floats
    Amount of offset in degrees.  If a string is given, it is
    converted to float.

  @param tbdata : SDFITS table

  @return: Boolean
    True if the offset direction was recognized; otherwise False
  """
  success = True
  if offset_dir == "AZ":
    tbdata.field('BEAMAOFF')[row] = float(offset_amt)
    tbdata.field('BEAMXOFF')[row] = 0.
    tbdata.field('BEAMEOFF')[row] = 0.
  elif offset_dir == "XEL":
    tbdata.field('BEAMAOFF')[row] = 0.
    tbdata.field('BEAMXOFF')[row] = float(offset_amt)
    tbdata.field('BEAMEOFF')[row] = 0.
  elif offset_dir == "EL":
    tbdata.field('BEAMAOFF')[row] = 0.
    tbdata.field('BEAMXOFF')[row] = 0.
    tbdata.field('BEAMEOFF')[row] = float(offset_amt)
  elif offset_dir == "AZEL":
    tbdata.field('BEAMAOFF')[row] = float(offset_amt[0])
    tbdata.field('BEAMXOFF')[row] = 0.
    tbdata.field('BEAMEOFF')[row] = float(offset_amt[1])
  elif offset_dir == "XELEL":
    tbdata.field('BEAMAOFF')[row] = 0.
    tbdata.field('BEAMXOFF')[row] = float(offset_amt[0])
    tbdata.field('BEAMEOFF')[row] = float(offset_amt[1])
  elif offset_dir == None:
    tbdata.field('BEAMAOFF')[row] = 0.
    tbdata.field('BEAMXOFF')[row] = 0.
    tbdata.field('BEAMEOFF')[row] = 0.
  elif offset_dir == "HA":
    tbdata.field('BEAMHOFF')[row] = float(offset_amt)
    tbdata.field('BEAMCOFF')[row] = 0.
    tbdata.field('BEAMDOFF')[row] = 0.
  elif offset_dir == "XDEC":
    tbdata.field('BEAMHOFF')[row] = 0.
    tbdata.field('BEAMCOFF')[row] = float(offset_amt)
    tbdata.field('BEAMDOFF')[row] = 0.
  elif offset_dir == "DEC":
    tbdata.field('BEAMHOFF')[row] = 0.
    tbdata.field('BEAMCOFF')[row] = 0.
    tbdata.field('BEAMDOFF')[row] = float(offset_amt)
  elif offset_dir == "HADEC":
    tbdata.field('BEAMHOFF')[row] = float(offset_amt[0])
    tbdata.field('BEAMCOFF')[row] = 0.
    tbdata.field('BEAMDOFF')[row] = float(offset_amt[1])
  elif offset_dir == "XDECDEC":
    tbdata.field('BEAMHOFF')[row] = 0.
    tbdata.field('BEAMCOFF')[row] = float(offset_amt[0])
    tbdata.field('BEAMDOFF')[row] = float(offset_amt[1])
  else:
    return False
