"""
Opens SDFITS file for analysis

The command line arguments --project, --date and --DSS will open the files in
the correct session directory but does not check if the file is actually from
that session.  The corresponding values in the file are used.

For each SDFITS file found it creates an item item in 'header' and an item in
'examiner' keyed with the number of the order in which they were found, which
is usually alphanumerically in file name.

Examples::
  In [8]: header[0]
  Out[8]: 
  SIMPLE  =                    T / conforms to FITS standard                      
  BITPIX  =                    8 / array data type                                
  NAXIS   =                    0 / number of array dimensions                     
  EXTEND  =                    T                                                  
  BLOCKED = 'T       '                                                            
  DATE    = '2017/08/02'                                                          
  ORIGIN  = 'FITSfile.__init__'                                                   
  TELESCOP= 'DSS-43  '                                                            
  SITELONG=        211.019942862 / degrees west of Greenwich                      
  SITELAT =        -35.403983527 / degrees                                        
  SITEELEV=              688.867 / meters                                         
  OBSGEO-X=         -4460894.585 / meters                                         
  OBSGEO-Y=          2682361.554 / meters                                         
  OBSGEO-Z=          -3674748.58 / meters 
  In [9]: examiner[0].header
Out[9]: 
XTENSION= 'BINTABLE'           / binary table extension                         
BITPIX  =                    8 / array data type                                
NAXIS   =                    2 / number of array dimensions                     
NAXIS1  =              3670834 / length of dimension 1                          
NAXIS2  =                    2 / length of dimension 2                          
PCOUNT  =                    0 / number of group parameters                     
GCOUNT  =                    1 / number of groups                               
TFIELDS =                   53 / number of table fields                         
EXTNAME = 'SINGLE DISH'        / required keyword value                         
NMATRIX =                    1 / one DATA column                                
VELDEF  = 'FREQ-OBS'           / raw receiver frequency                         
PROJID  = 'TAMS    '                                                            
OBSERVER= 'Horiuchi'                                                            
TIMESYS = 'UTC     '           / DSN standard time                              
TELESCOP= 'DSS-43  '                                                            
SITELONG=        211.019942862 / degrees west of Greenwich                      
SITELAT =        -35.403983527 / degrees                                        
SITEELEV=              688.867 / meters                                         
OBSGEO-X=         -4460894.585 / meters                                         
OBSGEO-Y=          2682361.554 / meters                                         
OBSGEO-Z=          -3674748.58 / meters                                         
BACKEND = 'SAO spectrometer'                                                    
MAXIS1  =                32768 / length of DATA axis 1                          
FREQRES =                  0.0                                                  
CTYPE1  = 'FREQ-OBS'           / channel frequency in telescope frame           
CTYPE2  = 'RA---GLS'           / RA J2000.0                                     
MAXIS2  =                    1                                                  
CTYPE3  = 'DEC--GLS'           / decl. J2000                                    
MAXIS3  =                    1                                                  
CTYPE4  = 'STOKES  '           / polarization code: -8 ... 4                    
MAXIS4  =                    2                                                  
CTYPE5  = 'TIME    '           / Python time.time() value                       
MAXIS5  =                    7                                                  
CTYPE6  = 'BEAM    '           / beam 1 or 2                                    
MAXIS6  =                    2                                                  
...                  
TTYPE53 = 'SPECTRUM'                                                            
TFORM53 = '917504E '                                                            
TDIM53  = '(32768,1,1,2,7,2)'                                                       
"""
# temporary for astropy and numpy mismatch
import IPython
IPython.version_info = IPython.release.version.split('.')
IPython.version_info.append('')

import glob
import logging
import numpy
import os
import pyfits
import time

from pylab import *
from matplotlib.font_manager import FontProperties

# import from local modules
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from Data_Reduction.SLATool import SessionAnalyzer
from Data_Reduction.FITS.SDFITSplotter import DSNFITSplotter
from Radio_Astronomy import rms_noise
from support.lists import unique
from support.logs import initiate_option_parser, init_logging, get_loglevel
      
fontP = FontProperties()
fontP.set_size('x-small')      

def show_sources(examiners):
  """
  """
  sources = []
  keys = examiners.keys()
  for key in keys:
    sources += examiners[key].get_sources()
  return unique(sources)
  
def get_average(examiners, source='67P_CG_201'):
  """
  """
  rowfmt = "%03d %5s   %d  %5.1f  %6.1f  %6.4f  %6.4f   %4.2f"
  first_spectrum = {0:True, 1:True}
  sum_y = {0:0, 1:0}
  sum_Tsys = {0:0, 0:1}
  sum_intgr = {0:0, 1:0}
  for exkey in examiners.keys():
    ex = examiners[exkey]
    print "FITS file", os.path.basename(ex.hdulist.filename())
    for tbkey in ex.tables.keys():
      tb = ex.tables[tbkey]
      table_source = tb.sources[0] # for TAMS datasets
      if source in table_source:
        rows = tb.get_rows('OBJECT', table_source)
        x, y, rms, Tsys, intgr = tb.reduce_line(rows=rows)
        print    "                              r.m.s. noise"
        print    "DOY  time  Pol  Tsys   int.   meas'd  expect  ratio"
        print    "--- -----  --- -----  ------  ------  ------  -----"
        for polkey in range(2):
          exp_rms = rms_noise(Tsys[polkey], 1020e6/32768., intgr[polkey])
          ratio = rms[polkey]/exp_rms
          print rowfmt % (tb.DOY, tb.timestr, polkey+1, Tsys[polkey], 
                          intgr[polkey], rms[polkey], exp_rms, ratio)
          if first_spectrum[polkey] and rms[polkey] != numpy.nan:
            sum_y[polkey] = intgr[polkey]*y[polkey]
            sum_Tsys[polkey] = Tsys[polkey]*intgr[polkey]
            sum_intgr[polkey] = intgr[polkey]
            len_x = len(x)
            first_spectrum[polkey] = False
          elif rms[polkey] != numpy.nan:
            if len(x) != len_x:
              print "X array size mismatch for pol",str(polkey+1)
              continue
            sum_y[polkey] += intgr[polkey]*y[polkey]
            sum_Tsys[polkey] += Tsys[polkey]*intgr[polkey]
            sum_intgr[polkey] += intgr[polkey]
      else:
        print source,"is not in table"
  for polkey in range(2):
    sum_y[polkey] /= sum_intgr[polkey]
    sum_Tsys[polkey] /= sum_intgr[polkey]
  return x, sum_y, sum_Tsys, sum_intgr

def figure_rows_and_columns(nspectra):
  """
  """
  cols = int(sqrt(nspectra))
  rows = int(ceil(nspectra/rows))
  return rows, cols
  
  
examples = """
Examples
========
Display diagnostic spectra from SDFITS files; prompt for session
  run interactive.py
Process a specific observing session.
  run interactive.py --console_loglevel=debug --DSS=43 --project=ISM_RRL 
                        --date=2016/013
"""
if __name__ == "__main__":
  p = initiate_option_parser(__doc__, examples)
  p.usage='data_summary.py [options] [data_directory]'

  p.add_argument('--date',
               dest = 'date',
               type = str,
               default = None,
               help = 'Date of observation as YEAR/DOY string')
  p.add_argument('-D', '--dss',
               dest = 'dss',
               type = int,
               default = None,
               help = 'DSN station number')
  p.add_argument('-n', '--name_pattern',
              dest = 'name_pattern',
              type = str,
              default = '*.fits',
              help = 'Pattern required in datafile file name. Default: *')
  p.add_argument('-p', '--project',
               dest = 'project',
               type = str,
               default = None,
               help = "Project code")
  args = p.parse_args()

  logging.basicConfig(level=logging.INFO)
  mylogger = logging.getLogger()
  # change default console logging
  if get_loglevel(args.console_loglevel) > get_loglevel("info"):
    loglevel = "info"
  else:
    loglevel = args.console_loglevel
  init_logging(mylogger,
                 loglevel = get_loglevel(args.file_loglevel),
                 consolevel = get_loglevel(loglevel),
                 logname = args.logpath+"FITS_interactive.log")

  if args.date:
    yearstr, doystr = args.date.split('/')
    year = int(yearstr)
    DOY = int(doystr)
  else:
    year = None
    DOY = None
  sa = SessionAnalyzer(project=args.project, dss=args.dss, year=year, DOY=DOY)

