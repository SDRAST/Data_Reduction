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
import glob
import logging
import pyfits
import time

# import from local modules
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from Data_Reduction.FITS.plot import DSNFITSplotter
#from Data_Reduction.FITS.SDFITSexaminer import DSNFITSexaminer
from support.logs import initiate_option_parser, init_logging, get_loglevel


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
  p.add_argument('-D', '--DSS',
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

  # Get a list of datafiles
  projectdatapath, projworkpath, datapath = get_obs_dirs(
    *get_obs_session(dss=args.dss, project=args.project, date=args.date,
                     datafmt="FITS"))
                                                   
  # get the datafiles to be processed     
  datafiles = glob.glob(datapath+args.name_pattern)
  datafiles.sort()

  header = {}; examiner = {}; n_scans = {}; scan_keys = {}; sttm = {}
  dfindex = 0
  for datafile in datafiles:
    hdulist = pyfits.open(datafile)
    header[dfindex] = hdulist[0].header
    examiner[dfindex]  = DSNFITSplotter(hdulist[1])
    # redefine project work path
    #   eventually this should put output with the correct project but right
    #   now, the 67P project was misnamed at observe time
    project = examiner[dfindex].header['PROJID']
    n_scans[dfindex] = len(examiner[dfindex].dataset)
    scan_keys[dfindex] = range(n_scans[dfindex])
    first_scan = scan_keys[dfindex][0]
    # This will go at the top of each figure
    sttm[dfindex] = time.gmtime(
      examiner[dfindex].dataset[first_scan]['UNIXtime'][0,0,0,0,0,0])
    title = "%4d/%02d/%02d (%03d) %02d:%02d" % (sttm[dfindex].tm_year,
                                                sttm[dfindex].tm_mon,
                                                sttm[dfindex].tm_mday,
                                                sttm[dfindex].tm_yday,
                                                sttm[dfindex].tm_hour,
                                                sttm[dfindex].tm_min)
    source = examiner[dfindex].dataset['OBJECT'][first_scan]
    mylogger.info("Dataset %d is for project %s on %s with source %s",
                  dfindex, project, title, source)
    dfindex += 1

