"""
Opens SDFITS file for analysis
"""
import glob
import logging
import pyfits
import time

# import from local modules
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from support.logs import initiate_option_parser, init_logging, get_loglevel
from SDFITSexaminer import DSNFITSexaminer

examples = """
Examples
========
Display diagnostic spectra from SDFITS files; prompt for session
  run interactive.py
Process a specific observing session.
  run interactive.py --console_loglevel=debug --DSS=43 --project=RRL 
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

  header = {}; examiner = {}; n_scans = {}; scan_keys = {}; sttm = {}; dfindex = 0
  for datafile in datafiles:
    hdulist = pyfits.open(datafile)
    header[dfindex] = hdulist[0].header
    dss = int(header[dfindex]['TELESCOP'].split('-')[1])
    examiner[dfindex]  = DSNFITSexaminer(hdulist[1], dss=dss)
    # redefine project work path
    #   eventually this should put output with the correct project but right
    #   now, the 67P project was misnamed at on=bserve time
    project = examiner[dfindex].header['PROJID']
    #parts = projworkpath.split('/')
    #parts[4] = project
    #projworkpath = '/'.join(parts)
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

