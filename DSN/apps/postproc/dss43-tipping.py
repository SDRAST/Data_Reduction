"""
dss43-tipping: analyzing SAO/TAMS tipping curve data

This needs to be run on the JPL network or VPN to access the server with the
NMC log data.

Creates a dict 'tables' which are TIPPING CURVE extensions derived from the
tipping curve files for the specified date.

Notes
=====
This can be improved by using the support.tunneling module to check for a
tunnel and creating one if need.  For now, a tunnel must be opened to crux.
"""
#import astropy.io.fits as pyfits
import logging
#import math
#import numpy
import os

from pylab import *
#from scipy.optimize import curve_fit

#from Data_Reduction import get_obs_session
from Data_Reduction.TAMS_tipping import TidTipFinder
from Data_Reduction.DSN.SAO.H5_to_FITS import TipTableMaker
#from Data_Reduction.FITS.DSNFITS import FITSfile
#from Data_Reduction.tipping import airmass, fit_tipcurve_data
from local_dirs import projects_dir
from support.logs import initiate_option_parser, init_logging, get_loglevel

logger = logging.getLogger(__name__)

examples = """
Examples
========
fill in
"""
if __name__ == "__main__":
  p = initiate_option_parser(__doc__, examples)
  p.usage='dss43-tipping.py [options]'

  p.add_argument('--date',
               dest = 'date',
               type = str,
               default = None,
               help = 'Date of observation as YEAR/DOY string')
  p.add_argument('-d', '--dss',
               dest = 'dss',
               type = int,
               default = 43,
               help = 'DSN station number')
  p.add_argument('-p', '--project',
               dest = 'project',
               type = str,
               default = '67P',
               help = "Project code")
  args = p.parse_args()

  mylogger = logging.getLogger()
  # change default console logging
  if get_loglevel(args.console_loglevel) > get_loglevel("info"):
    loglevel = "info"
  else:
    loglevel = args.console_loglevel
  init_logging(mylogger,
                 loglevel = get_loglevel(args.file_loglevel),
                 consolevel = get_loglevel(loglevel),
                 logname = args.logpath+"TidTipAnalyzer.log")
  mylogger.debug("dss43-tipping args: %s", args)
  
  
  # object for analyzing tipping data
  finder = TidTipFinder(date=args.date)
    
  # object for creating table
  ptm = TipTableMaker(args.dss, finder.year, finder.DOY, finder.datefiles)
  
  # make tipping curve tables
  tables = {}
  tbindex = 0
  for datafile in finder.datefiles:
    unixtime = float(os.path.basename(datafile)[12:-4])
    el, Tsys = finder.get_tipping_data(datafile)
    tables[tbindex] = ptm.make_table(unixtime, el, Tsys)
    tbindex += 1
    

  
