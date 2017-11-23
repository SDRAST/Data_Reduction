import calendar
import glob
import logging
import math
import numpy
import os
import time

from pylab import *

from Data_Reduction import get_obs_session
from DatesTimes import calendar_date
from support.logs import initiate_option_parser, init_logging, get_loglevel

logger = logging.getLogger(__name__)

class TidTipAnalyzer(object):
  """
  """
  def __init__(self,
               directory="/home/kuiper/mnt/crux/home/ops/roach_data/sao_test_data/log_dir/",
               date="2015/204"):
    """
    """
    self.logger = logging.getLogger(logger.name+".TidTipAnalyzer")
    if date:
      yearstr, doystr = args.date.split('/')
      year = int(yearstr)
      DOY = int(doystr)
      year,month,day = calendar_date(year, DOY)
      UNIXdate = calendar.timegm((year,month,day, 0,0,0, 0,0))
      self.dir = directory
      if year in [2014,2015,2016]:
        self.dir += "2014-2016/"
      pattern = self.dir+"tipping_data"+str(UNIXdate)[:5]+"*"
      filtered_files = glob.glob(pattern)
      self.am = {}
      self.tsys = {}
      index = 0
      for f in filtered_files:
        self.logger.debug("__init__: opening %d: %s", index, os.path.basename(f))
        data = numpy.load(f)
        self.am[index] = 1/numpy.sin(math.pi*data[0].astype(float)/180)
        self.tsys[index] = numpy.concatenate(data[1], axis=1).transpose()
        index += 1
      

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

  mylogger = logging.getLogger()
  mylogger.setLevel(logging.DEBUG)
  
  analyzer = TidTipAnalyzer(date=args.date)
  for tip in analyzer.am.keys():
    for IF in range(4):
      plot(analyzer.am[tip], analyzer.tsys[tip][:,IF], '.', label="tip %d IF %d" % (tip,IF))
  legend(loc="upper left", numpoints=1)
  grid()
  title(args.date)
  xlabel("airmass")
  show()
