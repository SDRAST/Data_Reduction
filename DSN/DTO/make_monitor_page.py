"""
make session index page

At the moment this only displays the results for the I channel because the
Q channel gets its data from the same input
"""
import cPickle
import glob
import h5py
import logging
import os
import sys

from Data_Reduction import get_obs_dirs
from Data_Reduction.DSN import path_to_remote
from support.lists import unique
from support.logs import get_loglevel, initiate_option_parser, init_logging

logger = logging.getLogger(__name__)

def make_grid(title="Table Title", files=None):
  """
  grid thumbnails
  """
  html = "<CENTER><H2>"+title+"</H2></CENTER>\n"
  html += "<TABLE>\n"
  nrows = len(files)/4
  if nrows % 4:
    nrows +=1
  for row in range(nrows):
    html += "<TR>"
    for col in range(4):
      index = 4*row + col
      if index < len(files):
        filename = os.path.basename(files[index])
        html += "<TD><A HREF='"+filename+"'><IMG SRC='thumbnails/"+filename+"'</TD>"
    html += "\n</TR>"
  html += "</TABLE>\n"
  return html

def find_good_signals(datapath, session_path):
  """
  select I inputs only
  
  This takes a minute or so to load the data from the pkl file.
  
  @param filename : name of a session's 'spectra*' file
  @type  filename : str
  """
  logger.debug("find_good_signals: in %s", datapath+"mon-*.hdf5")
  session_data = glob.glob(datapath+"mon-*.hdf5")
  logger.debug("find_good_signals: found %s", session_data)
  good_signals = []
  for filename in session_data:
    data = h5py.File(filename)
    if data.attrs['channel'] == "I":
      logger.debug("find_good_signals: signals are %s", data.attrs['signal'])
      good_signals.append(data.attrs['signal'])
    data.close()
  return good_signals

if __name__ == "__main__":
  logging.basicConfig(level=logging.DEBUG)
  mylogger=logging.getLogger()

  p = initiate_option_parser("Kurtosis monitor","")
  p.usage = "python kurtspec_monitor.py <kwargs>"
  # Add other options here
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
  p.add_argument('-p', '--project',
                 dest = 'project',
                 type = str,
                 default = "PESD",
                 help = "Project code")
  p.add_argument('-u', '--use_only',
                 dest = 'use_only',
                 type=str,
                 default = None,
                 help = 'Use only')
  p.add_argument('-w', '--workstation',
                 dest = 'workstation',
                 type = str,
                 default = None,
                 help = "Workstation for post-processing")
  
  args = p.parse_args()
  
  # This cannot be delegated to another module or class
  mylogger = init_logging(logging.getLogger(),
                          loglevel   = get_loglevel(args.file_loglevel),
                          consolevel = get_loglevel(args.console_loglevel),
                          logname    = args.logpath+__name__+".log")
  mylogger.debug("arguments: %s", args)
  mylogger.debug(" Handlers: %s", mylogger.handlers)

  if args.date:
    yearstr, doystr = args.date.split('/')
    year = int(yearstr)
    doy = int(doystr)
  else:
    mylogger.error(" 'date' is a required argument")
    sys.exit(1)
  if args.dss:
    pass
  else:
    mylogger.error(" 'dss' is a required argument")
    sys.exit(1)
  if args.workstation:
    pass
  else:
    mylogger.error(" 'workstation' is a required argument")
    sys.exit(1)
  
  # get the session path
  projectdatapath, projworkpath, datapath = \
                    get_obs_dirs("PESD", args.dss, year, doy)
  mylogger.debug("project data path: %s", projectdatapath)
  mylogger.debug("project work path: %s", projworkpath)
  datapath = path_to_remote(args.workstation, projectdatapath)
  mylogger.debug("data path: %s", datapath)
  sessionpath = path_to_remote(args.workstation, projworkpath)
  
  # select signals to display
  if args.use_only:
    good_signals = args.use_only.split(',')
  else:
    good_signals = find_good_signals(datapath, projworkpath)
  
  # extract the signal information
  #    passband files
  pbfiles = glob.glob(sessionpath+"thumbnails/passband*.png")
  pbfiles.sort()
  #    get observed signals
  signames = []
  for filename in pbfiles:
    signame = os.path.basename(filename).split('-')[1]
    if signame in good_signals:
      signames.append(signame)
  signames = unique(signames)
  signames.sort()
  mylogger.debug("signals for this session: %s", signames)
  
  # page header
  html = "<HTML>\n"
  html += "<CENTER><H1>Session Summary for %s on %4d/%03d</H1></CENTER>\n" % \
         (args.project, year, doy)
  for signame in signames:
    html += "<HR/>" 
    #   get the images
    pfiles = glob.glob(sessionpath+"thumbnails/specgram-power-"+signame+"*.png")
    pfiles.sort()
    mylogger.debug("found: %s", pfiles)
    pbfiles = glob.glob(sessionpath+"thumbnails/passband-"+signame+"*.png")
    pbfiles.sort()
    mylogger.debug("found: %s", pbfiles)
    kfiles = glob.glob(sessionpath+"thumbnails/specgram-kurtosis-"+signame+"-*.png")
    kfiles.sort()
    mylogger.debug("found: %s", kfiles)
    # get the spectra
    ktfiles = glob.glob(sessionpath+"thumbnails/kurtosis-"+signame+"*.png")
    ktfiles.sort()
  
    # generate the signal images
    band, dss, pol = signame[0],signame[1:3],signame[3]
    if band == "K":
      # refine this to include K-hi
      band == "Ku"
    html += "<CENTER>\n"
    html += "<H2>DSS-%2s %s-band %sCP for %s on %4d/%03d</H1>\n" % \
            (dss, band, pol, args.project, year, doy)
    html += "<P><A HREF='PESD/background'>Background</A></P>\n"
    html += "<P><A HREF='PESD/aboutKurtosis'>About Kurtosis</A></P>\n"
    html += "<P>Click on any image for full resolution.</P>\n"
    html += "<TABLE BORDER=1>\n"
    
    html += "<TR>\n"
    html += "<TD>\n"
    html += make_grid(title=band+"-band Kurtosis", files=kfiles)
    html += "</TD>\n"
    html += "</TR>\n"
      
    html += "<TR>\n"
    html += "<TD>\n"
    html += make_grid(title=band+"-band Frequency Averaged Kurtosis", files=ktfiles)
    html += "</TD>\n"
    html += "</TR>\n"
      
    html += "<TR>\n"
    html += "<TD>\n"
    html += make_grid(title=band+"-band Time Averaged Passbands", files=pbfiles)
    html += "</TD>\n"
    html += "</TR>\n"
      
    html += "<TR>\n"
    html += "<TD>\n"
    html += make_grid(title=band+"-band Total Power", files=pfiles)
    html += "</TD>"
    html += "</TR>\n"
    
    html += "</CENTER>\n"
    html += "</TABLE>\n"
  
  html += "</HTML>\n"
  if os.path.exists(sessionpath):
    pass
  else:
    os.makedirs(sessionpath)
  fullpath = sessionpath+"index.html"
  htfile = open(fullpath,"w")
  htfile.write(html)
  htfile.close()
  

