# -*- coding: iso-8859-15 -*-
"""
Reads output from WVSR post-processing to make spectroscopy datasets

The columns in the .fft files are::
  freq      - frequency in Hz
  IF1-ps    - pol 1 power
  IF2-ps    - pol 2 power
  IF1-phase - pol 1 voltage phase
  IF2-phase - pol 2 voltage phase
  I         - Stokes I, which should be the same as the mean of
              IF1-ps and IF2-ps
  Q         - Stokes Q, which for linearly polarized radiation should be 
              ±(IF1-ps - IF2-ps) for the polarized component if any. The sign
              can be used to decided which is pol X (or V) and which pol Y
              (or H) if the value is non-zero.
  U         - Stokes U
  V         - Stokes V
  P         - polarization angle in radians
  count     - number of seconds of data (number of spectra averaged)?
  index     - probably meaningless

The spectral line dataset structure looks like this::
In [5]: specdata.keys()
Out[5]: ['header', 'chan 02', 'chan 01']
In [6]: specdata['header'].keys()
Out[6]: ['fftdir']
In [7]: specdata['chan 01'].keys()
Out[7]: 
[1, 2, 3, ... 165, 166, 167, 168, 'freq', 'header']
In [8]: specdata['chan 01']['header']
Out[8]: {'bandwidth': 8000000.0, 'bits/samp': 8, 'OBSFREQ': 8309000000.0}
In [9]: specdata['chan 01']['freq']
Out[9]: array([-4000000., -3999940., -3999880., ...,
                3999820.,  3999880.,  3999940.], dtype=float32)
In [10]: specdata['chan 01'][113].keys()
Out[10]: ['header', 'data']
In [11]: specdata['chan 01'][113]['header'].keys()
Out[11]: ['start', 'datafile', 'end', 'source', 'EXPOSURE']
In [12]: specdata['chan 01'][113]['data'].shape
Out[12]: (131072, 1, 1, 4)
"""
import logging
import numpy
import cPickle as pickle

from numpy import loadtxt
from os import symlink
from os.path import basename, exists

from Automation.WVSR import get_Naudet_FFT_dir, parse_scan_files

from Automation.WVSR.stokes import get_source_averages
from Data_Reduction.DSN.WVSR import get_last_scr_file, get_WVSR_parameters

logger = logging.getLogger(__name__)

def get_channel_IDs(metadata):
  """
  get the channel IDs
  """
  chan_keys = []
  first_scan = metadata.keys()[0]
  for key in metadata[first_scan]:
    if type(key) == str:
      if key[:4] == "chan":
        chan_keys.append(key)
  chan_keys.sort()
  if len(chan_keys):
    logger.debug("get_channel_IDs: found %d channels from %s to %s",
                 len(chan_keys), chan_keys[0], chan_keys[-1])
  else:
    logger.debug("get_channel_IDs: found no channels")
  return chan_keys
  
def get_FFT_data(obsdir, datadir):
  """
  Get spectra and metadata from observing and .fft files
  
  @param obsdir : directory with observation scripts
  @type obsdir : str
  """
  fftdir = get_Naudet_FFT_dir(obsdir)
  if exists(fftdir):
    pass
  else:
    raise RuntimeError("get_FFT_data: %s does not exist" % fftdir)
    
  # get the metadata from the scans files
  metadata = parse_scan_files(obsdir)
  if metadata == None:
    return None, None
  metadata['fftdir'] = fftdir
  logger.debug("get_FFT_data: FFT files are in %s", fftdir)
  
  # get the scan numbers
  scannums = metadata.keys()
  scannums.remove('fftdir')
  scannums.sort()
  logger.debug("get_FFT_data: found %d scans from %s to %s",
                 len(scannums), scannums[0], scannums[-1])
  # get channel IDs
  chan_keys = get_channel_IDs(metadata)
  
  # create the data dict for the science data and for the diagnostics
  scandata = {'header': {'fftdir': fftdir}}
  scandata['header']['dtype'] = ['FREQ-OBS', 'RA--GLS', 'DEC--GLS',
                                 ['I', 'Q', 'U', 'V']]
  IFdata   = {'header': {'fftdir': fftdir}}
  IFdata['header']['dtype'] = ['FREQ-OBS', 'RA--GLS', 'DEC--GLS',
                               ['IF1', 'IF2', 'phase1', 'phase2', 'deg_pol']]
  for chan in chan_keys:
    scandata[chan] = {'header': {}}
    IFdata[chan]   = {'header': {}}
    for scan in scannums:
      scandata[chan][scan] = {}
      scandata[chan][scan]['header'] = {}
      scandata[chan][scan]['header']['start']    = metadata[scan]['start']
      scandata[chan][scan]['header']['end']      = metadata[scan]['end']
      scandata[chan][scan]['header']['OBJECT']   = metadata[scan]['source']
      IFdata[chan][scan] = {}
      IFdata[chan][scan]['header'] = {}
      IFdata[chan][scan]['header']['start']    = metadata[scan]['start']
      IFdata[chan][scan]['header']['end']      = metadata[scan]['end']
      try:
        scandata[chan][scan]['header']['datafile'] = metadata[scan][chan]
        IFdata[chan][scan]['header']['datafile'] = metadata[scan][chan]
      except KeyError:
        # no file for this channel
        scandata[chan][scan]['header']['datafile'] = ''
        IFdata[chan][scan]['header']['datafile'] = ''
  # get the metadata from the .scr files
  WVSRdata = {}
  try:
    scrfile = open(get_last_scr_file(obsdir, "LRCP"),"r")
  except:
    for pol in ['LCP', 'RCP']:
      scrfile = open(get_last_scr_file(obsdir, pol),"r")
      scrdata = scrfile.readlines()
      scrfile.close()
      WVSRdata[pol] = get_WVSR_parameters(scrdata)
      logger.debug("get_FFT_data: got WVSR data for %s", pol) 
  # the WVSR metadata should be the same for both pols so ignore RCP
  for chan in chan_keys:
    channum = int(chan[-2:])
    scandata[chan]['header']['OBSFREQ'] = WVSRdata['LCP'][0][channum]
    # The next is a temporary ad hoc fix.
    scandata[chan]['header']['RESTFREQ'] = scandata[chan]['header']['OBSFREQ']
    scandata[chan]['header']['BANDWID'] = WVSRdata['LCP'][1][channum]
    scandata[chan]['header']['bits/samp'] = WVSRdata['LCP'][2][channum]
    IFdata[chan]['header']['OBSFREQ']  = WVSRdata['LCP'][0][channum]
    IFdata[chan]['header']['BANDWIDT'] = WVSRdata['LCP'][1][channum]
    IFdata[chan]['header']['bits/samp'] = WVSRdata['LCP'][2][channum]
    # now get the data
    logger.debug("get_FFT_data: getting data for channel %s", chan)
    # for each scan of that channel
    for scan in scannums:
      thisdata = get_data(scandata, chan, scan)
      if type(thisdata) != numpy.ndarray:
        # no data returned
        logger.warning("get_FFT_data: no data for %s scan %d", chan, scan) 
        continue
      if scandata[chan].has_key('freq'):
        # we already have freq data from previous scan
        pass
      else:
        # first scan with freq data
        scandata[chan]['freq'] = thisdata['freq']
        IFdata[chan]['freq']   = thisdata['freq']
      try:
        scandata[chan][scan]['header']['EXPOSURE'] = thisdata['count']
        IFdata[chan][scan]['header']['EXPOSURE']   = thisdata['count']
        scandata[chan][scan]['data'] = numpy.empty((131072,1,1,4))
        IFdata[chan][scan]['data']   = numpy.empty((131072,1,1,5))
        scandata[chan][scan]['data'][:,0,0,0] = thisdata['I']
        scandata[chan][scan]['data'][:,0,0,1] = thisdata['Q']
        scandata[chan][scan]['data'][:,0,0,2] = thisdata['U']
        scandata[chan][scan]['data'][:,0,0,3] = thisdata['V']
        IFdata[chan][scan]['data'][:,0,0,0] = thisdata['IF1-ps']
        IFdata[chan][scan]['data'][:,0,0,1] = thisdata['IF2-ps']
        IFdata[chan][scan]['data'][:,0,0,2] = thisdata['IF1-phase']
        IFdata[chan][scan]['data'][:,0,0,3] = thisdata['IF2-phase']
        IFdata[chan][scan]['data'][:,0,0,4] = thisdata['P']
        logger.debug("get_FFT_data: got data for channel %s scan %d", chan, scan)
      except KeyError:
        logger.debug("get_FFT_data: no data for channel %s scan %d", chan, scan)
  # save the data
  logger.debug("get_FFT_data: writing Stokes data...")
  pklfile = open(datadir+"specdata.pkl","w")
  pickle.dump(scandata, pklfile)
  pklfile.close()
  logger.debug("get_FFT_data: ...done")
  logger.debug("get_FFT_data: writing VLBI data...")
  pklfile = open(datadir+"vlbidata.pkl","w")
  pickle.dump(scandata, pklfile)
  pklfile.close()
  logger.debug("get_FFT_data: ...done")
  # this should now work for scan 113
  # plot(scandata['chan 01'][113][0], scandata['chan 01'][113][1][:,0,0,0])
  return scandata, IFdata

def read_FFT_file(datafilename):
  """
  """
  try:
    datafile = open(datafilename,"r")
    logger.debug("read_FFT_file: getting data file %s", datafile)
  except IOError, details:
    # Probably missing FFT files
    logger.warning("read_FFT_file: no FFT data for because %s", str(details))
    return None
  firstline = datafile.readline()
  datafile.close()
  #logger.debug("get_data: first line: %s", firstline)
  if firstline == "" or firstline.split()[0] == "Error":
    logger.warning("read_FFT_file: post-processor could not open WVSR file")
    return None
  try:
    data = loadtxt(datafilename, 
                      dtype={'names': ('freq', 'IF1-ps', 'IF2-ps', 'IF1-phase',
                                       'IF2-phase', 'I', 'Q', 'U', 'V', 'P',
                                       'count','index'),
                             'formats': 10*['f4']+2*['i4']})
    logger.debug("read_FFT_file: loaded data from %s", datafilename)
    return data
  except ValueError:
    logger.warning("read_FFT_file: could not convert %s", datafilename)
    return None
  

def get_data(scans, chan, scan):
  """
  This loads data from the FFT files associated with 'scan'
  
  @param scans : scan metadata
  @type  scans : dict
  
  @param scan : scan number
  @type  scan : int
  """
  datafilename = scans['header']['fftdir'] + \
                 scans[chan][scan]['header']['datafile']
  logger.debug("get_data: looking for %s scan %d in %s: %s",
                chan, scan, scans['header']['fftdir'],
                datafilename)
  try:
    datafile = open(datafilename,"r")
    logger.debug("get_data: getting data file %s", datafile)
  except IOError, details:
    # Probably missing FFT files
    logger.warning("get_data: no FFT file for chan %s scan %d because %s",
                   chan, scan, str(details))
    return None
    
  firstline = datafile.readline()
  datafile.close()
  #logger.debug("get_data: first line: %s", firstline)
  if firstline == "" or firstline.split()[0] == "Error":
    logger.warning("get_data: post-processor could not open WVSR file")
    return None
  try:
    data = loadtxt(datafilename, 
                      dtype={'names': ('freq', 'IF1-ps', 'IF2-ps', 'IF1-phase',
                                       'IF2-phase', 'I', 'Q', 'U', 'V', 'P',
                                       'count','index'),
                             'formats': ['f8']+9*['f8']+2*['i4']})
    logger.debug("get_data: loaded data for channel %s", chan)
    return data
  except ValueError:
    logger.warning("get_data: could not convert %s", datafilename)
    return None

def get_sources(scans, noref=True):
  """
  Returns a list of sources and the first scan number
  
  @param scans : metadata for all the scans
  @type  scans : dict
  
  @param noref : ignore the reference position data
  @type  noref :
  
  @return: tuple (list of source names, dict with first scan for each source}
  """
  def add_source_if_new(sources, source, scan, first_scan, last_scan):
    """
    """
    try:
      # source exists in the list
      sources.index(source)
      last_scan[source] = scan
    except ValueError:
      # this is a new source
      logger.debug("get_sources.add_source: adding %s at scan %d", source, scan)
      sources.append(source)
      first_scan[source] = scan
    return sources, first_scan, last_scan
  
  # the sources should be the same in all channels
  scannums = scans['chan 01'].keys()
  scannums.remove('freq')
  scannums.remove('header')
  scannums.sort()
  sources = []
  first_scan = {}
  last_scan = {}
  logger.debug("get_sources: scans: %s", scannums)
  for scan in scannums:
    # sources should be the same in both channels
    try:
      source = scans['chan 01'][scan]['header']['OBJECT']
      logger.debug("get_sources: scan %d is for %s", scan, source)
      if source[-4:] == "-ref":
        last_scan[source[:-4]] = scan
        if noref:
          pass # do not add to list
        else:
          sources, first_scan, last_scan = \
                add_source_if_new(sources, source, scan, first_scan, last_scan)
      else:
        sources, first_scan, last_scan = \
                add_source_if_new(sources, source, scan, first_scan, last_scan)
    except ValueError:
      logger.warning("get_sources: no header data for scan %d", scan)
  # last source if last source is a -ref
  last_scan[source] = scan        
  return sources, first_scan, last_scan

