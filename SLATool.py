"""
Spectral Line Analysis Tool

Environment for analyzing spectral line from multiple observatories. Manages
various data reduction programs
"""
import astropy.io.fits as pyfits
import copy
import glob
import logging
import numpy
import os
import os.path
import time

import Data_Reduction as DR
import Data_Reduction.FITS.SDFITSexaminer as FITSex
import Data_Reduction.tipping as DRtip
import Radio_Astronomy

polnames = ["L", "R"]

logger = logging.getLogger(__name__)

class Spectrum(object):
  """
  Generic container for a spectrum and its metadata
  
  Public Attributes::
    * bandwidth
    * frame
    * origin
    * spectrum - a 2D ndarray (2 pols, N pixels)
    * obsfreq
    * refpix
  """
  def __init__(self):
    """
    """
    pass
  
class SessionAnalyzer(object):
  """
  Tool for reducing multiple data reduction sessions
  
  Attributes::
    datapath
    DOY
    DSS
    examiners
    logger
    project
    projectdatapath
    projworkpath
    year
    
  Example::
    In [1]: from Data_Reduction.SLATool import SessionAnalyzer
    In [2]: sa = SessionAnalyzer(project='67P', year=2015, DOY=204)
    In [3]: x, sum_y, sum_Tsys, sum_intgr = sa.get_average()
  """
  def __init__(self, project=None, dss=None, year=None, DOY=None):
    """
    initiate a SessionAnalyzer
    
    @param project : name as in /usr/local/projects directory
    @type  project : str
    
    @param dss : DSN station
    @type  dss : int
    
    @param year : of observing session
    @type  year : int
    
    @param DOY : of observing session
    @type  DOY : int
    """
    self.logger = logging.getLogger(logger.name+".SessionAnalyzer")
    
    # get the session details
    #if year and DOY:
    #  date = "%4d/%03d" % (year, DOY)
    #  date = None
    #  self.project, self.DSS, self.year, self.DOY = DR.get_obs_session(dss=dss,
    #                                                project=project, date=date)
    #else:
    if True:
      self.project = project
      self.DSS = dss
      self.year = year
      self.DOY = DOY
    self.projectdatapath, self.projworkpath, self.datapath = \
      DR.get_obs_dirs(self.project, self.DSS, self.year, self.DOY, datafmt="FITS")
    # get the datafiles to be processed     
    datafiles = glob.glob(self.datapath+"*.fits")
    if datafiles:
      datafiles.sort()
      self.logger.info("__init__: found %d datafiles %4d/%03d",
                       len(datafiles), self.year, self.DOY)
      self.examiners = self.open_datafiles(datafiles)
      if self.examiners:
        self.get_sources()
    else:
      self.logger.warning("__init__: no datafiles for DSS-%2d on DOY %4d/%03d",
                          self.DSS, self.year, self.DOY)
      self.examiners = None
  
  def open_datafiles(self, datafiles):
    """
    Opens the files for the session in /usr/local/RA_data/FITS
    """
    dfindex = 0
    examiners = {}
    for datafile in datafiles:
      self.logger.info("open_datafiles:opening %s", os.path.basename(datafile))
      examiners[dfindex] = FITSex.DSNFITSexaminer(parent=self, FITSfile=datafile)
      for table_key in list(examiners[dfindex].tables.keys()):
        examiners[dfindex].tables[table_key].report_table()
      dfindex += 1
    if examiners:
      self.examiner_keys = list(examiners.keys())
      self.examiner_keys.sort()
      self.logger.info("open_datafiles: started %d examiners" % 
                     len(list(examiners.keys())))
    else:
      self.logger.warning("open_datafiles: no files found")
    return examiners

  def make_session_names(self):
    """
    creates the directory path to the data from the session
    
    It is crafted from the name of the first FITS file
    """
    first_ex = self.examiners[self.examiner_keys[0]]
    sttm = time.gmtime(first_ex.tables[0].get_first_good_value('UNIXtime'))
    figtitle = "%4d/%02d/%02d (%03d)" % (sttm.tm_year, sttm.tm_mon,
                                         sttm.tm_mday, sttm.tm_yday)
    savename = "_".join(os.path.splitext(os.path.basename(first_ex.file))[0].split('-')[:2])
    savepath = os.path.dirname(first_ex.file)+"/"+savename
    return figtitle, savepath

  def get_sources(self):
    """
    Gets a list of sources from all the datasets in the session
    """
    self.sources = {}
    keys = list(self.examiners.keys())
    for key in keys:
      source_names = self.examiners[key].get_sources()
      for name in source_names:
        if name in self.sources:
          self.sources[name].append(key)
        else:
          self.sources[name] = [key]
    return list(self.sources.keys())

  def get_good_weather_data(self, examiner_keys=None):
    """
    Get data from specific datasets for analyzing environmental conditions.
      
    Returns a dict with keys 'TAMBIENT', 'WINDDIRE', 'UNIXtime', 'TSYS',
    'HUMIDITY', 'PRESSURE', 'ELEVATIO', 'WINDSPEE'.  The data asociated with
    each key is a dict with numpy array for (SIG state) True and for False.
    The 'TSYS' array has four axes representing::
    
      time index   - 0-based sequence in order of matplotlib datenum 
      subchannel   - CYCLE value
      beam         - 1-based number sequence
      IF           - 1-based number sequence, usually representing pol
      
    The other keys have only a time axis.
    """
    good_data = {}
    if examiner_keys:
      pass
    else:
      examiner_keys = self.examiner_keys
      
    for exkey in examiner_keys:
      examiner = self.examiners[exkey]
      self.logger.debug("get_good_weather_data: data file is %s",examiner.file)
      # assume multiple tables are in time order
      tablekeys = list(examiner.tables.keys())
      tablekeys.sort
      self.logger.debug("get_good_weather_data: table keys %s", tablekeys)
      for tablekey in tablekeys:
        table = examiner.tables[tablekey]
        wx_data = table.get_wx_datacubes()
        param_keys = list(wx_data.keys())
        self.logger.debug("get_good_weather_data: parameters: %s", param_keys)
        for param_key in param_keys:
          for switch_state in [True, False]:
            # the state may already be defined so check  
            if param_key in good_data:
              # parameter already exists in the output data
              if switch_state in wx_data[param_key]:
                # this state exists in the input data
                if switch_state in good_data[param_key]:
                  # this state also exists in the output data so append
                  good_data[param_key][switch_state] = \
                         numpy.append(good_data[param_key][switch_state],
                                      wx_data[param_key][switch_state], axis=0)
                  self.logger.debug("get_good_weather_data: add %s for SIG=%s",
                                    param_key, switch_state)
                else:
                  # this state does not exist in the output so create
                  good_data[param_key][switch_state] = \
                                               wx_data[param_key][switch_state]
                  self.logger.debug(
                               "get_good_weather_data: create state %s for %s",
                               switch_state, param_key)
              else:
                # state does not exist in the input data
                pass
            else:
              # parameter not defined in the output data so create
              self.logger.debug("get_good_weather_data: got %s for SIG=%s",
                                param_key, switch_state)
              good_data[param_key] = \
                               {switch_state: wx_data[param_key][switch_state]}
    return good_data

  def fit_Tsys_to_airmass(self, weather_data=None, examiner_keys=None,
                          Tatm=250, linear=True):
    """
    Fit tipping curve data implicit in sessions elev and Tsys data
    
    Returns numpy arrays with indices for sig/ref state, subchannel, beam, IF.
    The first returned value is the zero airmass intercept and the second its
    standard deviation.  The second returned value is K/airmass
    
    @param weather_data : consolidated environmental data
    @type  weather_data : dict
    
    @param examiner_keys : keys of files from this date to be included
    @type  examiner_keys : list of int
    
    @param Tatm : air temperature along line of sight
    @type  Tatm : float
  
    @param linear : use the linear (low tau) approximation
    @type  linear : True
    """
    # get the data
    if weather_data:
      pass
      self.logger.debug("fit_Tsys_to_airmass: using provided data")
    elif examiner_keys:
      self.logger.debug("fit_Tsys_to_airmass: extracting data from tables %s",
                        examiner_keys)
      weather_data = self.get_good_weather_data(examiner_keys=examiner_keys)
    else:
      self.logger.debug("plot_elev_and_Tsys: using all data for this date")
      weather_data = self.get_good_weather_data()
    
    # get Tsys data structure  
    if 'TSYS' in weather_data:
      nrows, num_cy, num_bm, num_pl = weather_data['TSYS'][True].shape
    else:
      self.logger.error("fit_Tsys_to_airmass: no system temperature data")
      return None
    states = list(weather_data['TSYS'].keys())
    num_st = len(states)
    
    if 'ELEVATIO' in weather_data:
      pass
    else:
      self.logger.error("fit_Tsys_to_airmass: no elevation data")
      return None
      
    # fit the data
    param_shape = (num_st, num_cy, num_bm, num_pl)
    Trx = numpy.zeros(param_shape); sigTrx = numpy.zeros(param_shape)
    tau = numpy.zeros(param_shape); sigtau = numpy.zeros(param_shape)
    for sig in states:              # sig/ref state first
      for subch in range(num_cy):   # subchannel second
        for beam in range(num_bm):  # beam third
          for pol in range(num_pl): # pol fourth
            el = weather_data['ELEVATIO'][sig]
            Tsys = weather_data['TSYS'][sig][:,subch,beam,pol]
            Trx[sig,subch,beam,pol], sigTrx[sig,subch,beam,pol], \
            tau[sig,subch,beam,pol], sigtau[sig,subch,beam,pol] \
            = DRtip.fit_tipcurve_data(el, Tsys, Tatm=Tatm, linear=linear)
    return Trx, sigTrx, tau, sigtau
  
  def Tsys_scaling(self, weather_data=None, examiner_keys=None, Tatm=250):
    """
    factors to rescale the system powers or temperatures to K based on FE noise
    
    @param weather_data : consolidated environmental data
    @type  weather_data : dict
    
    @param examiner_keys : keys of files from this date to be included
    @type  examiner_keys : list of int
    
    @param Tatm : air temperature along line of sight
    @type  Tatm : float
    """
    # use default linear approximation
    Prx, sigPrx, Ppam, sigPpam = self.fit_Tsys_to_airmass(
                                                   weather_data=weather_data,
                                                   examiner_keys=examiner_keys,
                                                   Tatm=Tatm)
    # estimate gain (or Tsys correction)
    num_st, num_cy, num_bm, num_pl = Prx.shape
    param_shape = (num_st, num_cy, num_bm, num_pl)
    gain = numpy.zeros(param_shape)    # K per powerunit
    sigT = numpy.zeros(param_shape)    # K
    Kpam = numpy.zeros(param_shape)    # K per airmass
    sigKpam = numpy.zeros(param_shape)
    for sig in range(num_st):       # sig/ref state first
      for subch in range(num_cy):   # subchannel second
        for beam in range(num_bm):  # beam third
          for pol in range(num_pl): # pol fourth
            # assume same FE for all datasets in the session
            Trx = self.examiners[0].tables[0].FE.Tsys_vacuum(beam=beam+1,
                                                             pol=polnames[pol])
            G = Trx/Prx[sig,subch,beam,pol]
            gain[sig,subch,beam,pol]    = G
            sigT[sig,subch,beam,pol]    = G*sigPrx[sig,subch,beam,pol]
            Kpam[sig,subch,beam,pol]    = G*Ppam[sig,subch,beam,pol]
            sigKpam[sig,subch,beam,pol] = G*sigPpam[sig,subch,beam,pol]
    return gain, sigT, Kpam, sigKpam
  
  def rescale_Tsys(self, gain, siggain, Kpam, sigKpam):
    """
    rescale system temperatures from tipping curve fits
    """
    
    for exkey in list(self.examiners.keys()):
      for tbkey in list(self.examiners[exkey].tables.keys()):
        tb = self.examiners[exkey].tables[tbkey]
        Tsys = tb.data['TSYS']
  
  def get_average(self, source='67P_CG_201'):
    """
    Computes average spectrum for a source in the session
    
    Prints r.m.s. noise for each dataset and all datasets together
    
    @param source : source for which averaging is done
    """
    rowfmt = "%03d %5s   %d  %5.1f  %6.1f  %6.4f  %6.4f   %4.2f"
    first_spectrum = {0:True, 1:True}
    sum_y = {0:0, 1:0}
    sum_Tsys = {0:0, 0:1}
    sum_intgr = {0:0, 1:0}
    for exkey in list(self.examiners.keys()):
      ex = self.examiners[exkey]
      print(("FITS file", os.path.basename(ex.file)))
      for tbkey in list(ex.tables.keys()):
        tb = ex.tables[tbkey]
        table_source = tb.sources[0] # for TAMS datasets
        if source in table_source:
          rows = tb.get_rows('OBJECT', table_source)
          x, y, rms, Tsys, intgr = tb.reduce_line(rows=rows)
          print ("                              r.m.s. noise")
          print ("DOY  time  Pol  Tsys   int.   meas'd  expect  ratio")
          print ("--- -----  --- -----  ------  ------  ------  -----")
          for polkey in range(2):
            exp_rms = Radio_Astronomy.rms_noise(Tsys[polkey], 1020e6/32768., intgr[polkey])
            ratio = rms[polkey]/exp_rms
            print(rowfmt % (tb.DOY, tb.timestr, polkey+1, Tsys[polkey], 
                            intgr[polkey], rms[polkey], exp_rms, ratio))
            if first_spectrum[polkey] and rms[polkey] != numpy.nan:
              sum_y[polkey] = intgr[polkey]*y[polkey]
              sum_Tsys[polkey] = Tsys[polkey]*intgr[polkey]
              sum_intgr[polkey] = intgr[polkey]
              len_x = len(x)
              first_spectrum[polkey] = False
            elif rms[polkey] != numpy.nan:
              if len(x) != len_x:
                print(("X array size mismatch for pol",str(polkey+1)))
                continue
              sum_y[polkey] += intgr[polkey]*y[polkey]
              sum_Tsys[polkey] += Tsys[polkey]*intgr[polkey]
              sum_intgr[polkey] += intgr[polkey]
        else:
          print((source,"is not in table"))
    for polkey in range(2):
      sum_y[polkey] /= sum_intgr[polkey]
      sum_Tsys[polkey] /= sum_intgr[polkey]
    return x, sum_y, sum_Tsys, sum_intgr

