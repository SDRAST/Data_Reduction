# -*- coding: utf-8 -*-
"""
Functions for reducing Solar Patrol data

A lot of modules under this need fixing up because of change from package 
Observatory to MonitorControl and algorithmic to object-oriented programming.
"""
import datetime
import glob
import numpy
import os.path
# numpy has a function 'copy'.  Override it.
import copy
import matplotlib.dates
import logging

import Astronomy
import support
import Data_Reduction as DR
import DatesTimes as DT
#from Data_Reduction.GAVRT.solar import create_metadata_sheet

logger = logging.getLogger(__name__)

class Observation(DR.Observation):
  """
  DSS-28 observations from before the ``tlog`` database table became available.
  """
  def __init__(self, name=None, dss=28, date=None, project='SolarPatrol'):
    """
    initialize observation from t-files
  
    @param date : required YYYY/DDD
    @type  date : str
  
    @param dss : optional station number; default 28
    @type  dss : int
    
    @param channel_names : optional subset of channels to be processed
    @type  channel_names : list of str
    """
    self.logger = logging.getLogger(logger.name+".Observation")
    DR.Observation.__init__(self, name=date, date=date, dss=dss, 
                            project=project)
    
    #self.obs =Astronomy.Ephem.DSS(dss)
    #y,d = date.split('/')
    #self.year = int(y); self.DOY = int(d)
    #projdatapath, self.sessionpath, rawdatapath = \
    #                          DR.get_obs_dirs(project, dss, self.year, self.DOY,
    #                                          datafmt=None)
  
  def extended_init(self, channels_names=None):
    datapath = self.sessionpath+'t'+str(self.year)[-2:]+("%03d" % self.DOY)
    if channel_names:
      pass
    else:
      channel_names = self.get_channel_names(datapath)
    channel_names.sort()
    self.logger.debug("processing channels %s", channel_names)
    # this has created Channel objects which we now want to treat as
    # superclasses for the Channel class defined here
    for key in list(self.channel.keys()):
        del(self.channel[key])
        self.channel[key] = self.Channel(self, key)
    self.get_obs_times(datapath)
    self.logger.info("loading data from files; this takes a while")
    self.data = {'counts': {}}
    for key in list(self.channel.keys()):
        self.channeldata = self.channel[key].load_tfile_data(datapath, dss=dss,
                                              start=self.start,  stop=self.end)
        if 'unixtime' in self.data:
          # then the basic data have read in from the first channel
          pass
        else:
          self.data['unixtime'] = self.channeldata[0]['Epoch']
          self.data['az'] = self.channeldata[0]['Az']
          self.data['el'] = self.channeldata[0]['El']
          self.data['diode'] = self.channeldata[0]['Diode']
          self.data['integr'] = self.channeldata[0]['Int']
        self.data['counts'][key] = self.channeldata[0]['Tsys']
        # del(self.channeldata)  
 
    
  def get_channel_names(self, datapath):
    """
    """
    self.logger.debug("get_channel_names: for %s", datapath)
    names = glob.glob(datapath+".*")
    self.logger.debug("get_channel_names: from %s", names)
    channel_names = []
    for name in names:
      channel_names.append("%02d" % int(os.path.splitext(name)[1][1:]))
    return channel_names
  
  def get_obs_times(self, datapath):
    """
    returns: (datetime.datetime, datetime.datetime)
    """
    # first get start and stop times
    starts = []; ends = []
    for channel in list(self.channel.keys()):
      filename = datapath+"."+str(int(channel))
      self.logger.debug("get_channel times: from %s", filename)
      fd = open(filename,'r')
      lines = fd.readlines()
      fd.close()
      starts.append(float(lines[1].split()[3]))
      ends.append(float(lines[-1].split()[3]))
    start = DT.UnixTime_to_datetime(numpy.array(starts).max())
    end   = DT.UnixTime_to_datetime(numpy.array(ends).min())
    self.start = start.replace(tzinfo=datetime.timezone.utc)
    self.end   = end.replace(tzinfo=datetime.timezone.utc)

  class Channel(DR.Observation.Channel):
    """
    A GAVRT DSS-28 receiver channel if the signal for one freq, pol, and IFtype
    """
    def __init__(self, parent, name, freq=None, bw=None, pol=None):
      """
      """
      DR.Observation.Channel.__init__(self, parent, name,
                                      freq=freq, bw=bw, pol=pol)
        
    def load_tfile_data(self, path, start=None, stop=None, dss=28):
      """
      Load data from a GAVRT t-file

      @param filename : name of data file (t.YYDDD.C)
      @type  filename : str

      @param start : optional start time, default: beginning of file
      @type  start : datetime

      @param stop : optional stop time, default: end of file
      @type  stop : datetime

      @return: (indexed numpy data array,
                column labels for data array,
                list of date numbers, one for each row,
                list of right ascensions (hrs), one for each row,
                list of declinations (degs),
                list of Tsys counts (kHz))
      """
      filename = path+"."+str(int(self.name))
      datafile = open(filename,"r")
      labels = datafile.readline().strip().split()
      logger.debug("load_tfile_data: labels: %s", labels)
      datafile.close()
      labels.insert(0,'DOY')
      labels.insert(0,'Year')
      logger.debug("load_tfile_data: new labels: %s", labels)
    
      # colums are: Year DOY UTC Epoch Chan Tsys Int Az El Diode Level CryoTemp
      #              i4   i4  S8   f8    S2  f4   f4 f4 f4   i4    f4    f4
      data = numpy.loadtxt(filename,skiprows=1,
                     dtype = {'names': tuple(labels),
                              'formats': ('i4','i4','S8','f8',
                                     'S2','f4','f4','f4','f4','i4','f4', 'f4')})
      date_nums = []
      ras = []
      decs = []
      azs = []
      elevs = []
      tsys = []
      for index in range(len(data)):
        time = data[index]['Epoch']
        dt = datetime.datetime.utcfromtimestamp(time).replace(
                                                   tzinfo=datetime.timezone.utc)
        if ((start and stop) and (start <= dt and dt <= stop)) \
            or start == None or stop == None:
          time_tuple = (dt.year,
                        DT.day_of_year(dt.year,dt.month,dt.day)
                        + (  dt.hour
                           + dt.minute/60.
                           + dt.second/3600.
                           + dt.microsecond/3600./1e6)/24.)
          azimuth = data[index]['Az']
          elevation = data[index]['El']
          ra,dec = Astronomy.AzEl_to_RaDec(azimuth,
                                           elevation,
                                           self.parent.latitude,
                                           self.parent.longitude,
                                           time_tuple)
          date_nums.append(matplotlib.dates.date2num(dt))
          ras.append(ra)
          decs.append(dec)
          azs.append(float(azimuth))     # because openpyxl doesn't like
          elevs.append(float(elevation)) # <type 'numpy.float32'>
          tsys.append(data[index]['Tsys'])
      return data, labels, date_nums, ras, decs, azs, elevs, tsys
   
def zenith_gain(freq):
  """
  Antenna response vs frequency for DSS-28

  @param freq : frequency in MHz
  @type  freq : float

  @return: float
  """
  parfile = open(project_path
                 + "DSS-28_technical/efficiency_vs_freq_pars.pkl","r")
  pars = cPickle.load(parfile)
  parfile.close()
  effic = {}
  avg_effic = 0
  for key in list(pars.keys()):
    effic[key] = pars[key](freq)/100.
    avg_effic += effic[key]
  # right now I don't know what Pol A and Pol B are
  avg_effic /= len(list(pars.keys()))
  return avg_effic

def DSS28_beamwidth(freq):
  """
  Estimate DSS-28 beamwidth until we have something better.

  Consistent with {14: 0.045, 8.45: 0.064, 6: 0.0873, 3.1: 0.1725}
  @param freq : GHz
  @type  freq : float

  @return: beamwidth in deg (float)
  """
  return 0.54/freq

def get_session_t_files(project_path, date_info):
  """
  """
  session_path = project_path+"Observations/dss28/%4d/%03d/" % (date_info[1],date_info[4])
  logger.debug("get_session_t_files: session path is %s", session_path)
  wb_name = "%4d-%03d.xlsx" % (date_info[1],date_info[4])
  logger.debug("get_session_t_files: workbook name is %s", wb_name)
  wb_file = session_path+"/"+wb_name
  logger.debug("get_session_t_files: filename is %s", wb_file)
  f1 = glob.glob(session_path+"/t12*.?") # for chans 2, 4, 6, 8
  f1.sort()
  f2 = glob.glob(session_path+"/t12*.??") # for chans 10, 12, 14, 16
  f2.sort()
  files = f1+f2
  logger.debug("get_session_t_files: found files %s", files)
  return files

def extract_data(datatype, wb, start, stop, meta_column, files):
  """
  Extract map from the data files and create/fill a metadata worksheet

  First the first t-file is loaded to get the start and stop time.  These are
  used to create the meta-data sheet name.  If the sheet exists it is loaded,
  else it is created.  Creation begins with the copying the session meta-data
  sheet.  The 'create_metadata_sheet' adds additional empty columns.  Then for
  each row the appropriate t-file is read and the meta-data are obtained.

  @param wb : observation spreadsheet
  @type  wb : workbook instance

  @param start : datetime of first data point
  @type  start : mpl datenum

  @param stop : datetime of last data point
  @type  stop : mpl datenum

  @return: metadata worksheet instance
  """
  # Get metadata from first data file
  data, labels, date_nums, ras, decs, azs, elevs, tsys = \
    load_tfile_data(files[0])
  # get the index for the start time and the end time
  start_hr  = num2date(start).hour
  start_min = num2date(start).minute
  stop_hr   = num2date(stop).hour
  stop_min  = num2date(stop).minute
  startstr = "%02d%02d" % (start_hr,start_min)
  stopstr  = "%02d%02d" % (stop_hr, stop_min)
  if datatype == None:
    sheetname = "Other-"
  elif datatype.lower() == 'map':
    sheetname = "Map-"
  elif datatype.lower() == 'time-series':
    sheetname = "Plot-"
  elif datatype == 'boresight':
    sheetname = "Bore-"
  sheetname += startstr+"-"+stopstr
  logger.debug("Looking for sheet %s", sheetname)
  metadata_sheet = wb.get_sheet_by_name(sheetname)
  if metadata_sheet == None:
    # Doesn't exist, make it.
    logger.debug("Attempting to create worksheet %s", sheetname)
    try:
      metasheet = wb.get_sheet_by_name('Metadata')
    except Exception as details:
      logger.error("Could not get spreadsheet metadata", exc_info=True)
      return None
    try:
      wb.add_sheet(copy.deepcopy(metasheet))
    except Exception as details:
      logger.error("Could not copy metadata sheet", exc_info=True)
      return None
    else:
      logger.debug("New sheets: %s", str(wb.get_sheet_names()))
      metadata_sheet = wb.get_sheet_by_name('Metadata')
      logger.debug("Active sheet: %s", metadata_sheet.title)
      metadata_sheet.title = sheetname
      logger.debug("Active sheet was renamed to %s", metadata_sheet.title)

    #logger.debug("Creating worksheet %s named %s",
    #               metadata_sheet,sheetname)
    #metadata_sheet = create_metadata_sheet(metadata_sheet,sheetname)
  else:
    logger.debug("%s already exists",metadata_sheet.title)
  for filename in files:
    bname = os.path.basename(filename)
    # find the row for this file or create it
    row = support.excel.get_row_number(metadata_sheet,meta_column['File'],bname)
    logger.debug("%s meta data will go into row %s", bname,row)
    if row == None:
      row = metadata_sheet.get_highest_row()
      logger.debug("%s meta data are not yet in %s", bname,sheetname)
      logger.debug("%s meta data will now go into row %d",bname,row)
    # get data for this file
    data, labels, date_nums, ras, decs, azs, elevs, tsys = \
          load_tfile_data(filename)
    start_index = support.nearest_index(date_nums,start)
    stop_index = support.nearest_index(date_nums,stop)
    if start_index == -1 or stop_index == -1:
      wb.remove_sheet(metadata_sheet)
      return None
    else:
      metadata_sheet.cell(row=row,
                        column=meta_column['File']).value = bname
      metadata_sheet.cell(row=row,
                        column=meta_column["First"]).value = start_index
      metadata_sheet.cell(row=row,
                        column=meta_column["Last"] ).value = stop_index
      metadata_sheet.cell(row=row,
                        column=meta_column['Start']).value = num2date(start)
      metadata_sheet.cell(row=row,
                        column=meta_column['Stop'] ).value = num2date(stop)
  set_column_dimensions(metadata_sheet)
  return metadata_sheet
