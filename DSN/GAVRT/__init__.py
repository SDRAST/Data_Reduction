# -*- coding: utf-8 -*-
"""
A lot of modules under this need fixing up because of change from package 
Observatory to MonitorControl
"""
import datetime
from glob import glob
from os.path import basename
from pylab import *
# numpy has a function 'copy'.  Override it.
import copy
from matplotlib.pylab import date2num
import logging

import Astronomy
from support.excel import get_row_number
from Data_Reduction import nearest_index
from Data_Reduction.DSN.GAVRT.solar import create_metadata_sheet
from DSMS.excel import set_column_dimensions

mylogger = logging.getLogger("__main__."+__name__)

def load_tfile_data(filename, start=None, stop=None, dss=28):
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
  longitude,latitude,altitude,tz,name, diam = \
      Astronomy.get_geodetic_coords(dss=dss)
  datafile = open(filename,"r")
  labels = datafile.readline().strip().split()
  datafile.close()
  labels.insert(0,'DOY')
  labels.insert(0,'Year')

  # colums are: Year DOY UTC Epoch Chan Tsys Int Az El Diode Level
  #              i4   i4  S8   f8    S2  f4   f4 f4 f4   i4    f4
  data = loadtxt(filename,skiprows=1,
                 dtype = {'names': tuple(labels),
                          'formats': ('i4','i4','S8','f8',
                                      'S2','f4','f4','f4','f4','i4','f4')})
  date_nums = []
  ras = []
  decs = []
  azs = []
  elevs = []
  tsys = []
  for index in range(len(data)):
    time = data[index]['Epoch']
    dt = datetime.datetime.utcfromtimestamp(time)
    if ((start and stop) and (start <= dt and dt <= stop)) \
        or start == None or stop == None:
      time_tuple = (dt.year,
                    Astronomy.day_of_year(dt.year,dt.month,dt.day)
                    + (  dt.hour
                       + dt.minute/60.
                       + dt.second/3600.
                       + dt.microsecond/3600./1e6)/24.)
      azimuth = data[index]['Az']
      elevation = data[index]['El']
      ra,dec = Astronomy.AzEl_to_RaDec(azimuth,
                                       elevation,
                                       latitude,
                                       longitude,
                                       time_tuple)
      date_nums.append(date2num(dt))
      ras.append(ra)
      decs.append(dec)
      azs.append(float(azimuth))     # because openpyxl doesn't like
      elevs.append(float(elevation)) # <type 'numpy.float32'>
      tsys.append(data[index]['Tsys'])
  return data, labels, date_nums, ras, decs, azs, elevs, tsys

def plot_tsys(fig, date_nums, tsys, label=None, picker=3):
  """
  Plot system temperatures

  @param fig : figure which contains current axes
  @type  fig : figure() instance

  @param date_nums : datetime for each point
  @type  date_nums : list of mpl date numbers

  @param tsys : system temperature counts
  @type  tsys : list of float
  """
  #plot_date(date_nums, tsys, '-')
  plot_date(date_nums, tsys, '-', picker=picker)
  grid()
  text(1.05,0.5,label,
       horizontalalignment='left',
       verticalalignment='center',
       transform = gca().transAxes)
  title = num2date(date_nums[0]).strftime("%Y %j")
  fig.autofmt_xdate()
  return title

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
  for key in pars.keys():
    effic[key] = pars[key](freq)/100.
    avg_effic += effic[key]
  # right now I don't know what Pol A and Pol B are
  avg_effic /= len(pars.keys())
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
  mylogger.debug("get_session_t_files: session path is %s", session_path)
  wb_name = "%4d-%03d.xlsx" % (date_info[1],date_info[4])
  mylogger.debug("get_session_t_files: workbook name is %s", wb_name)
  wb_file = session_path+"/"+wb_name
  mylogger.debug("get_session_t_files: filename is %s", wb_file)
  f1 = glob(session_path+"/t12*.?") # for chans 2, 4, 6, 8
  f1.sort()
  f2 = glob(session_path+"/t12*.??") # for chans 10, 12, 14, 16
  f2.sort()
  files = f1+f2
  mylogger.debug("get_session_t_files: found files %s", files)
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
  mylogger.debug("Looking for sheet %s", sheetname)
  metadata_sheet = wb.get_sheet_by_name(sheetname)
  if metadata_sheet == None:
    # Doesn't exist, make it.
    mylogger.debug("Attempting to create worksheet %s", sheetname)
    try:
      metasheet = wb.get_sheet_by_name('Metadata')
    except Exception, details:
      mylogger.error("Could not get spreadsheet metadata", exc_info=True)
      return None
    try:
      wb.add_sheet(copy.deepcopy(metasheet))
    except Exception, details:
      mylogger.error("Could not copy metadata sheet", exc_info=True)
      return None
    else:
      mylogger.debug("New sheets: %s", str(wb.get_sheet_names()))
      metadata_sheet = wb.get_sheet_by_name('Metadata')
      mylogger.debug("Active sheet: %s", metadata_sheet.title)
      metadata_sheet.title = sheetname
      mylogger.debug("Active sheet was renamed to %s", metadata_sheet.title)

    #mylogger.debug("Creating worksheet %s named %s",
    #               metadata_sheet,sheetname)
    #metadata_sheet = create_metadata_sheet(metadata_sheet,sheetname)
  else:
    mylogger.debug("%s already exists",metadata_sheet.title)
  for filename in files:
    bname = basename(filename)
    # find the row for this file or create it
    row = get_row_number(metadata_sheet,meta_column['File'],bname)
    mylogger.debug("%s meta data will go into row %s", bname,row)
    if row == None:
      row = metadata_sheet.get_highest_row()
      mylogger.debug("%s meta data are not yet in %s", bname,sheetname)
      mylogger.debug("%s meta data will now go into row %d",bname,row)
    # get data for this file
    data, labels, date_nums, ras, decs, azs, elevs, tsys = \
          load_tfile_data(filename)
    start_index = nearest_index(date_nums,start)
    stop_index = nearest_index(date_nums,stop)
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
