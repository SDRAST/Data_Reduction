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
import time

from matplotlib.font_manager import FontProperties
from pylab import *

from Data_Reduction import get_obs_dirs, get_obs_session
from Data_Reduction.FITS.DSNFITS import FITSfile
from Data_Reduction.FITS.SDFITSplotter import DSNFITSplotter,make_legend_labels
from DSMS import DSN_complex_of
from MonitorControl import Observatory, Telescope
from Radio_Astronomy import rms_noise

fontP = FontProperties()
fontP.set_size('x-small')

colors = ["b", "g", "r", "c", "m", "y"]

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
    self.logger = logging.getLogger(logger.name+".DSNFITSplotter")
    
    # get the session details
    if year and DOY:
      date = "%4d/%03d" % (year, DOY)
    else:
      date = None
    self.project, self.DSS, self.year, self.DOY = get_obs_session(dss=dss,
                                                    project=project, date=date)
    self.projectdatapath, self.projworkpath, self.datapath = \
      get_obs_dirs(self.project, self.DSS, self.year, self.DOY, datafmt="FITS")
    # get the datafiles to be processed     
    datafiles = glob.glob(self.datapath+"*.fits")
    datafiles.sort()
    self.logger.info("__init__: found DOY %03d %d datafiles",
                     self.DOY, len(datafiles))
    self.examiners = self.open_datafiles(datafiles)
    self.get_sources()
  
  def open_datafiles(self, datafiles):
    """
    Opens the files for the session in /usr/local/RA_data/FITS
    """
    dfindex = 0
    examiner = {}
    for datafile in datafiles:
      self.logger.info("open_datafiles: %s", os.path.basename(datafile))
      examiner[dfindex]  = DSNFITSplotter(parent=self, FITSfile=datafile)
      examiner[dfindex].datafile = datafile
      for table_key in examiner[dfindex].tables.keys():
        examiner[dfindex].tables[table_key].report_table()
      dfindex += 1
    self.examiner_keys = examiner.keys()
    self.examiner_keys.sort()
    self.logger.info("open_datafiles: started %d examiners" % 
                     len(examiner.keys()))
    return examiner

  def get_sources(self):
    """
    Gets a list of sources from all the datasets in the session
    """
    self.sources = {}
    keys = self.examiners.keys()
    for key in keys:
      source_names = self.examiners[key].get_sources()
      for name in source_names:
        if self.sources.has_key(name):
          self.sources[name].append(key)
        else:
          self.sources[name] = [key]
    return self.sources.keys()

  def get_good_weather_data(self):
    """
    finds usable weather data in all examiners
    
    True means that some table has good data but not necessarily all.
    """
    have_elev = False
    have_Tsys = False
    have_Tambient = False
    have_pressure = False
    have_humidity = False
    have_windspeed = False
    have_winddirec = False
    for exkey in self.examiner_keys:
      examiner = self.examiners[exkey]
      self.logger.debug("tsys_summary: data file is %s", examiner.datafile)
      # assume multiple tables are in time order
      tablekeys = examiner.tables.keys()
      tablekeys.sort
      for tablekey in tablekeys:
        table = examiner.tables[tablekey]
        good_data = table.get_good_rows()
        if good_data.has_key('elev'):
          have_elev = True
        if good_data.has_key('TSYS'):
          have_Tsys = True
        if good_data.has_key('Tambient'):
          have_Tambient = True
        if good_data.has_key('pressure'):
          have_pressure = True
        if good_data.has_key('humidity'):
          have_humidity = True
        if good_data.has_key('windspeed'):
          have_windspeed = True
        if good_data.has_key('winddirec'):
          have_winddirec = True      
    return good_data
  
  def get_average(self, source='67P_CG_201'):
    """
    Print r.m.s. noise for each dataset and all datasets together
    
    @param source : source for which averaging isdone
    """
    rowfmt = "%03d %5s   %d  %5.1f  %6.1f  %6.4f  %6.4f   %4.2f"
    first_spectrum = {0:True, 1:True}
    sum_y = {0:0, 1:0}
    sum_Tsys = {0:0, 0:1}
    sum_intgr = {0:0, 1:0}
    for exkey in self.examiners.keys():
      ex = self.examiners[exkey]
      print "FITS file", os.path.basename(ex.file)
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
    
  def plot_elev_and_Tsys(self, figtitle, good_wx_data, savepath):
    """
    """
    fig1, ax = subplots(nrows=1, ncols=2, squeeze=True)
    ax2 = ax[0].twinx()
    fig1.suptitle(figtitle)
    fig1.set_size_inches(12,5, forward=True)
    
    # axes 0: elevation and Tsys (or avg pwr) time
    if good_wx_data.has_key('elev') or good_wx_data.has_key('TSYS'):
      ax[0].xaxis.set_major_formatter( DateFormatter('%H:%M') )
    num_ex = len(self.examiner_keys)
    for key in self.examiner_keys:
      table_keys = self.examiners[key].tables.keys()
      num_tb = len(table_keys)
      table_keys.sort()
      for tkey in table_keys:
        table = self.examiners[key].tables[tkey]
        num_cy = table.props['num cycles']
        num_bm = table.props['num beams']
        num_pl = table.props['num IFs']
        good_wx_data = table.get_good_rows()
        if good_wx_data.has_key('elev'):
          # plot elevation vs time
          ax[0].plot_date(good_wx_data['mpltime'],
                        good_wx_data['elev'],"-k",label="elevation")
        # axes 0, right side, plot system temperature or average power vs time
        if good_wx_data.has_key('TSYS'):
          for subch_idx in range(num_cy):  # subchannels first
            for beam in range(num_bm): # beams second
              for pol in range(num_pl): #pols third
                label = make_legend_labels(dskeys=self.examiner_keys,
                                           tbkeys=table_keys,
                                           sckeys=range(num_cy),
                                           bmkeys=range(num_bm),
                                           plkeys=range(num_pl),
                                           dskey=key,
                                           tbkey=tkey, 
                                           sckey=subch_idx,
                                           bmkey=beam,
                                           plkey=pol)
                color_index = (subch_idx*num_bm + beam)*num_pl + pol
                #color_index = num_bm*(num_bm*subch_idx + beam) + pol
                ax2.plot_date(good_wx_data['mpltime'][subch_idx::table.props['num cycles']],
                        good_wx_data['TSYS'][subch_idx][beam][pol], marker='.',
                        color=colors[color_index], label=label)
        # right axes: plot tsys or average power vs airmass
        if good_wx_data.has_key('elev') and good_wx_data.has_key('TSYS'):
          for subch in range(num_cy):
            for beam in range(num_bm):
              for pol in range(num_pl):
                label = make_legend_labels(dskeys=self.examiner_keys,
                                           tbkeys=table_keys,
                                           sckeys=range(num_cy),
                                           bmkeys=range(num_bm),
                                           plkeys=range(num_pl),
                                           dskey=key,
                                           tbkey=tkey,
                                           sckey=subch,
                                           bmkey=beam,
                                           plkey=pol)
                color_index = (subch*num_bm + beam)*num_pl + pol
                #color_index = table.props['num beams']*(table.props['num cycles']*subch_idx + beam) + pol
                ax[1].plot(1/sin(pi*array(good_wx_data['elev'])/180.)[subch::table.props['num cycles']],
                       good_wx_data['TSYS'][subch][beam][pol],
                       color=colors[color_index],
                       label=label)
    if good_wx_data.has_key('elev') == False:
      ax[0].text(0.5, 0.6,'bad elevation data',
               horizontalalignment='center',
               verticalalignment='center',
               transform = ax[0].transAxes)
    ax[0].set_xlabel("UT")
    ax[0].set_ylabel("Elevation (deg)")
    ax[0].grid()
    if good_wx_data.has_key('TSYS'):
      ax2.set_ylabel("T$_{sys}$ (K)")
    else:
      ax2.set_ylabel("average power")
    ax2.grid()
    fig1.autofmt_xdate()
    if good_wx_data.has_key('elev') == False:
      ax[1].text(0.5, 0.5,'bad elevation data',
               horizontalalignment='center',
               verticalalignment='center',
               transform = ax[1].transAxes)         
    ax[1].set_xlabel("Airmass")
    ax[1].yaxis.set_major_formatter(NullFormatter())
    ax[1].grid()
    lines, labels = ax[1].get_legend_handles_labels()
    fig1.legend(lines, labels, loc="upper right", ncol=1, prop = fontP)
    fig1.show()
    fig1.savefig(savepath+"-elev.png")

  def plot_weather(self, good_wx_data, savepath):
    """
    """
    fig2, wax = subplots(nrows=3, ncols=1, squeeze=True)
    fig2.suptitle("Weather")
    fig2.set_size_inches(6,3, forward=True)
    fig2.subplots_adjust(hspace=0) # no space between plots in a column
    fig2.subplots_adjust(left=0.15)
    for key in self.examiner_keys:
      table_keys = self.examiners[key].tables.keys()
      table_keys.sort()
      for tkey in table_keys:
        table = self.examiners[key].tables[tkey]
        good_wx_data = table.get_good_rows()
        # left axes: temperature
        if good_wx_data.has_key('Tambient'):
          wax[0].plot_date(good_wx_data['mpltime'], good_wx_data['Tambient'],"-k")
        # middle axes: pressure
        if good_wx_data.has_key('pressure'):
          wax[1].plot_date(good_wx_data['mpltime'], good_wx_data['pressure'],"-k")
        # right axes: humidity
        if good_wx_data.has_key('humidity'):
          wax[2].plot_date(good_wx_data['mpltime'], good_wx_data['humidity'],"-k")
    if good_wx_data.has_key('Tambient') or good_wx_data.has_key('pressure'):
      wax[2].xaxis.set_major_formatter( DateFormatter('%H:%M') )
    if good_wx_data.has_key('Tambient') == False:
      wax[0].text(0.5,0.5,'bad ambient temperature data',
                horizontalalignment='center',
                verticalalignment='center',
                transform = wax[0].transAxes)
    wax[0].set_ylabel("Temp (C)")
    wax[0].grid(True)
    if good_wx_data.has_key('pressure'):
      wax[1].yaxis.set_major_formatter( FormatStrFormatter('%6.2f') )
    else:
      wax[1].text(0.5,0.5,'bad pressure data',
                horizontalalignment='center',
                verticalalignment='center',
                transform = wax[1].transAxes)
    wax[1].grid(True)
    wax[1].set_ylabel("Pres (mb)")
    if good_wx_data.has_key('pressure'):
      wax[1].yaxis.set_major_formatter( FormatStrFormatter('%6.2f') )
    if good_wx_data.has_key('humidity') == False:
      wax[2].text(0.5,0.5,'bad humidity data',
                horizontalalignment='center',
                verticalalignment='center',
                transform = wax[2].transAxes)
    wax[2].grid(True)
    wax[2].set_ylabel("Humidity")
    if good_wx_data.has_key('Tambient') or good_wx_data.has_key('pressure') or \
      good_wx_data.has_key('humidity'):
      fig2.autofmt_xdate()
    fig2.show()
    fig2.savefig(savepath+"-weather.png")

  def plot_wind(self, good_wx_data, savepath):
    """
    """
    fig3, wax3 = subplots(nrows=2, ncols=1, squeeze=True)
    fig3.suptitle("Wind")
    fig3.set_size_inches(6,4.5, forward=True)
    fig3.subplots_adjust(hspace=0) # no space between plots in a column
    for key in self.examiner_keys:
      table_keys = self.examiners[key].tables.keys()
      table_keys.sort()
      for tkey in table_keys:
        table = self.examiners[key].tables[tkey]
        good_wx_data = table.get_good_rows()
        # left axes: windspeed
        if good_wx_data.has_key('windspeed'):
          wax3[0].plot_date(good_wx_data['mpltime'],good_wx_data['windspeed'],"-k")
        # right axes: wind direction
        if good_wx_data.has_key('winddirec'):
          wax3[1].plot_date(good_wx_data['mpltime'], good_wx_data['winddirec'],"-k")
    if good_wx_data.has_key('windspeed') or good_wx_data.has_key('winddirec'):
      wax3[1].xaxis.set_major_formatter( DateFormatter('%H:%M') )
    if good_wx_data.has_key('windspeed') == False:
      wax3[0].text(0.5,0.5,'bad wind speed data',
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = wax3[0].transAxes)
    wax3[0].set_ylabel("Speed (km/h)")
    wax3[0].grid(True)
    if good_wx_data.has_key('winddirec') == False:
      wax3[1].text(0.5,0.5,'bad wind direction data',
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = wax3[1].transAxes)
    wax3[1].set_ylabel("Direction")
    wax3[1].grid(True)
    if good_wx_data.has_key('windspeed') or good_wx_data.has_key('winddirec'):
      fig3.autofmt_xdate()
    fig3.show()
    fig3.savefig(savepath+"-wind.png")
      
  def consolidate(self, tables=[], sources=[], bands=['22U']):
    """
    Merge selected data into a new SDFITS file
    """
    if tables:
      exkeys = tables.sort()
    else:
      exkeys = self.examiners.keys()
    exkeys.sort()
    numrows = 0 # number of rows in the new table
    numrecs = 1 # maximum number of records along the time axis
    # start a new file with one table and verify consistency
    for exkey in exkeys:
      ex_idx = exkeys.index(exkey)
      if ex_idx:
        pass
      else:
        # first examiner (ex_idx=0); start a FITS file
        station = self.examiners[exkey].header['TELESCOP'].strip()
        complx = DSN_complex_of(station)
        obs = Observatory(complex)
        dss = self.examiners[exkey].tables[0].dss
        tel = Telescope(obs, dss=dss)
        newfile = FITSfile(tel)
        tb = self.examiners[exkeys[ex_idx]].tables[0]
        num_chans = tb.props['num chans']
        num_IFs = tb.props['num IFs']
        time_axis = tb.props['time axis']
        num_beams = tb.props['num beams']
      # process all the original FITS files, check for consistency
      for tbkey in self.examiners[exkey].tables.keys():
        tb = self.examiners[exkey].tables[tbkey]
        numrows += len(tb.data)
        for key in tb.props['num records'].keys():
          numrecs = max(numrecs, tb.props['num records'][key])
        if tb.props['num chans'] != num_chans or \
           tb.props['num IFs']   != num_IFs   or \
           tb.props['time axis'] != time_axis or \
           tb.props['num beams'] != num_beams:
          self.logger.error("consolidate: SPECTRUM mismatch at ex %d", exkey)
    newfile.header = newfile.make_basic_header() # table header
    newfile.make_basic_columns() # includes make_offset_columns()
    newfile.add_time_dependent_columns(numrecs)
    # make the data column
    axis = 1; newfile.make_data_axis(newfile.header, newfile.columns, axis,
                                     num_chans, 'FREQ-OBS', 'D', unit='Hz',
                                comment="channel frequency in telescope frame")
    fmt_multiplier = newfile.header['MAXIS1']
    dimsval = "("+str(newfile.header['MAXIS1'])+","
    axis +=1; newfile.make_data_axis(newfile.header, newfile.columns, axis,
                                     1, 'RA---GLS', 'D', unit='deg',
                                     comment="RA J2000.0")
    fmt_multiplier *= newfile.header['MAXIS2']
    dimsval += str(newfile.header['MAXIS2'])+","
    axis +=1; newfile.make_data_axis(newfile.header, newfile.columns, axis,
                                     1, 'DEC--GLS','D', unit='deg',
                                     comment="decl. J2000")
    fmt_multiplier *= newfile.header['MAXIS3']
    dimsval += str(newfile.header['MAXIS3'])+","
    # this ends the required columns
    if num_IFs > 1 or time_axis:
      axis+=1; newfile.make_data_axis(newfile.header, newfile.columns, axis,
                                   num_IFs, 'STOKES',  'I',
                                   comment="polarization code: -8 ... 4")
      fmt_multiplier *= newfile.header['MAXIS4']
      dimsval += str(newfile.header['MAXIS4'])+"," 
    if time_axis:
      axis+=1; newfile.make_data_axis(newfile.header, newfile.columns, axis,
                                   numrecs, 'TIME', 'E', unit='s',
                                   comment='Python time.time() value')
      fmt_multiplier *= newfile.header['MAXIS5']
      dimsval += str(newfile.header['MAXIS5'])+","
    if num_beams > 1:
      axis+=1; newfile.make_data_axis(newfile.header, newfile.columns, axis,
                                   num_beams, 'BEAM',    'I',
                                                comment="beam 1 or 2")
      fmt_multiplier *= newfile.header['MAXIS6']
      dimsval += str(newfile.header['MAXIS6'])+","
    # this ends the optional columns; now make the data column
    dimsval = dimsval[:-1]+")"
    data_format = str(fmt_multiplier)+"E"
    newfile.columns.add_col(pyfits.Column(name='SPECTRUM', format=data_format,
                         dim=dimsval))
    # copy missing columns and header entries
    for exkey in exkeys:
      for tbkey in self.examiners[exkey].tables.keys():
        tb = self.examiners[exkey].tables[tbkey]
        # check for missing columns
        for name in tb.columns.names:
          if name in newfile.columns.names:
            pass
          else:
            self.logger.debug("consolidate: ex%d tb%d has %s",
                              exkey, tbkey, name)
            name_idx = tb.data.columns.names.index(name)
            colfmt = tb.data.columns.formats[name_idx]
            self.logger.debug("consolidate: column format is %s", colfmt)
            tdimkey = 'TDIM'+str(name_idx+1)
            if tdimkey in tb.header.keys():
              # this column has multidimensional data
              self.logger.debug("consolidate: processing %s", tdimkey)
              coldim = tb.header[tdimkey]
              # the column dimension 'coldim' is a string in FITS order
              self.logger.debug("consolidate: got DIM %s (type %s)",
                                coldim, type(coldim))
              dimparts = coldim.split(',')
              timesize = int(dimparts[4])
              dimparts[4] = str(numrecs)
              newdim = ','.join(dimparts)
              self.logger.debug("consolidate: new DIM %s (type %s)",
                                newdim, type(newdim))
              # resize the format to be compatible with new number or records
              fmtsize = int(colfmt[:-1])
              fmttype = colfmt[-1]
              recsize = fmtsize/timesize
              newfmt = str(numrecs*recsize)+fmttype
              newfile.columns.add_col(pyfits.Column(name=name, format=newfmt,
                                                    dim=newdim))
            else:
              self.logger.debug("consolidate: adding column %s", name)
              newfile.columns.add_col(pyfits.Column(name=name, format=colfmt))
        # copy missing header entries
        stripped_header = tb.header.copy(strip=True)
        for keyword in stripped_header:
          if keyword in newfile.header.keys():
            # maybe check if same?
            pass
          else:
            newfile.header.append(tb.header.cards[keyword])
            self.logger.debug('consolidate: copied header %s', keyword)
    # now the new table has all the columns
    # create the BINTABLE HDU (only one)
    FITSrec = pyfits.FITS_rec.from_columns(newfile.columns, nrows=numrows)
    newfile.tables[0] = pyfits.BinTableHDU(data=FITSrec, header=newfile.header, name="SINGLE DISH")
    # copy the rows
    newrow=0
    for exkey in exkeys:
      for tbkey in self.examiners[exkey].tables.keys():
        tb = self.examiners[exkey].tables[tbkey]
        for row in range(len(tb.data)):
          self.logger.debug("consolidate: copying ex%d tb%d row %d to %d",
                            exkey, tbkey, row, newrow)
          # pad the SPECTRUM array along the TIME axis
          inrow = tb.data[row]
          outrow = newfile.tables[0].data[newrow] # copy empty row to'outrow'
          # copy each column into 'outrow'
          for name in tb.data.columns.names:
            self.logger.debug("consolidate: copying %s", name)
            if type(inrow[name]) == numpy.ndarray:
              insize = inrow[name].shape[-5] # time is the fifth axis
              outsize = newfile.tables[0].data[newrow][name].shape[-5]
              self.logger.debug("consolidate: %d records in, %d records out",
                                insize, outsize)
              pad = outsize - insize
              outarray = numpy.pad(inrow[name], 
                                   ((0,0),(0,pad),(0,0),(0,0),(0,0),(0,0)),
                                   'constant', constant_values=0)
              # update the cell
              newfile.tables[0].data[newrow][name] = outarray
            else:
              newfile.tables[0].data[newrow][name] = inrow[name]
          newrow += 1
    hdulist = pyfits.HDUList([newfile.prihdu] + newfile.tables.values())
    t = time.gmtime()
    fname_root = self.projworkpath + "mrg" + time.strftime("_%Y-%j-%H%M%S", t)
    savename = fname_root+".fits"
    self.logger.info('consolidate: writing %s', savename)
    try:
      hdulist.writeto(savename, clobber=False)
      saved = True
    except:
      self.logger.error('consolidate: failed; try new name')
      saved = False
    new_exkey = exkeys[-1]+1
    self.examiners[new_exkey]  = DSNFITSplotter(parent=self, hdulist=hdulist)
    if saved:
      self.examiners[new_exkey].file = savename
    self.logger.info("consolidate: created new examiner/plotter %d", new_exkey)
    return self.examiners[new_exkey]

