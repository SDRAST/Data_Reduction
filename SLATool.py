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

from Data_Reduction import get_obs_dirs, get_obs_session
from Data_Reduction.FITS.DSNFITS import FITSfile
from Data_Reduction.FITS.SDFITSplotter import DSNFITSplotter
from DSMS import DSN_complex_of
from MonitorControl import Observatory, Telescope
from Radio_Astronomy import rms_noise

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
      for table_key in examiner[dfindex].tables.keys():
        examiner[dfindex].tables[table_key].report_table()
      dfindex += 1
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
                                      
  def fit_mean_power_to_airmass(self, Tvac_func, replace=False):
    """
    fits the mean power data vs airmass to the radiative transfer equation
    
    This assumes that every IF has a way of measuring power. The measured power
    is a single value along the first axis of the data array (or last index in
    a C/Python array).  If there are multiple records then they will be
    averaged.
    
    @param Tvac_func - a function for system temperature with no atmosph or CBR
    @type  Tvac_func - function(beam,pol)
    """
    def opacity_fitting(x, a, tau):
      x_rad = numpy.deg2rad(x)
      x_sec = 1/numpy.sin(x_rad)
      return a + tau*x_sec
    
    for beam_idx in range(self.props['num beams']):
      for IFidx in range(self.props['num IFs']):
        self.logger.debug('fit_mean_power_to_airmass: processing IF%d', IFidx+1)
        Tvac = Tvac_func(beam=beam_idx ,pol=IFidx)
        pol = ['L','R'][IFidx]
        msg = "estimated %sCP zenith Tsys in vacuum= %6.2f" % (pol, Tvac)
        for ex in self.examiners.keys():
          self.header.add_history(msg)
          # average records here if needed.
          subchannels = len(self.props['num cycles'])
          self.logger.debug("fit_mean_power_to_airmass: subchannels: %s",
                          subchannels)
          for subchannel in subchannels:
            subch = subchannels.index(subchannel)
            # Get the data for this subchannel.
            #   the following can be expressed as
            #     mean_power = self.data['avgpower'][subch::2,IFidx,0,0,0]
            #   or
            #     mean_power = self.data['avgpower'][:,IFidx,0,0,0][subch::2]
            #   or
            #     mean_power = self.data[subch::2]['avgpower'][:,IFidx,0,0,0]
            # assuming the first row is the same as the subchannel
            mean_power = self.data['avgpower'][subch::2,IFidx,0,0,0]
            elevation  = self.data['ELEVATIO'][subch::2]
            # get an elevation array with the 'nan' and zero values removed
            mask = ~(numpy.isnan(elevation) | numpy.equal(elevation, 0))
            # these are the indices in the data for this subchannel
            #     'where' and its friends return a tuple
            indices = numpy.where(mask)[0]
            self.logger.debug(
             "fit_mean_power_to_airmass: good data rows for subchannel %d: %s",
             subchannel, indices)
            elv = elevation[mask]
            # remove the items with the same indices from mean_power
            pwr = mean_power[mask]
            #self.logger.debug('fit_mean_power_to_airmass: elevations: %s', elv)
            self.logger.debug('fit_mean_power_to_airmass: subch %d pwr shape: %s',
                          subchannel, pwr.shape)
            # fit the data
            popt, pcov = curve_fit(opacity_fitting, elv, pwr, p0=[0, 0])
            intercept, slope = popt[0], popt[1]
            self.logger.debug(
                         "fit_mean_power_to_airmass: intercept, slope: %f, %f",
                         intercept, slope)
            msg = \
            "IF%d, subch%d gain=%9.3e counts, gain_slope=%9.3e counts/airmass"\
                % (IFidx+1, subchannel, intercept, slope)
            self.header.add_history(msg)
            if replace:
              gain = Tvac/intercept
              K_per_am = gain*slope
              self.logger.debug(
           "fit_mean_power_to_airmass: convert power to Tsys for subch%s %sCP",
                          subch+1, pol)
              new_indices = numpy.where(self.data['CYCLE'] == subchannel)[0]
              self.logger.debug(
                          "fit_mean_power_to_airmass: table rows for Tsys: %s",
                          new_indices)
              self.logger.debug(
                          "fit_mean_power_to_airmass: destination shape is %s",
                          self.data['TSYS'][new_indices,IFidx,0,0,0].shape)
              self.data['TSYS'][new_indices,IFidx,0,0,0] = gain * pwr
            else:
              self.logger.warning(
                        "fit_mean_power_to_airmass: failed; Tsys not computed")

