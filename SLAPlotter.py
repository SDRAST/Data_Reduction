"""
Session plotting methods

Plotting component of an environment for analyzing spectral line from multiple
observatories. Manages various data reduction programs.
"""
import logging

from matplotlib.dates import DateFormatter
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FormatStrFormatter
from os.path import basename, dirname, splitext
from pylab import *

from Data_Reduction.FITS.SDFITSplotter import make_legend_labels, DSNFITSplotter
from Data_Reduction.SLATool import SessionAnalyzer
from Data_Reduction.tipping import airmass

logger = logging.getLogger(__name__)

fontP = FontProperties()
fontP.set_size('x-small')

brown = '#630000'
teal = '#0084a5'
colors = ["b", "g", "r", "c", "m", "y", brown, teal]
sigref = {False: "r", True:"s"}

class SessionPlotter(SessionAnalyzer):
  """
  """
  def __init__(self, project=None, dss=None, year=None, DOY=None):
    """
    """
    SessionAnalyzer.__init__(self, project=project, dss=dss,
                                   year=year, DOY=DOY)
    self.logger = logging.getLogger(logger.name+".SessionAnalyzer")
    self.logger.info("__init__: getting FITS examiners for DSS%2d on %4d/%03d",
                     dss, year, DOY)
    if self.examiners:
      for ekey in list(self.examiners.keys()):
        self.examiners[ekey].plotter = {}
        for tkey in list(self.examiners[ekey].tables.keys()):
          self.examiners[ekey].plotter[tkey] = DSNFITSplotter.Plotter(self,
                                             self.examiners[ekey].tables[tkey])
    else:
      self.logger.warning("__init__: no FITS examiners created")
    
  def plot_elev_and_Tsys(self, figtitle=None, weather_data=None,
                         examiner_keys=None, savepath=None):
    """
    Plots Tsys vs time and vs elevation.
    
    The data asociated with each key of 'weather_data' is a dict with numpy
    array for (SIG state) True and for False.    The 'TSYS' array has four
    axes representing::
      time index   - 0-based sequence in order of matplotlib datenum 
      subchannel   - CYCLE value
      beam         - 1-based number sequence
      IF           - 1-based number sequence, usually representing pol
    The other keys have only a time axis.
    
    @param figtitle : figure title
    @type  figtitle : str
    
    @param weather_data : consolidated environmental data
    @type  weather_data : dict
    
    @param examiner_keys : keys of files from this date to be included
    @type  examiner_keys : list of int
    
    @param savepath : path to directory to save figure
    @type  savepath : str
    """
    if weather_data:
      pass
      self.logger.debug("plot_elev_and_Tsys: using provided data")
    elif examiner_keys:
      self.logger.debug("plot_elev_and_Tsys: extracting data from tables %s",
                        examiner_keys)
      weather_data = self.get_good_weather_data(examiner_keys=examiner_keys)
    else:
      self.logger.debug("plot_elev_and_Tsys: using all data for this date")
      weather_data = self.get_good_weather_data()
      
    fig1, ax = subplots(nrows=1, ncols=2, squeeze=True)
    if figtitle:
      fig1.suptitle(figtitle)
    else:
      fig1.suptitle("%4d/%03d DSS-%2d" % (self.year, self.DOY, self.DSS))
    fig1.set_size_inches(12,5, forward=True)
    ax2 = ax[0].twinx()
    
    if 'ELEVATIO' in weather_data or 'TSYS' in weather_data:
      ax[0].xaxis.set_major_formatter( DateFormatter('%H:%M') )
    if 'TSYS' in weather_data:
      nrows, num_cy, num_bm, num_pl = weather_data['TSYS'][True].shape
    states = list(weather_data['TSYS'].keys())
    num_st = len(states)
    for sig in states:
      # axes 0: left side, elevation vs time
      if 'ELEVATIO' in weather_data:
        # plot elevation vs time
        mpltime = epoch2num(weather_data['UNIXtime'][True])
        ax[0].plot_date(mpltime, weather_data['ELEVATIO'][True], "-k",
                        label=sigref[sig])
      # axes 0, right side, plot system temperature or average power vs time
      if 'TSYS' in weather_data:
        for subch_idx in range(num_cy):  # subchannels first
          for beam in range(num_bm): # beams second
            for pol in range(num_pl): #pols third
              label = make_legend_labels(sckeys=list(range(num_cy)),
                                         bmkeys=list(range(num_bm)),
                                         plkeys=list(range(num_pl)),
                                         sckey=subch_idx,
                                         bmkey=beam,
                                         plkey=pol)
              color_index = num_st*((subch_idx*num_bm + beam)*num_pl + pol) \
                            + 1-int(sig)
              tsys_len = len(weather_data['TSYS'][sig][:,subch_idx,beam,pol])
              mpltime = epoch2num(weather_data['UNIXtime'][sig])[:tsys_len]
              ax2.plot_date(mpltime,
                            weather_data['TSYS'][sig][:,subch_idx,beam,pol],
                            marker='.',
                            color=colors[color_index]) # , label=label+"s")
      # right axes: plot tsys or average power vs airmass
      if 'ELEVATIO' in weather_data and 'TSYS' in weather_data:
        for subch in range(num_cy):
          for beam in range(num_bm):
            for pol in range(num_pl):
              label = make_legend_labels(sckeys=list(range(num_cy)),
                                         bmkeys=list(range(num_bm)),
                                         plkeys=list(range(num_pl)),
                                         sckey=subch,
                                         bmkey=beam,
                                         plkey=pol)
              color_index = num_st*((subch*num_bm + beam)*num_pl + pol) \
                            + 1-int(sig)
              tsys_len = len(weather_data['TSYS'][sig][:,subch_idx,beam,pol])
              am = airmass(weather_data['ELEVATIO'][sig])[:tsys_len]
              ax[1].plot(am,  weather_data['TSYS'][sig][:,subch,beam,pol],
                         '.', color=colors[color_index],
                         label=label+sigref[sig])
    # set axis labels or report missing data
    if ('ELEVATIO' in weather_data) == False:
      ax[0].text(0.5, 0.6,'bad elevation data',
               horizontalalignment='center',
               verticalalignment='center',
               transform = ax[0].transAxes)
    ax[0].set_xlabel("UT")
    ax[0].set_ylabel("Elevation (deg)")
    ax[0].grid()
    if 'TSYS' in weather_data:
      if ax2.get_ylim()[1] > 1000:
        ax2.set_ylabel("power (count)")
      elif ax2.get_ylim()[1] > 10:
        ax2.set_ylabel("T$_{sys}$ (K)")
      else:
        ax2.set_ylabel("power (count)")
    ax2.grid()
    for tick in ax[0].get_xticklabels():
      tick.set_rotation(30)

    if ('ELEVATIO' in weather_data) == False:
      ax[1].text(0.5, 0.5,'bad elevation data',
               horizontalalignment='center',
               verticalalignment='center',
               transform = ax[1].transAxes)         
    ax[1].set_xlabel("Airmass")
    ax[1].yaxis.set_major_formatter(NullFormatter())
    ax[1].grid()
    lines, labels = ax[1].get_legend_handles_labels()
    fig1.legend(lines, labels, numpoints=1, loc="upper right", ncol=1,
                prop = fontP)
    if savepath:
      fig1.savefig(savepath+"-elev.png")
      self.logger.info("plot_elev_and_Tsys: saved to %s", savepath+"-elev.png")
    else:
      fig1.savefig(self.datapath+"-elev.png")
      self.logger.info("plot_elev_and_Tsys: saved to %s",
                       self.datapath+"-elev.png")
    #fig1.show()

  def plot_weather(self, figtitle=None, weather_data=None, examiner_keys=None,
                   savepath=None):
    """
    Plots temperature, humidity and pressure from get_good_weather_data()
    """
    if weather_data:
      pass
      self.logger.debug("plot_elev_and_Tsys: using provided data")
    elif examiner_keys:
      self.logger.debug("plot_elev_and_Tsys: extracting data from tables %s",
                        examiner_keys)
      weather_data = self.get_good_weather_data(examiner_keys=examiner_keys)
    else:
      self.logger.debug("plot_elev_and_Tsys: using all data for this date")
      weather_data = self.get_good_weather_data()
    fig2, wax = subplots(nrows=3, ncols=1, squeeze=True)
    if figtitle:
      fig2.suptitle(figtitle)
    else:
      fig2.suptitle("Weather")
    fig2.set_size_inches(6,3, forward=True)
    fig2.subplots_adjust(hspace=0) # no space between plots in a column
    fig2.subplots_adjust(left=0.15)
    mpltime = epoch2num(weather_data['UNIXtime'][True])
    # left axes: temperature
    if 'TAMBIENT' in weather_data:
      wax[0].plot_date(mpltime, weather_data['TAMBIENT'][True],"-k")
    # middle axes: pressure
    if 'PRESSURE' in weather_data:
      wax[1].plot_date(mpltime, weather_data['PRESSURE'][True],"-k")
    # right axes: humidity
    if 'HUMIDITY' in weather_data:
      wax[2].plot_date(mpltime, weather_data['HUMIDITY'][True],"-k")
    if 'TAMBIENT' in weather_data or 'PRESSURE' in weather_data:
      wax[2].xaxis.set_major_formatter( DateFormatter('%H:%M') )
    if ('TAMBIENT' in weather_data) == False:
      wax[0].text(0.5,0.5,'bad ambient temperature data',
                horizontalalignment='center',
                verticalalignment='center',
                transform = wax[0].transAxes)
    wax[0].set_ylabel("Temp (C)")
    wax[0].grid(True)
    if 'PRESSURE' in weather_data:
      wax[1].yaxis.set_major_formatter( FormatStrFormatter('%6.2f') )
    else:
      wax[1].text(0.5,0.5,'bad pressure data',
                horizontalalignment='center',
                verticalalignment='center',
                transform = wax[1].transAxes)
    wax[1].grid(True)
    wax[1].set_ylabel("Pres (mb)")
    if 'PRESSURE' in weather_data:
      wax[1].yaxis.set_major_formatter( FormatStrFormatter('%6.2f') )
    if ('HUMIDITY' in weather_data) == False:
      wax[2].text(0.5,0.5,'bad humidity data',
                horizontalalignment='center',
                verticalalignment='center',
                transform = wax[2].transAxes)
    wax[2].grid(True)
    wax[2].set_ylabel("Humidity")
    if 'TAMBIENT' in weather_data or 'PRESSURE' in weather_data or \
      'HUMIDITY' in weather_data:
      fig2.autofmt_xdate()
    #fig2.show()
    fig2.savefig(savepath+"-weather.png")

  def plot_wind(self, figtitle=None, weather_data=None, examiner_keys=None,
                savepath=None):
    """
    Plots wind velocity and direction
    """
    if weather_data:
      pass
      self.logger.debug("plot_elev_and_Tsys: using provided data")
    elif examiner_keys:
      self.logger.debug("plot_elev_and_Tsys: extracting data from tables %s",
                        examiner_keys)
      weather_data = self.get_good_weather_data(examiner_keys=examiner_keys)
    else:
      self.logger.debug("plot_elev_and_Tsys: using all data for this date")
      weather_data = self.get_good_weather_data()
    fig3, wax3 = subplots(nrows=2, ncols=1, squeeze=True)
    if figtitle:
      fig3.suptitle(figtitle)
    else:
      fig3.suptitle("Wind")
    fig3.set_size_inches(6, 4.5, forward=True)
    fig3.subplots_adjust(hspace=0) # no space between plots in a column
    mpltime = epoch2num(weather_data['UNIXtime'][True])
    # left axes: windspeed
    if 'WINDSPEE' in weather_data:
      wax3[0].plot_date(mpltime, weather_data['WINDSPEE'][True],"-k")
      wax3[0].set_ylabel("Speed (km/h)")
      wax3[0].grid(True)
    else:
      wax3[0].text(0.5,0.5,'bad wind speed data',
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = wax3[0].transAxes)
    # right axes: wind direction
    if 'WINDDIRE' in weather_data:
      wax3[1].plot_date(mpltime, weather_data['WINDDIRE'][True],"-k")
      wax3[1].xaxis.set_major_formatter( DateFormatter('%H:%M') )
      wax3[1].set_ylabel("Direction")
      wax3[1].grid(True)
    else:
      wax3[1].text(0.5,0.5,'bad wind direction data',
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = wax3[1].transAxes)
    if 'WINDSPEE' in weather_data or 'WINDDIRE' in weather_data:
      fig3.autofmt_xdate()
    #fig3.show()
    fig3.savefig(savepath+"-wind.png")

  def plot_passband(self, figtitle=None):
    """
    Plots the passbands as dynamic spectra for the whole session.
    
    This plots the passbands for each SINGLE DISH table in the session.
    Generally, there is only one.
    """
    for dfindex in self.examiner_keys:
      session_name = splitext(basename(self.examiners[dfindex].file))[0]
      savefile = self.datapath + session_name   
      for tablekey in list(self.examiners[dfindex].tables.keys()):
        table = self.examiners[dfindex].tables[tablekey]
        plotter = self.examiners[dfindex].plotter[tablekey]
        if len(table.obsmodes) > 1:
          raise RuntimeError("multiple observing modes not yet supported")
        plotter.show_passband(
                      figtitle=session_name.replace("_"," ")+"-"+str(tablekey),
                             savepath=savefile+"-"+str(tablekey)+"_specgm.png")  
      
  def plot_bmsw_diff(self, figtitle=None, savefig=True):
    """
    Plot the difference between the two beams for each pol
    
    This is an indicator of receiver gain stability.
    """
    for dfindex in self.examiner_keys:
      examiner = self.examiners[dfindex]
      for tablekey in list(examiner.tables.keys()):
        table = self.examiners[dfindex].tables[tablekey]
        plotter = self.examiners[dfindex].plotter[tablekey]
        if len(table.obsmodes) > 1:
          raise RuntimeError("multiple observing modes not yet supported")
        num_scans = len(table.scan_keys)
        num_cycles = len(table.cycle_keys)
        num_rows = len(table.acs_rows)
        num_beams = table.props['num beams']
        num_pols = table.props["num IFs"]
        # collect the diagnostic spectra
        if table.props["full Stokes"]:
          # when SPECTRA are Stokes parameters, plot IF power for two pol modes
          if table.props["num IFs"] == 2:
            datasource = "IFSPECTR"
            num_chans = table.props['num IFspec chans']
        else:
          if "SPECTRUM" in table.data.columns.names:
            datasource = "SPECTRUM"
          else:
            datasource = "DATA"
          num_chans = table.props['num chans']
        self.logger.debug(" data source is %s", datasource)
      
        # create arrays of the right dimensions
        spectra = table.prepare_summary_arrays(num_chans)

        # get data statistics for scaling plots
        ymin, ymax, ymean, ystd = table.get_data_stats()
      
        # get the data
        cycles = table.cycle_keys
        for scan in table.scan_keys:
          scan_idx = table.scan_keys.index(scan) # scan numbers can start anywhere
          for cycle in cycles:
            subch_idx = cycle - 1
            for beam_idx in range(table.props["num beams"]):
              beam = beam_idx+1
              for IF_idx in range(table.props["num IFs"]):
                pol = IF_idx+1
                #self.logger.debug("plot_bmsw_diff: processing scan %d, subch %d, beam %d, pol %d",
                #             scan, cycle, beam, pol)
                if table.props["time axis"]:
                  # average scan spectrum and include record spectra in image
                  # assume there is a scan number equal to the cycle number
                  image, spectrum = \
                      table.average_records(scan, cycle, beam, pol)
                  # this is the average spectrum for the scan
                  spectra[scan_idx][subch_idx][beam_idx][IF_idx] = spectrum
                else:
                  # no time axis
                  spec_indices = table.get_indices(scan=scan, cycle=cycle, 
                                            beam=beam, pol=pol)
                  spectra[scan_idx][subch_idx][beam_idx][IF_idx] = \
                                           table.data[datasource][spec_indices]

        # number of summaries
        #     we want a summary for every beam, polarization and subchannel
        if table.props["num beams"] < 2:
          self.logger.error("plot_bmsw_diff: %d beams is not enough", 
                            table.props["num beams"])
          return False
        num_summar = table.props["num beams"]/2 \
                    *table.props["num IFs"] * num_cycles

        # One row for dynamic spectra of scans or records
        if figtitle:
          pass
        else:
          figtitle = basename(self.examiners[0].file)[4:-5]
        # for beam-1 minus beam-2 differences
        # THIS DOES NOT YET HANDLE OBSMODE CHANGES
        if table.data[0]['OBSMODE'] == 'LINEPBSW':
          fig, ax = plotter.init_multiplot(
                            figtitle+"-"+str(tablekey)
                            +" beam 1 - beam 2, summed over records",
                            nrows=1, ncols=num_summar)
          col = 0
          labels = {}
          for subch in range(num_cycles):
            # subchannel index, not subchannel number
            for beam in range(0,table.props["num beams"],2):
              # beams taken by pairs
              for pol in range(table.props["num IFs"]):
                label = make_legend_labels(sckeys=list(range(num_cycles)),
                                           plkeys=list(range(num_pols)),
                                           sckey=subch,
                                           plkey=pol)                           
                # beam-1 minus beam-2 differences
                ax[col].set_title(label)
                for scan in range(len(table.scan_keys)):
                  # scan index, not SCAN number
                  beamdiff = spectra[scan][subch][0][pol] - \
                             spectra[scan][subch][1][pol]
                  ax[col].plot(beamdiff, label=str(scan))
                ax[col].grid(True)
                ax[col].set_xlim(0,table.props['num chans'])
                for tick in ax[col].get_xticklabels():
                  tick.set_rotation(45)
                fig.subplots_adjust(top=0.88)
                fig.subplots_adjust(bottom=0.15)
                col += 1
          last_col = len(ax)-1
          lines, labels = ax[last_col].get_legend_handles_labels()
          fig.legend(lines, labels, loc="upper right", ncol=2, prop = fontP)
          datasetID = splitext(basename(examiner.file))[0]+"-"+str(tablekey)
          fig.savefig(self.datapath+datasetID+"_beam_diff.png")
          show()
        else:
          self.logger.warning("plot_bmsw_diff: examiner %d table %d mode is %s",
                              dfindex, tablekey, table.data[0]['OBSMODE'])
      # end table loop
    # end examiner loop

  def plot_possw_diff(self, figtitle=None, savefig=True):
    """
    Plots the difference between a SIG=True scan and the next SIG=False scan.
    
    This eliminates receiver systematics.
    """
    for dfindex in self.examiner_keys:
      session_name = splitext(basename(self.examiners[dfindex].file))[0]
      savefile = self.datapath + session_name
      self.logger.debug("plot_possw_diff: saving as %s", savefile)
      for tablekey in list(self.examiners[dfindex].tables.keys()):
        plotter = self.examiners[dfindex].plotter[tablekey]
        if len(plotter.obsmodes) > 1:
          raise RuntimeError("multiple observing modes not yet supported")
        obsmode = plotter.get_first_good_value('OBSMODE')
        if obsmode == 'LINEPSSW' or obsmode == 'LINEPBSW':
          if savefig:
            plotter.plot_PSSW_spectra(
                      figtitle=session_name.replace("_"," ")+"-"+str(tablekey),
                             savepath=savefile+"-"+str(tablekey)+"_on-off.png")
          else:
            plotter.plot_PSSW_spectra(
                      figtitle=session_name.replace("_"," ")+"-"+str(tablekey))
        self.logger.warning("plot_possw_diff: nothing to do for OBSMODE %s",
                            obsmode)
      # end table loop
    show()

