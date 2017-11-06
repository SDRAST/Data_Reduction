"""
Provides a subclass to DSNFITSexaminer for plotting data.

Possible symbols::
    ================    ===============================
    character           description
    ================    ===============================
    ``'o'``             circle marker
    ``'v'``             triangle_down marker
    ``'^'``             triangle_up marker
    ``'<'``             triangle_left marker
    ``'>'``             triangle_right marker
    
    ``'s'``             square marker
    ``'p'``             pentagon marker
    ``'D'``             diamond marker
    ``'h'``             hexagon1 marker
    ``'H'``             hexagon2 marker
    
    ``'1'``             tri_down marker
    ``'2'``             tri_up marker
    ``'3'``             tri_left marker
    ``'4'``             tri_right marker
    ``'*'``             star marker
    
    ``'+'``             plus marker
    ``'x'``             x marker
    ``'d'``             thin_diamond marker
    ``'|'``             vline marker
    ``'_'``             hline marker
    
    ``'.'``             point marker
    ``','``             pixel marker
    ``'-'``             solid line style
    ``'--'``            dashed line style
    ``'-.'``            dash-dot line style
    ``':'``             dotted line style
    ================    ===============================
"""

import logging
import numpy as np
import re

from matplotlib.font_manager import FontProperties
from os.path import basename
from pylab import *

from Data_Reduction import trim_extremes
from Data_Reduction.FITS.SDFITSexaminer import DSNFITSexaminer
from DatesTimes import UnixTime_to_MPL
from support.graphics import get_screen_resolution
from support.text import clean_TeX

logger = logging.getLogger(__name__)
plotcolors = ['b','g','r','m','c']

plotsymbols = ['o','v','^','<','>',
               's','p','D','h','H',
               '1','2','3','4','*',
               '+',"x","d","|","_"]

fontP = FontProperties()
fontP.set_size('x-small')      

class DSNFITSplotter(DSNFITSexaminer):
  """
  """
  def __init__(self, parent=None, FITSfile=None, hdulist=None):
    """
    Create a new DSNFITSplotter object from an HDU list
    
    If invoked from within another object, then parent should be 'self'.
    Either a FITS file or an HDU list must be provided.
    """
    mylogger = logging.getLogger(logger.name+".DSNFITSplotter")
    DSNFITSexaminer.__init__(self, parent=parent, FITSfile=FITSfile,
                             hdulist=hdulist)
    self.logger = mylogger
    self.plotter = {}
    for key in self.tables.keys():
      table = self.tables[key]
      self.logger.debug("__init__: processing %s", table)
      self.plotter[key] = self.Plotter(self, table)
    self.logger.debug("__init__: completed")
  
  class Plotter(DSNFITSexaminer.Table):
    """
    """
    def __init__(self, parent, table):
      """
      Initialization was already done for Table superclass
      """
      for attr in table.__dict__.keys():
        self.__setattr__(attr, table.__getattribute__(attr))
      self.logger = logging.getLogger(parent.logger.name+".Plotter")

    #------------------------ Plotter methods ---------------------------------
    
    def figure_rows_and_columns(self, naxes, **kwargs):
      """
      Computes number of rows and columns for 'subplots'
      
      @param naxes : number of subplots (axes) needed
      @type  naxes : int
      """
      # get maximum number of rows and columns
      screensize = get_screen_resolution()
      widthIN, heightIN  = screensize['width']/100., screensize['height']/100.
      freewidth, freeheight = 0.95*widthIN, 0.95*heightIN
      self.logger.debug("figure_rows_and_columns: width available = %f in",
                        freewidth)
      self.logger.debug("figure_rows_and_columns: height available = %f in",
                        freeheight)
      if kwargs.has_key("width"):
        width = kwargs["width"] # inches
      else:
        width = 4
      if kwargs.has_key("heigth"):
        height = kwargs["width"] # inches
      else:
        height = 4
      max_cols = int(round(freewidth/width))
      max_rows = int(round(freeheight/height))
      self.logger.debug("figure_rows_and_columns: max cols and rows: %d,%d",
                        max_cols, max_rows)
      max_graphs = max_cols*max_rows
      aspect = float(max_cols)/max_rows
      # how many figures do we need?
      num_figs = 1+naxes/max_graphs
      if num_figs > 1:
        num_cols = max_cols
        num_rows = max_rows
      else:
        num_cols = int(ceil(sqrt(aspect*naxes)))
        num_rows = int(ceil(float(naxes)/num_cols))
      self.logger.debug("figure_rows_and_columns: %d rows, %d columns",
                        num_rows, num_cols)
      return num_figs, num_rows, num_cols, (width,height)
      
    def init_multiplot(self, title, nrows, ncols, size=(4,4)): 
      """
      create a figure with multiple plots sharing common X and Y axis
      
      The sublots have no space between them.
      
      When the number of graphs is large then multiple figures need to be 
      created with 'init_multiplots'
      
      @param title : figure title
      @type  title : str
      
      @param nrows : number of rows of subplots
      @type  nrows : int
      
      @param ncols : number of columns of subplots
      @type  ncols : int
      
      @param width : optional keyword argument for figure width (in)
      @type  width : float
      
      @param heigth : optional keyword argument for figure heigth (in)
      @type  height : float
      """
      width, height = size
      fig, ax = subplots(nrows=nrows, ncols=ncols) # , sharex=True, sharey=True)
      self.logger.debug("init_multiplot: %d rows, %d columns, size: %s",
                        nrows, ncols, size)
      # could also be fig.set_size_inches(size, forward=True)
      fig.set_size_inches(ncols*width, nrows*height, forward=True)
      fig.subplots_adjust(wspace=0, hspace=0) # no space between plots in a row
      fig.suptitle(title)
      return fig, ax

    def init_multiplots(self, title, nfigs, nrows, ncols, size=(4,4)):
      """
      create multiple figures with multiple plots
      """
      self.logger.debug("init_multiplots: %d figs, %d rows, %d columns",
                        nfigs, nrows, ncols)
      figs = {}
      axs = {}
      for fignum in range(nfigs):
        figs[fignum], axs[fignum] = self.init_multiplot(title, nrows, ncols,
                                                        size=size)
      return figs, axs
      
    def show_all_spectra(self, rows=None, IFspectra=False):
      """
      Plot spectra from a Table
      
      In each subplot are all the spectra for each beam and polarization from
      one row. If there are multiple records in a row, all records are plotted
      (not the average over records)
      
      @param rows : optional list of table rows; default: all
      @type  rows : list of int
      
      @param IFspectra : plot IFSPECTR if there is one; default: False
      @type  IFspectra : bool
      """
      # gets spectra from SPECTRUM column with RA and dec indices removed
      if rows == None:
        spectra = self.get_spectra(self.row_keys)
      else:
        spectra = self.get_spectra(rows)
      num_graphs = len(spectra)/self.props['num cycles']
      nfigs, nrows, ncols, size = self.figure_rows_and_columns(num_graphs)
      figs, axs = self.init_multiplots("Basic Spectra",nfigs, nrows, ncols, size)
      self.logger.debug("show_spectra: figure keys: %s", figs.keys())
      self.logger.debug("show_spectra: axes   keys: %s", axs.keys())
      for fignum in figs.keys():
        fig = figs[fignum]
        ax = axs[fignum]
        for row in range(nrows):
          for col in range(ncols):
            specidx = nrows*ncols*fignum + ncols*row + col
            if specidx >= len(spectra):
              break
            scan = self.data['SCAN'][specidx]
            cycle = self.data['CYCLE'][specidx]
            cycle_idx = self.cycle_keys.index(cycle)
            self.logger.debug("show_spectra: doing row %d column %d spectrum %d",
                              row, col, specidx)
            spec = spectra[specidx] # returns spectra from one row
            if self.props['full Stokes']:
              if "IFSPECTR" in self.data.columns.names and IFspectra:
                nchans = self.props['num IFspec chans']
                npols = self.props['num IFs']
              else:
                nchans = self.props['num chans']
                npols = spec.shape[0]
            else:
              nchans = self.props['num chans']
              npols = self.props['num IFs']
            nbeams = self.props['num beams']
            self.logger.debug("show_spectra: %d channels, %d pols", nchans, npols)
            nrecs = self.props['num records'][cycle]
            for rec in range(nrecs):
              record = rec+1
              symbol = plotsymbols[rec % 20]
              for beam_idx in range(nbeams):
                beam = beam_idx+1
                for pol_idx in range(npols):
                  pol = pol_idx+1
                  # indices without row (0), RA (-1), dec (-2)
                  if len(spec.shape) == 4:
                    indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                               beam=beam, record=record,
                                               trimmed=True)[1:]
                  elif len(spec.shape) == 2:
                    indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                               trimmed=True)[1:]
                  self.logger.debug("show_all_spectra: indices: %s", indices)
                  color = plotcolors[beam_idx*npols+pol_idx % 5]
                  #trimmed = trim_extremes(spec[indices])
                  trimmed = spec[indices]
                  label = self.make_legend_labels(sckey=cycle_idx,
                                                  bmkey=beam_idx,
                                                  plkey=pol_idx)
                  ax[row][col].plot(trimmed, color+',', label=label)
            # end loop over records
            ax[row][col].grid(True)
            ax[row][col].text(0.5, 0.95, 'scan '+str(scan),
                              transform=ax[row][col].transAxes,
                              horizontalalignment='center', 
                              verticalalignment='top')
            if row == nrows-1:
              for tick in ax[row][col].get_xticklabels():
                tick.set_rotation(45)
            if col == 0:
              ax[row][col].set_ylabel("Power (counts)")
          # end loop over col
        # end loop over row
        lines, labels = ax[0][0].get_legend_handles_labels()
        fig.legend(lines, labels, loc="upper right", ncol=2, prop = fontP)
      # end loop over fig
      show()
    
    def make_legend_labels(self, dskey=None, tbkey=None, sckey=None,
                           bmkey=None, plkey=None):
      dskeys=[]
      tbkeys=[]
      sckeys = range(self.props['num cycles'])
      bmkeys = range(self.props['num beams'])
      if self.props['full Stokes']:
        plkeys = range(4)
      else:
        plkeys = range(self.props['num IFs'])
      return make_legend_labels(dskeys=dskeys, tbkeys=tbkeys, sckeys=sckeys,
                                bmkeys=bmkeys, plkeys=plkeys,
                                dskey=dskey, tbkey=tbkey, sckey=sckey,
                                bmkey=bmkey, plkey=plkey)
                                
    def plot_BPSW_spectra(self, spectra=None, rows=None):
      """
      plot reduced beam and position switched spectra
      
      @param spectra : optional dict of reduced spectra from 'BPSW_spectra()'
      @type  spectra : numpy array
      
      @param rows : optional list of table rows to compute spectra; deault: all
      @type  rows : list of int
      """
      if spectra:
        npairs = len(spectra)
      elif rows :
        rows = self.row_keys
        spectra = self.BPSW_spectra(rows)
        npairs = len(spectra)
      else:
        rows = self.row_keys
        spectra = self.BPSW_spectra(rows)
        npairs = len(spectra)
      pairs = range(npairs)
      nrows, ncols = self.figure_rows_and_columns(len(spectra))
      fig, ax = self.init_multiplot("BPSW Spectra", nrows, ncols)
      for row in range(nrows):
        for col in range(ncols):
          pair = row*ncols + col
          spec = spectra[pair]
          nrecs, npols, nchans = spec.shape
          for rec in range(nrecs):
            symbol = plotsymbols[rec % 20]
            for pol in range(npols):
              color = plotcolors[pol % 5]
              trimmed = trim_extremes(spec[rec, pol])
              ax[row][col].plot(trimmed, color+',',
                          label="rec"+str(rec)+" P"+str(pol+1))
          ax[row][col].grid(True)
          ax[row][col].text(0.5, 0.95, 'pair '+str(pair+1),
                            transform=ax[row][col].transAxes,
                            horizontalalignment='center', 
                            verticalalignment='top')
          if row == nrows-1:
            for tick in ax[row][col].get_xticklabels():
              tick.set_rotation(45)
          if col == 0:
            ax[row][col].set_ylabel("Normalized Power")
      lines, labels = ax[0][0].get_legend_handles_labels()
      fig.legend(lines, labels, loc="upper right", ncol=2, prop = fontP)
      show()
          
    def plot_PSSW_spectra(self, scans=[]):
      """
      Plot position switched spectra
      
      We need to check self.data['OBSMODE']
      
      @param scans : list of scan numbers, ons and offs
      @type  scans : list of int
      """
      if scans == []:
        scans = self.scan_keys
      for cycle in self.cycle_keys:
        figure()
        for scan in scans[::2]:
          try:
            indices = self.get_indices(scan, cycle)
          except ValueError, details:
            self.logger.warning("plot_PS_spectra: %s", str(details))
            continue
          row = indices[0]
          v = self.compute_X_axis(row, frame='RADI-LSR')
          on  = self.data['SPECTRUM'][indices]
          try:
            off =  self.data['SPECTRUM'][self.get_indices(scan+1, cycle)]
          except ValueError, details:
            self.logger.warning("plot_PS_spectra: %s", str(details))
            continue
          spectrum = (on-off)/off
          label = "("+str(scan)+"-"+str(scan+1)+")/"+str(scan+1)
          plot(v, spectrum, label=label)
        grid()
        xlabel("$V_{lsr}$ (km/s)")
        legend()
        heading = "DSS-"+str(self.dss)+" "+self.data['DATE-OBS'][row]
        heading += " %8.3f MHz" % (self.data['OBSFREQ'][row]/1e6)
        title(heading)
        show()

    def plot_line(self, rows=[], window=(-100,100),
                    frame="RADI-OBJ", source='67P', savepath=None):
      """
      plot reduced averaged spectral lines (from 'reduce_line()) for both pols
      
      @param rows : optional rows to include; default: all
      @type  rows : list of int
      
      @param window : optional left and right limits of the X axis (-100,100)
      @type  window : tuple of float
      
      @param frame : optional reference frame for Doppler calc ("RADI-OBJ")
      @type  frame : str
      
      @param source : source name for Doppler calculation; default: 67P
      @type  source : str
      
      @param savepath : optional path for saved image; default: no save
      @type  savepath : str
      """
      try:
        x, y, rms, Tsys, intgr = self.reduce_line(rows=rows, window=window,
                                   frame=frame, source=source)
      except TypeError:
        self.logger.error("plot_line: nothing to plot")
      else:
        figure()
        plot(x, y[0], label="Pol1")
        plot(x, y[1], label="Pol2")
        grid()
        xlim(*window)
        legend(prop=fontP)
        titlestr = source+" DSS-"+str(self.dss)+" "+self.datestr+" "+self.timestr
        title(titlestr)
        xlabel("$V_{LSR}$ (km s$^{-1}$")
        ylabel("$T_{ant}$ (K)")
        if savepath:
          fname = titlestr.replace("/","-").replace(" ","_")+".png"
          savefig(savepath+fname)
        return {"rms": rms, "Tsys": Tsys, "intgr": intgr}
    
    def plot_Tsys(self):
      """
      Plots of system tempearture for all subchannels, beams and pols
      
      If there is no time axis then there is a time and a Tsys for every scan,
      cycle and pol, that is, four plots over rows.
      
      If there is a time axis then there is also a time and Tsys for every 
      record
      """
      # Output the run of system temperatures as a diagnostic:
      for cycle in self.cycle_keys:
        cycle_idx = self.cycle_keys.index(cycle)
        fig, ax = subplots(self.props['num beams'],self.props['num IFs'])
        cindex = 0 # count for color selection
        for beam in range(self.props['num beams']):
          for pol in range(self.props['num IFs']):
            # just get the shape of the data cube
            indices = self.get_indices()
            if len(indices) == 6:
              # row, beam, record, IF, dec, RA, time
              plottimes = \
               UnixTime_to_MPL(self.data['UNIXtime'][:,beam,:,cycle_idx,0,0,0])
              tsys = self.data['TSYS'][:,beam,:,cycle_idx,0,0,0]
            elif len(indices) == 4:
              # row, IF, dec, RA, time
              plottimes = UnixTime_to_MPL(self.data['UNIXtime'][:])
              tsys = self.data['TSYS'][:,cycle_idx,0,0,0]
            else:
              self.logger.error("plot_Tsys: %d axes is invalid", len(indices))
              raise RuntimeError("invalid number of data axes")
            # draw the line
            ax[beam][pol].plot_date(plottimes[:len(tsys)], tsys, linestyle='-')
            # show the samples
            if beam == 0 and pol == 0:
              ax[beam][pol].plot_date(plottimes, tsys,
                        label="Pol "+str(pol) + " Beam "+str(beam))
            else:
              ax[beam][pol].plot_date(plottimes, tsys)
            cindex += 1
            ylabel(r"T$_{sys}$ (K)")
            legend(loc='best', fontsize='xx-small', numpoints=1)
            grid()
            titlestr = clean_TeX(self.data['DATE-OBS'][0])
            title(titlestr)
        fig.autofmt_xdate()
    
    def plot_Tsys_vs_am(self, axes):
      """
      plot system temperature versus airmass
      
      @param axes : subplot to use
      @type  axes : matplotlib.axes.Axes object
      """
      good_wx_data = self.get_good_rows()
      if good_wx_data.has_key('elev') and good_wx_data.has_key('TSYS'):
        for subch in range(self.props['num cycles']):
          for beam in range(self.props['num beams']):
            for pol in range(self.props['num IFs']):
              label = make_legend_labels(
                                        sckeys=range(self.props['num cycles']),
                                        bmkeys=range(self.props['num beams']),
                                        plkeys=range(self.props['num IFs']),
                                        sckey=subch, bmkey=beam, plkey=pol)
              color_index = self.props['num beams']*(self.props['num cycles'] \
                           *subch + beam) + pol
              axes.plot(1/sin(pi*array(good_wx_data['elev'])/180.)[
                                             subch::self.props['num cycles']],
                        good_wx_data['TSYS'][subch][beam][pol],
                        color=plotcolors[color_index],
                        label=label)
      axes.set_xlabel("Airmass")
      axes.grid()

  #--------------------------- DSNFITSplotter methods -------------------------
  
  def plot_average(self, frame='RADI-LSR', source=None,
                   xlimits=None, ylimits=None):
    """
    plot an averaged DSN FITS spectrum
    
    A DSN FITS spectrum is multi-dimensional with axes::
      [[beam,] [time,]], pol, [dec, RA], frequency-like
    [] indicates the axes may not be present.  During data manipulations the
    [ra, dec,] axes are the first to be eliminated.
    
    Note that the frame of the provided spectrum is not changed.  Only the
    X-axis is recomputed before plotting.
    
    @param row : row to be used to get observing parameters
    @type  row : int
    
    @param frame : frame in which to plot the data
    @type  frame : str
    
    @param xlimits : minimum and maximum X-axis values
    @type  xlimits : (float, float)
    
    @param ylimits : minimum and maximum Y-axis values
    @type  ylimits : (float, float)
    
    @param average : combine the two polarizations
    @type  average : bool
    
    """
    
    if type(spectrum) == DSNFITSexaminer.Table.Data:
      # this is a window into the original spectrum
      x = spectrum.x
      y = spectrum.y
      frame = spectrum.frame
      # change frame if desired
      if frame == spectrum.frame:
        self.logger.debug("plot_spectrum: keeping frame %s", frame)
      else:
      # change to the requested frame
        self.logger.debug("plot_spectrum: changing to frame %s", frame)
        new_x = spectrum.compute_X_axis(row, frame)
        x = new_x[window.channels[0]:window.channels[1]]
        self.logger.debug("plot_spectrum: plotting in frame %s", self.frame)
    else:
      # get SPECTRUM from current row
      self.logger.info("plot_spectrum: data is not a class Data instance")
      spectrum = self.get_spectra(row)
      if xlimits:
        ch1, ch2 = xlimits[0], xlimits[1]
      if frame == "RADI-OBJ" and self.frame != "RADI-OBJ":
        # Change to the object's rest frame
        source = self.dataset[row]["OBJECT"]
        if re.match('67P', source):
          self.logger.debug("plot_spectrum: using 67P frame")
          vspline = self.load_ephemeris('67P')
          x = spectrum.compute_X_axis(row=row, frame="RADI-OBJ",
                                                      vspline=vspline)[ch1:ch2]
          y = spectrum[ch1:ch2]
        else:
          self.logger.error("plot_spectrum: no ephemeris for %s", self.object)
          raise RuntimeException("cannot compute velocity of %s" % source)
      else:
        x = spectrum.compute_X_axis(row=row)[ch1:ch2]
        y = spectrum[ch1:ch2]
    self.logger.debug("plot_spectrum: x = %s", x)
    self.logger.debug("plot_spectrum: y = %s", y)
    
    # make the plot
    figure()
    rms = {}
    for pol in [0,1]:
      plot(x, y[pol], label="pol "+str(pol)+" ("+("%d"%self.tau[pol])+"s)")
      rms[pol] = y[pol].std()
    if average:
      ave = (y[0]+y[1])/2
      plot(x, ave, label="both ("+("%d"%(self.tau[0]+self.tau[1]))+"s)")
      rms['avg'] = ave.std()
    if xlimits:
      xlim(*xlimits)
    if ylimits:
      ylim(*ylimits)
    if self.frame == "CHAN-OBS":
      xlabel("Channel")
    elif self.frame[:5] == ("OPTI-" or "RADI-" or "RELA-"):
      if self.frame[4:] == "-OBS":
        xlabel(r"$V_{obs} (\mbox{km s}^{-1})$")
    else:
      xlabel("Frequency (MHz)")
    ylabel(r"Antenna Temperature (K)")
    grid()
    legend()
    titlestr = clean_TeX(str(ds0.year)+"/"+str(ds0.doy))
    title(titlestr)
    show()
    self.logger.info("plot_spectrum: pol0, pol1, avg: %s", rms)
    return rms
    
  def plot_Tsys(self):
    """
    """
    # Output the run of system temperatures as a diagnostic:
    if self.dataset == {}:
      self.get_datasets()
    ds0 = self.dataset[0]
    fig = figure()
    for dskey in self.dataset.keys():
      ds = self.dataset[dskey]
      self.logger.debug("plot_Tsys: dataset %d scan keys: %s", dskey, ds.scan_keys)
      for key in ds.scan_keys:
        index = ds.scan_keys.index(key)
        plottimes = UnixTime_to_MPL(ds.header['time'][index])
        cindex = 0
        for pol in [0,1]:
          for beam in [0,1]:
            tsys = ds.header['TSYS'][index][pol,:len(plottimes),beam]
            plot_date(plottimes[:len(tsys)], tsys, linestyle='-', 
                      color=plotcolors[cindex], marker=plotsymbols[dskey])
            if key == ds.scan_keys[0]:
              plot_date(plottimes[0], tsys[0],
                        color=plotcolors[cindex],
                        marker=plotsymbols[dskey],
                        label=clean_TeX(basename(ds.file))+", Pol "+str(pol) \
                        +" Beam "+str(beam))
            else:
              plot_date(plottimes[0], tsys[0], color=plotcolors[cindex],
                        marker=plotsymbols[dskey])
            cindex += 1
    ylabel(r"T$_{sys}$ (K)")
    legend(loc='best', fontsize='xx-small', numpoints=1)
    grid()
    titlestr = clean_TeX(str(ds0.year)+"/"+str(ds0.doy))
    title(titlestr)
    fig.autofmt_xdate()

#----------------------------------- module functions -------------------------

def make_legend_labels(dskeys=[], tbkeys=[], sckeys=[], bmkeys=[], plkeys=[],
                       dskey=None, tbkey=None, sckey=None, bmkey=None, plkey=None):
  """
  @param dskeys : all datafile or examiner keys
  @param tbkeys : all table keys
  @param sckeys : all subchannel keys
  @param bmkeys : all beam keys
  @param plkeys : all polarization keys
  @param dskey : datafile or examiner key
  @param tbkey : table key
  @param sckey : subchannel key
  @param bmkey : beam key
  @param plkeys :polarization key
  """
  label = ""
  if dskey != None and len(dskeys) > 1:
    label += "ds"+str(dskey+1)
  if tbkey != None and len(tbkeys) > 1:
    label += " tb"+str(tbkey+1)
  if sckey != None and len(sckeys) > 1:
    label += " sc"+str(sckey+1)
  if bmkey != None and len(bmkeys) > 1:
    label += " B"+str(bmkey+1)
  if plkey != None and len(plkeys) > 1:
    label += "P"+str(plkey+1)
  return label

