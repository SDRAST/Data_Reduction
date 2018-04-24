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
import numpy
import re

from matplotlib import rcParams
from matplotlib.dates import DateFormatter, epoch2num
from matplotlib.font_manager import FontProperties
from os.path import basename
from pylab import *

from Data_Reduction import trim_extremes
from Data_Reduction.FITS.SDFITSexaminer import DSNFITSexaminer, TidTipAnalyzer
from Data_Reduction.tipping import airmass
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

seconds_formatter = DateFormatter("%H:%M:%S")

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
    mylogger.debug("__init__: initializing superclass")
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
      self.logger = logging.getLogger(parent.logger.name+".Plotter")
      self.logger.debug("__init__: subclass of %s", table)
      for attr in table.__dict__.keys():
        self.__setattr__(attr, table.__getattribute__(attr))
        self.logger.debug("__init__: copying %s", attr)

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
      
    def init_multiplot(self, title, nrows, ncols, size=(4,4), sharey=False):
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
      
      @param size : figure width, figure height (in)
      @type  size : tuple of float
      
      @param sharey : same range on all Y axes (in)
      @type  sharey : bool
      """
      width, height = size
      fig, ax = subplots(nrows=nrows, ncols=ncols, sharey=sharey)
      self.logger.debug("init_multiplot: %d rows, %d columns, size: %s",
                        nrows, ncols, size)
      # could also be fig.set_size_inches(size, forward=True)
      fig.set_size_inches(ncols*width, nrows*height, forward=True)
      fig.subplots_adjust(wspace=0, hspace=0) # no space between plots in a row
      fig.suptitle(title)
      return fig, ax

    def init_multiplots(self, title, nfigs, nrows, ncols, size=(4,4),
                        sharey=False):
      """
      create multiple figures with multiple plots
      
      The sublots have no space between the axes in a figure. Use
      figure_rows_and_columns() to compute number of figures, rows, columns and
      figure size.
      
      @param title : figure title
      @type  title : str
      
      @param nfigs : number of figures
      @type  nfigs : int
      
      @param nrows : number of rows of subplots
      @type  nrows : int
      
      @param ncols : number of columns of subplots
      @type  ncols : int
      
      @param size : figure width, figure height (in)
      @type  size : tuple of float
      
      @param sharey : same range on all Y axes (in)
      @type  sharey : bool
      """
      self.logger.debug("init_multiplots: %d figs, %d rows, %d columns",
                        nfigs, nrows, ncols)
      figs = {}
      axs = {}
      for fignum in range(nfigs):
        figs[fignum], axs[fignum] = self.init_multiplot(title, nrows, ncols,
                                                        size=size,
                                                        sharey=sharey)
      return figs, axs
      
    def show_passband(self, figtitle=None, rows=None, savepath=None):
      """
      Plots the passbands from this file as dynamic spectra
      
      If there are multiple beams, there will be a figure column for each beam
      and each pol.  Else, if there is only one beam but multiple subchannels
      then there will be a figure column for each subchannel and each pol.
      Otherwise there is just a column for each pol.
      """
      if rows == None:
        spectra = self.get_spectra(self.row_keys)
      else:
        spectra = self.get_spectra(rows)
      # lay out the figure
      if self.props['num beams'] > 1:
        ncols = self.props['num beams']*self.props['num IFs']
      elif self.props['num cycles'] > 1:
        ncols = self.props['num cycles']*self.props['num IFs']
      else:
        ncols = self.props['num IFs']
      # collect the diagnostic spectra
      if self.props["full Stokes"]:
        # when SPECTRA are Stokes parameters, plot IF power for two pol modes
        if self.props["num IFs"] == 2:
          datasource = "IFSPECTR"
          num_chans = self.props['num IFspec chans']
          IFspectra = True
      else:
        if "SPECTRUM" in self.data.columns.names:
          datasource = "SPECTRUM"
        else:
          datasource = "DATA"
        IFspectra = False
        num_chans = self.props['num chans']
        self.logger.debug("show_passband: data source is %s", datasource)
      # prepare empty images
      images = self.prepare_summary_images(num_chans)
      # get data statistics for scaling plots
      ymin, ymax, ymean, ystd = self.get_data_stats()
      #   spectra are 2D arrays with frequency and scan as axes
      #   images, also 2D arrays, are only needed if there is a TIME axis in
      #   the data and then each record is a row in the array.
      cycles = self.cycle_keys
      # images with SCAN on the time axis have sub-images over record
      start_image = {}
      for cycle in cycles:
        subch_idx = cycle - 1
        start_image[subch_idx] = {}
        for beam_idx in range(self.props["num beams"]):
          start_image[subch_idx][beam_idx] = {}
          for IF_idx in range(self.props["num IFs"]):
            start_image[subch_idx][beam_idx][IF_idx] = True
      # get the data
      for scan in self.scan_keys:
        scan_idx = self.scan_keys.index(scan) # scan numbers can start anywhere
        for cycle in cycles:
          subch_idx = cycle - 1
          for beam_idx in range(self.props["num beams"]):
            beam = beam_idx+1
            for IF_idx in range(self.props["num IFs"]):
              pol = IF_idx+1
              #self.logger.debug("plot_passband: processing scan %d, subch %d, beam %d, pol %d",
              #               scan, cycle, beam, pol)
              if self.props["time axis"]:
                # average scan spectrum and include record spectra in image
                # assume there is a scan number equal to the cycle number
                image, spectrum = \
                      self.average_records(scan, cycle, beam, pol)
                # this is the all-record sub-image for the scan
                if start_image[subch_idx][beam_idx][IF_idx]:
                  images[subch_idx][beam_idx][IF_idx] = image
                  start_image[subch_idx][beam_idx][IF_idx] = False
                else:
                  images[subch_idx][beam_idx][IF_idx] = \
                      numpy.append(images[subch_idx][beam_idx][IF_idx], image,
                                   axis=0)
              else:
                # no time axis
                spec_indices = self.get_indices(scan=scan, cycle=cycle, 
                                                beam=beam, pol=pol)
                image_line = self.data[datasource][spec_indices].reshape(
                                                       num_chans,1).transpose()
                if start_image[subch_idx][beam_idx][IF_idx]:
                  images[subch_idx][beam_idx][IF_idx] = image_line
                  start_image[subch_idx][beam_idx][IF_idx] = False
                else:
                  images[subch_idx][beam_idx][IF_idx] = \
                              numpy.append(images[subch_idx][beam_idx][IF_idx],
                                           image_line, axis=0)
      if figtitle:
        pass
      else:
        figtitle = basename(self.parent.file)[4:-5].replace('_',' ')
      fig, ax = self.init_multiplot(figtitle+"  Spectogram",
                                       nrows=1, ncols=ncols)
      # display the data
      # labels are updated only in the first time around a loop
      labels = {}
      num_beams = self.props['num beams']
      num_cycles = len(self.cycle_keys)
      num_pols = self.props["num IFs"]
      halfband = self.data['BANDWIDT'][0]/2.e6
      for subch in range(num_cycles):
        for beam in range(self.props["num beams"]):
          for pol in range(self.props["num IFs"]):
            label = make_legend_labels(sckeys=range(num_cycles),
                                       bmkeys=range(num_beams),
                                       plkeys=range(num_pols),
                                       sckey=subch,
                                       bmkey=beam,
                                       plkey=pol)
            if self.props['num beams'] > 1:
              col = 2*beam + pol
            elif self.props['num cycles'] > 1:
              col = 2*subch + pol
            else:
              col = pol
            # dynamic spectra of IF power
            ax[col].set_title(label)
            height, width = images[subch][beam][pol].shape
            ax[col].imshow(images[subch][beam][pol], aspect="auto",
                           extent=(-halfband,halfband,0,height))
            if col == 0:
              ax[col].set_ylabel("Cumulative record number")
              ax[col].set_xlabel("Frequency (MHz)")
            fig.subplots_adjust(top=0.88)
            fig.subplots_adjust(bottom=0.15)
      last_col = len(ax)-1
      lines, labels = ax[last_col].get_legend_handles_labels()
      fig.legend(lines, labels, loc="upper right", ncol=2, prop = fontP)
      if savepath:
        fig.savefig(savepath)
      show()

    def show_all_spectra(self, IFspectra=False, sharey=False):
      """
      Plot spectra from a Table
      
      In each subplot are all the spectra for each beam and polarization from
      one row. If there are multiple records in a row, all records are plotted
      (not the average over records)
      
      @param rows : optional list of table rows; default: all
      @type  rows : list of int
      
      @param IFspectra : plot IFSPECTR if there is one; default: False
      @type  IFspectra : bool
      
      @param sharey : same range on all Y axes (in)
      @type  sharey : bool
      """
      # gets spectra from SPECTRUM column with RA and dec indices removed
      if rows == None:
        spectra = self.get_spectra(self.row_keys)
      else:
        spectra = self.get_spectra(rows)
      num_graphs = len(spectra)/self.props['num cycles']
      nfigs, nrows, ncols, size = self.figure_rows_and_columns(num_graphs)
      figs, axs = self.init_multiplots("Basic Spectra", nfigs, nrows, ncols, 
                                       size, sharey=sharey)
      self.logger.debug("show_all_spectra: figure keys: %s", figs.keys())
      self.logger.debug("show_all_spectra: axes   keys: %s", axs.keys())
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
            self.logger.debug(
                        "show_all_spectra: doing row %d column %d spectrum %d",
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
            nrecs = self.props['num records'][cycle]
            self.logger.debug(
                          "show_all_spectra: %d channels, %d pols, %d records",
                          nchans, npols, nrecs)
            for rec in range(nrecs):
              record = rec+1
              symbol = plotsymbols[rec % 20]
              for beam_idx in range(nbeams):
                beam = beam_idx+1
                for pol_idx in range(npols):
                  pol = pol_idx+1
                  # indices without row (0), RA (-1), dec (-2)
                  if len(spec.shape) == 4:
                    # with beam and time axis
                    indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                               beam=beam, record=record,
                                               trimmed=True)[1:]
                  elif len(spec.shape) == 3:
                    # with time axis
                    indices = self.get_indices(scan=scan, cycle=cycle, pol=pol,
                                               record=record,
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
          
    def plot_PSSW_spectra(self, figtitle=None, scans=[], savepath=None):
      """
      Plot position switched spectra
      
      We need to check self.data['OBSMODE']
      
      @param scans : list of scan numbers, ons and offs
      @type  scans : list of int
      """
      if scans == []:
        scans = self.scan_keys
      # lay out the figure
      if self.props['num beams'] > 1:
        ncols = self.props['num beams']*self.props['num IFs']
      elif self.props['num cycles'] > 1:
        ncols = self.props['num cycles']*self.props['num IFs']
      else:
        ncols = self.props['num IFs']
      if figtitle:
        pass
      else:
        first = basename(self.parent.file).index('_')+1
        figtitle = "DSS-"+str(self.dss) +" "+ \
                                           basename(self.parent.file)[first:-5]
      fig, ax = self.init_multiplot(figtitle+" (Sig-Ref)/Ref",
                                       nrows=1, ncols=ncols)
      # data source
      if 'IFSPECTR' in self.data.names:
        # this is for Stokes spectra
        datasource = 'IFSPECTR'
      elif 'SPECTRUM' in self.data.names:
        datasource = 'SPECTRUM'
      else:
        datasource = 'DATA'
      self.logger.debug("plot_PSSW_spectra: data column is %s", datasource)

      labels = {}
      num_beams = self.props['num beams']
      num_cycles = len(self.cycle_keys)
      num_pols = self.props["num IFs"]
      for subch_idx in range(num_cycles):
        cycle = subch_idx + 1
        for beam_idx in range(self.props["num beams"]):
          beam = beam_idx + 1
          for pol_idx in range(self.props["num IFs"]):
            pol = pol_idx+1
            label = make_legend_labels(sckeys=range(num_cycles),
                                       bmkeys=range(num_beams),
                                       plkeys=range(num_pols),
                                       sckey=subch_idx,
                                       bmkey=beam_idx,
                                       plkey=pol_idx)
            if self.props['num beams'] > 1:
              col = 2*beam_idx + pol_idx
            elif self.props['num cycles'] > 1:
              col = 2*subch_idx + pol_idx
            else:
              col = pol_idx
            # dynamic spectra of IF power
            ax[col].set_title(label)
            self.logger.debug("plot_PSSW_spectra: subch %d beam %d pol %d",
                              cycle, beam, pol)
            on_scans = numpy.unique(self.data['SCAN'][numpy.where(
                                  self.data['SIG'][self.row_keys] == True)[0]])
            of_scans = numpy.unique(self.data['SCAN'][numpy.where(
                                 self.data['SIG'][self.row_keys] == False)[0]])
            num_pairs = min(len(on_scans), len(of_scans))
            for index in range(num_pairs):
              on_scan = on_scans[index]
              of_scan = of_scans[index]
              # get the indices for the DATA cell in the ON position
              try:
                indices = self.get_indices(scan=on_scan, cycle=cycle, pol=pol,
                                           beam=beam)
                self.logger.debug("plot_PSSW_spectra: scan %d on indices: %s",
                                  on_scan, indices)
              except ValueError, details:
                self.logger.warning("plot_PSSW_spectra: %s", str(details))
                continue
              row = indices[0]
              # get the X-axis units for the ON position
              if datasource == 'IFSPECTR':
                v = self.compute_X_axis(row, frame='DELF-OBS',
                                      num_chans=self.props['num IFspec chans'])
              else:
                v = self.compute_X_axis(row, frame='DELF-OBS',
                                      num_chans=self.props['num chans'])
              on  = self.data[datasource][indices]
              # get the indices for the DATA cell in the OFF position
              try:
                ref_indices = self.get_indices(scan=of_scan, cycle=cycle,
                                               pol=pol, beam=beam)
                self.logger.debug("plot_PSSW_spectra: scan %d off indices: %s",
                                  of_scan, ref_indices)
                off =  self.data[datasource][ref_indices]
              except ValueError, details:
                self.logger.warning("plot_PSSW_spectra: %s", str(details))
                continue
              spectrum = (on-off)/off
              label = "("+str(on_scan)+"-"+str(of_scan)+")/"+str(of_scan)
              ax[col].plot(v, spectrum, label=label)
              if index == 0:
                ax[col].grid()
                ax[col].set_xlabel("Frequency (MHz)")
                #heading = "%8.3f MHz" % (self.data['OBSFREQ'][row]/1e6)
                #heading += " P"+str(pol)
                #ax[col].set_title(heading)
      last_col = len(ax)-1
      lines, labels = ax[last_col].get_legend_handles_labels()
      fig.legend(lines, labels, loc="upper right", ncol=2, prop = fontP)
      if savepath:
        savefig(savepath)
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
        titlestr = \
                 source+" DSS-"+str(self.dss)+" "+self.datestr+" "+self.timestr
        title(titlestr)
        xlabel("$V_{LSR}$ (km s$^{-1}$")
        ylabel("$T_{ant}$ (K)")
        if savepath:
          fname = titlestr.replace("/","-").replace(" ","_")+".png"
          savefig(savepath+fname)
        return {"rms": rms, "Tsys": Tsys, "intgr": intgr}
    
    def plot_all_Tsys(self):
      """
      Displays all Tsys values so user can select row range
      
      This works for WVSR data but needs to be elaborated for SAO data
      """
      figure()
      stride = self.props['num cycles']*self.props['num IFs']
      for cycle_idx in range(self.props['num cycles']):
        for beam_idx in range(self.props['num beams']):
          for pol_idx in range(self.props['num IFs']):
            label = self.make_legend_labels(sckey=cycle_idx, bmkey=beam_idx,
                                            plkey=pol_idx)
            for sig in [True, False]:
              if sig:
                lbl = label+" on"
                ls = "-"
                index = cycle_idx
              else:
                lbl = label+" off"
                ls = "--"
                index = 2+cycle_idx
              plot(self.data['TSYS'][index:len(self.row_keys):stride,
                                     pol_idx,0,0,0], ls=ls, label=lbl)
      grid()
      legend(loc='best')
      xlabel("row")
      ylabel("average power")
      title(self.datestr)
      show()
      
    def plot_Tsys(self, good_wx_data=None, X=None):
      """
      Plot average power versus time or airmass or list index
      
      Options for X are::
        "time"    - time of measurement
        "airmass" - 1/sin(elev)
        None      - list index
      """
      if good_wx_data:
        pass
      else:
        good_wx_data = self.get_wx_datacubes()
      heading = "DSS-%2d %s" % (self.dss, self.datestr)
      fig, ax = self.init_multiplot(heading, 1, self.props['num cycles'])
      for subch in self.cycle_keys:
        cycle_idx = subch - 1
        for beam in range(self.props['num beams']):
          for pol in range(self.props['num IFs']):
            for sig in [True, False]:
              tm = good_wx_data['UNIXtime'][sig]
              plottimes = epoch2num(tm)
              tsys = good_wx_data['TSYS'][sig][:,cycle_idx,beam,pol]
              el = good_wx_data['ELEVATIO'][sig]
              label = self.make_legend_labels(bmkey=beam, plkey=pol)
              color_idx = 2*int(sig) + pol
              if sig:
                lbl = label+" sig"
                ls = "-"
              else:
                lbl = label+" ref"
                ls = "--"
              if X == "time":
                ax[cycle_idx].plot_date(plottimes, tsys, linestyle=ls, marker='.',
                                    color=plotcolors[color_idx], label=lbl)
              elif X == "airmass":
                ax[cycle_idx].plot(1/sin(pi*array(el)/180.),
                               tsys, color=plotcolors[color_idx], marker='.',
                               ls=ls, label=lbl)
              else:
                ax[cycle_idx].plot(tsys, color=plotcolors[color_idx], marker='.',
                               ls=ls, label=lbl)
        ax[cycle_idx].grid(True)
        ax[cycle_idx].legend(loc='best', fontsize='xx-small', numpoints=1)
        if X == "time":
          ax[cycle_idx].xaxis.set_major_formatter(seconds_formatter)
          fig.autofmt_xdate()
        elif X == "airmass":
          ax[cycle_idx].set_xlabel("Airmass")
        else:
          ax[cycle_idx].set_xlabel("index")
        ax[cycle_idx].set_title("subchannel "+str(subch))
        fig.show()
                

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

#------------------- class for plotting TIPPING CURVE extension data ----------

class TidTipPlotter(TidTipAnalyzer):
  """
  class for plotting data from the TIPPING CURVE extensions
  """
  def __init__(self, extension):
    """
    """
    TidTipAnalyzer.__init__(self, extension)

  def plot_data():
    """
    """
    fig = figure()
    for IF in range(4):
      PM = IF+1
      rows = where(self.data['CYCLE'] == PM)
      tsys = self.data['TSYS'][rows]
      am = airmass(self.data['ELEVATIO'][rows])
      plot(am, tsys, plotcolors[IF]+'-')
      plot(am, tsys, plotcolors[IF]+'.', label="IF%d" % PM)
      # add fit to legend
      label = "IF%d  %5.1f$\pm$%3.1f  %5.3f$\pm$%5.3f" % (
                 IF, Trx[IF],sigTrx[IF], tau[IF],sigtau[IF])
      handles[IF].set_label(label)
    legend(loc="upper left", numpoints=1)
    grid()
    title(self.header['DATE-OBS'])
    xlabel("airmass")
    legend(loc="upper left", numpoints=1)
    show()
    if project:
      sessiondir = projects_dir+project+"/Observations/dss43/%4d/%03d/" % (
                                                           self.year, self.DOY)
      fig.savefig(sessiondir+"tipdatafit-"+str(index+1)+".png")

  
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
  
def get_power_range(table, subch, beam, pol):
            """
            calculate limits for power levels
            
            the idea here is to avoid huge spikes compressing the spectra
            will need some tuning
            
            @param subch : sub-channel number
            @param beam : beam number (1-based)
            @param pol   : polarization number (not NRAO pol code)
            """
            subch_idx = subch-1
            beam_idx = beam-1
            IF_idx = pol-1
            # gets the indices for scan 1
            indices = table.get_indices(cycle=subch, beam=beam, IF=pol)
            logger.debug("get_power_range: indices: %s", indices)
            # gets spectra from all rows
            spectrum = table.data['SPECTRUM'][:,indices[1:]]
            logger.debug("get_power_range: spectrum shape is %s", spectrum.shape)
            mean_pwr = spectrum.mean()
            pwr_std = spectrum.std()
            max_pwr = spectrum.max()
            min_pwr = spectrum.min()
            if max_pwr > mean_pwr + 4*pwr_std:
              ymax = mean_pwr + pwr_std
            else:
              ymax = mean_pwr + 4*pwr_std
            if min_pwr < mean_pwr - 4*pwr_std:
              ymin = mean_pwr - pwr_std
            else:
              ymin = mean_pwr - 4*pwr_std
            return ymin, ymax
