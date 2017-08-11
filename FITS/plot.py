"""
Provides a subclass to DSNFITSexaminer for plotting data.
"""
import logging
import numpy as np
import re

from pylab import *
from os.path import basename

from Data_Reduction.FITS.SDFITSexaminer import Data, DSNFITSexaminer
from DatesTimes import UnixTime_to_MPL
from Math.least_squares import savitzky_golay
from support.text import clean_TeX

logger = logging.getLogger(__name__)
plotcolors = ['b','g','r','m','c']
plotsymbols = ['o','v','^','<','>','s','p','D']
    
class DSNFITSplotter(DSNFITSexaminer):
  """
  """
  def __init__(self, FITSextension):
    """
    """
    mylogger = logging.getLogger(logger.name+".Plotter")
    DSNFITSexaminer.__init__(self, FITSextension)
    self.logger = mylogger
    self.logger.debug("__init__: completed")
      
  def plot_spectrum(self, spectrum=None, row=0, frame='RADI-LSR', xlimits=None, 
                          ylimits=None, average=False):
    """
    plot a DSN FITS spectrum
    
    A DSN FITS spectrum is multi-dimensional with axes::
      [[beam,] [time,]], pol, [dec, RA], frequency-like
    [] indicates the axes may not be present.  During data manipulations the
    [ra, dec,] axes are the first to be eliminated.
    
    Note that the frame of the provided spectrum is not changed.  Only the
    X-axis is recomputed before plotting.
    
    @param spectrum : raw or reduced spectra
    @type  spectrum : multi-dimensional SPECTRUM of DSN type
    
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
    self.logger.debug("plot_spectrum: spectrum is type %s", type(spectrum)
    if type(spectrum) == Data:
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
        x = 
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
    
    
     
  
  
