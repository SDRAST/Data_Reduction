"""
"""
import logging
import numpy as np
import re

from pylab import *
from os.path import basename

from DatesTimes import UnixTime_to_MPL
from Data_Reduction.DSN.SAO import Data, SAOexaminer
from Math.least_squares import savitzky_golay
from support.text import clean_TeX

logger = logging.getLogger(__name__)
plotcolors = ['b','g','r','m','c']
plotsymbols = ['o','v','^','<','>','s','p','D']
    
class SAOplotter(SAOexaminer):
  """
  """
  def __init__(self, project=None, dss=None, date=None):
    """
    """
    mylogger = logging.getLogger(logger.name+".Plotter")
    SAOexaminer.__init__(self, project=project, dss=dss, date=date)
    self.logger = mylogger
    self.logger.debug("__init__: completed")
      
  def plot_spectrum(self, data=None, frame=None, xlimits=None, ylimits=None,
                    average=True):
    """
    @param data : structure with dual-pol data
    @type  data : 
    """
    ds0 = self.dataset[0]
    # change frame if desired
    if frame == "RADI-OBJ" and self.frame == "RADI-OBJ":
      self.logger.debug("plot_spectrum: keeping frame RADI-OBJ")
    elif frame:
      # change to the requested frame
      self.logger.debug("plot_spectrum: changing to frame %s", frame)
      x = ds0.compute_X_axis(frame, ds0.dss)
      if data:
        data.frame = frame
    else:
      # no frame specified
      frame = "CHAN-OBS"
      if data:
        data.frame = self.frame
        self.logger.debug("plot_spectrum: frame of extract is %s", data.frame)
      x = ds0.compute_X_axis(frame, ds0.dss)
    self.frame = frame
    self.logger.debug("plot_spectrum: keeping frame %s", self.frame)
    
    # provide data or use avg_diff as default
    if data == None:
      # first get current frame, X-values and Y-values
      x = ds0.compute_X_axis(self.frame, ds0.dss)
      y = self.avg_diff
    elif type(data) != Data:
      self.logger.error("plot_spectrum: data is not a class Data instance")
    else:
      x = data.x
      y = data.y
      frame = data.frame
      if frame == "RADI-OBJ" and self.frame != "RADI-OBJ":
        # Change to the object's rest frame
        if re.match('67P', ds0.source):
          self.logger.debug("plot_spectrum: using 67P frame")
          vspline = self.load_ephemeris('67P')
          ch1,ch2 = data.channels
          x = ds0.compute_X_axis("RADI-OBJ", ds0.dss,
                                 vspline=vspline)[ch1:ch2]
          y = data.y
        else:
          self.logger.error("plot_spectrum: no ephemeris for %s", self.object)
    self.logger.debug("plot_spectrum: x = %s", x)
    self.logger.debug("plot_spectrum: y = %s", y)
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
    for dskey in list(self.dataset.keys()):
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
    
    
     
  
  
