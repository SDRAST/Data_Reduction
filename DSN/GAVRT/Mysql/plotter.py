"""
plots data from the DSS-28 EAC database
"""
import logging
import os
import pickle

from pylab import *

from Data_Reduction.DSN.GAVRT.Mysql.dss28db import DSS28db, Map, Session
from Data_Reduction.maps import plot_xdec_dec
logger = logging.getLogger(__name__)

_host,_user,_pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p",
                                     "rb" ))

class DBPlotter(DSS28db):
  """
  """
  #def __init__(self, host=_host, user=_user, pw=_pw,
  #                   name='dss28_eac', port=3306):
  def __init__(self):
    """
    """
    mylogger = logging.getLogger(logger.name+".DBPlotter")
    #DSS28db.__init__(self, host=host, user=user, pw=pw)
    DSS28db.__init__(self)
    self.logger = mylogger
 
  def get_session_plotter(self, year, doy):
    """
    get IDs for an observing session
    """
    if self.sessions.has_key(year) == False:
      self.sessions[year] = {}
    self.sessions[year][doy] = SessionPlotter(self, year, doy)
    return self.sessions[year][doy]


class SessionPlotter(Session):
  """
  """
  def __init__(self, parent, year, doy):
    """
    """
    mylogger = logging.getLogger(logger.name+".SessionPlotter")
    Session.__init__(self, parent, year, doy, plotter=True)
    self.logger = mylogger
    self.get_map_plotters()
   
  def get_map_plotters(self):
    """
    Returns the raster configuration IDs for the specified date
    """
    map_cfg_ids = self.db.get(
                         "select raster_cfg_id from raster_cfg where year = " +
                         str(self.year) + " and doy = " + str(self.doy) + 
                         ";")
    self.logger.debug("get_maps_plotters: map IDs: %s", map_cfg_ids)
    self.maps = {}
    for map_id in map_cfg_ids[:,0]:
      self.logger.debug("get_maps_plotters: getting %d", map_id)
      self.maps[map_id] = MapPlotter(self, map_id)
    return self.maps
    
    
class MapPlotter(Map):
  """
  """
  def __init__(self, parent, raster_cfg_id, name=None):
    """
    initialize a GAVRT data plotter
    """
    mylogger = logging.getLogger(logger.name+".MapPlotter")
    Map.__init__(self, parent, raster_cfg_id, name=name)
    self.logger = mylogger
  
  def plot_xdec_dec(self):
    """
    """
    UNIXtime, xdecs, decs, tsrc = self.get_raster_data()
    plot_xdec_dec(self.raster_data['xdec'],
                  self.raster_data['dec'], self.name)
    
  def contours(self, channel, width=1.0, height=1.0, contours=None):
    """
    """
    x,y,z = self.regrid(width=width, height=height)
    fig, ax = subplots(nrows=1, ncols=1)
    if contours:
      contour( x, y, z[channel], contours, linewidths=0.5, colors='k')
      contourf(x, y, z[channel], contours, cmap=plt.cm.jet)
    else:
      contourf(x, y, z[channel], cmap=plt.cm.jet)
    grid()
    ax.set_aspect('equal')
    ax.set_xlim(-width/2., width/2.)
    ax.set_ylim(-height/2., height/2.)
    ax.set_xlabel("Cross-declination offset")
    ax.set_ylabel("Declination offset")
    ax.set_title(self.name+" at " + str(self.rss_cfg[channel]['sky_freq']) + 
                 " MHz " + self.rss_cfg[channel]['pol'][0].upper())
    colorbar(shrink=0.6) # draw colorbar
    hh = self.cfg['utc'].seconds/3600
    mm = (self.cfg['utc'].seconds-3600*hh)/60
    fig.suptitle(str(int(self.cfg['year'])) +
                 "/"+str(int(self.cfg['doy'])) + (" %02d:%02d UT" % (hh,mm)))
    
    
    
