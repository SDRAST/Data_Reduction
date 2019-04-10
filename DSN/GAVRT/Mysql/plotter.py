"""
plots data from the DSS-28 EAC database
"""
import logging
import os
import pickle
import time

from pylab import *

from Data_Reduction.DSN.GAVRT.Mysql.dss28db import DSS28db, Map, Session
from Data_Reduction.maps import plot_xdec_dec
logger = logging.getLogger(__name__)

_host,_user,_pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p",
                                     "rb" ))

class DBPlotter(DSS28db):
  """
  Attributes::
    logger - logging.Logger object
  
  Attributes inherited from DSS28db
    receiver - receivers which provide data
    sessions - dict of sessions obtained with 'get_session'
  """
  def __init__(self):
    """
    """
    mylogger = logging.getLogger(logger.name+".DBPlotter")
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
    obs_dir = "/usr/local/projects/SolarPatrol/Observations/dss28/"
    self.session_dir = obs_dir + "%4d" % self.year +"/"+ "%03d" % self.doy +"/"
    if not os.path.exists(self.session_dir):
      os.makedirs(self.session_dir, mode=0777)
    self.get_map_plotters()
      
  def get_map_plotters(self):
    """
    Returns the raster configuration IDs for the specified date
    """
    map_cfg_ids = self.db.get(
                         "select raster_cfg_id from raster_cfg where year = " +
                         str(self.year) + " and doy = " + str(self.doy) + 
                         ";")
    self.logger.debug("get_map_plotters: cfg IDs: %s", map_cfg_ids)
    self.maps = {}
    if map_cfg_ids.shape == (0,):
      self.logger.error("get_map_plotters: no maps for this session")
    else:
      self.logger.debug("get_maps_plotters: map IDs: %s", map_cfg_ids[:,0])
      for map_id in map_cfg_ids[:,0]:
        self.logger.debug("get_maps_plotters: getting %d", map_id)
        self.maps[map_id] = MapPlotter(self, map_id)
    return self.maps

  def plot_centered_offsets(self, mapkeys=None):
    """
    show the positions at which data samples were taken
    
    @param mapkeys : numbers of the maps (default: all)
    @type  mapkeys : list of int
    """
    if mapkeys:
      self.logger.info("plot_centered_offsets: processing: %s", mapkeys)
    else:
      mapkeys = self.maps.keys()
      mapkeys.sort()
      self.logger.info("plot_centered_offsets: keys found: %s", mapkeys)
    good_maps = []
    for key in mapkeys:
      try:
        self.maps[key].map_data.keys()
      except AttributeError:
        # there is no map data yet
        self.logger.info("plot_centered_offsets: getting maps for %s", key)
        if self.maps[key].maps_from_tlogs():
          good_maps.append(key)
        else:
          self.logger.error("plot_centered_offsets: no valid tlog data")
          continue
    fig, ax = subplots(nrows=4, ncols=5)
    fig.set_size_inches(32,30,forward=True)
    width = 0.8; height = 0.8
    for row in range(4):
      for col in range(5):
        key = mapkeys[0] + 5*row + col
        try:
          Map = self.maps[key]
        except KeyError:
          # last map
          break
        if key in good_maps:
          pass
        else:
          continue
        self.maps[key].center_map()
        sca(ax[row][col])
        plot_xdec_dec(self.maps[key].map_data['xdec_offset'],
                                   self.maps[key].map_data['dec_offset'],
                                   self.maps[key].name)
    fig.savefig("raw_positions.png")

  def show_images(self, mapkeys=None):
    """
    show the images for the designated maps
    
    @param mapkeys : numbers of the maps (default: all)
    @type  mapkeys : list of int
    """
    if mapkeys:
      pass
    else:
      mapkeys = self.maps.keys()
    mapkeys.sort()
    first_map = mapkeys[0]
    last_map = mapkeys[-1]
    self.logger.debug("show_images: map keys: %s", mapkeys)
    self.logger.debug("show_images: first map is %d", first_map)
    self.logger.debug("show_images: last map is %d", last_map)
    last_map = mapkeys[-1]
    for key in mapkeys:
      self.logger.debug("show_images: for map %d", key)
      try:
        self.maps[key].map_data.keys()
      except AttributeError:
        if self.maps[key].raster_keys == None:
          self.logger.error("show_images: map has no data")
          continue
        else:
          self.maps[key].maps_from_tlogs()
    
    for mapno in mapkeys:
      try:
        num_chan = len(self.maps[mapno].channels)
      except AttributeError, details:
        # if there is no map data, probably
        self.logger.error("show_images: map %d failed due to attribute err %s",
                          mapno, str(details))
        continue
      if num_chan == 0:
        continue
      if num_chan <= 4:
        fig, axis = subplots(nrows=1, ncols=num_chan)
        ax = [axis] # needs to be a list to be compatible
      else:
        fig, ax = subplots(nrows=2, ncols=num_chan-4)
      fig.set_size_inches(16,10,forward=True)
      width = 0.8; height = 0.8
      self.logger.debug("show_images: processing map %d", mapno)
      Map = self.maps[mapno]
      try:
        channels = Map.channels
      except AttributeError:
        continue
      Map.center_map()
      x,y,z = Map.regrid(width=width, height=height)
      for chan in channels:
        if z.has_key(chan):
          index = channels.index(chan)
          col = index % 4
          row = index//4
          self.logger.debug("show_images: processing map %d channel %d", mapno, chan)
          ax[row][col].contourf(x, y, z[chan], cmap=plt.cm.jet)
          ax[row][col].grid(True)
          ax[row][col].set_aspect('equal')
          ax[row][col].set_xlim(-width/2., width/2.)
          ax[row][col].set_ylim(-height/2., height/2.)
          ax[row][col].set_xlabel("Cross-declination offset")
          ax[row][col].set_ylabel("Declination offset")
          ax[row][col].set_title(
                      Map.name+" at " + str(Map.rss_cfg[chan]['sky_freq'][0]) +
                               " MHz " + Map.rss_cfg[chan]['pol'][0].upper())
          self.logger.debug("show_images: finished row %d, column %d for map %d ch %d",
                          row, col, mapno, chan)
        else:
          continue
      #fig.suptitle(str(int(Map.cfg['year']))+"/"+str(int(Map.cfg['doy'])))
      fig.suptitle(str(self.year)+"/"+str(self.doy))
      fig.savefig(self.session_dir+"map"+ ("%04d" % mapno) +".png")
    show() 
  
  def show_boresights(self, channel=None):
    """
    """
    try:
      bs_keys = self.bs_data.keys()
    except AttributeError:
      self.get_boresight_data(chan=channel)
      bs_keys = self.bs_data.keys()
    bs_keys.sort()
    for key in bs_keys:
      axis = self.xpwr_metadata[:,4][where(self.xpwr_metadata[:,0]==str(key))][0]
      src_id = self.xpwr_metadata[:,3][where(self.xpwr_metadata[:,0]==str(key))][0]
      src_name = self.db.get_source_names([src_id])['source'][0]
      fig = figure()
      if channel:
        if type(channel) == int:
          channels = [channel]
      else:
        channels = unique(self.bs_data[key][:,6])
      for ch in channels:
        UNIXtime = self.bs_data[key][:,0][where(self.bs_data[key][:,6]==ch)]
        Top      = self.bs_data[key][:,1][where(self.bs_data[key][:,6]==ch)]
        plot(UNIXtime, Top, label=str(ch))
      title(src_name+' '+axis+' ('+str(key)+')')
      grid()
      legend()
      fig.savefig(self.session_dir+"boresight"+ ("%05d" % key) +".png")

            
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
    plots the xdec and dec positions at which data were taken
    
    These data come from table 'raster'. These should agree with the grid
    points computed for the channel selected for the Z-plot realtime display.
    For the other channels, the data taken from table 'tlogs' may differ
    because of timing differences in the voltage-to-frequency converters.
    The table records 'az' and 'el' for each 'chan' which can be used to
    compute 'RA' and 'dec' and then 'center_map' can generate relative offsets.
    """
    UNIXtime, xdecs, decs, tsrc = self.get_raster_data()
    plot_xdec_dec(self.raster_data['xdec'],
                  self.raster_data['dec'], self.name)
  
  def plot_centered_offsets(self):
    """
    shows relative positions computed from tlog az and el data
    """
    plot_xdec_dec(self.map_data['xdec_offset'], self.map_data['dec_offset'],
                  self.name)
    
  def contours(self, channel, width=1.0, height=1.0, contours=None):
    """
    make contour maps; this calls 'regrid'
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
    
    
    
