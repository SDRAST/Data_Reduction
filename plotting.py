import Data_Reduction as DR
import pylab as PL
            
class MapPlotter(DR.Map):
  """
  """
  def __init__(self, parent, raster_cfg_id, name=None):
    """
    initialize a GAVRT data plotter
    """
    mylogger = logging.getLogger(logger.name+".MapPlotter")
    Map.__init__(self, parent, raster_cfg_id, name=name)
    self.logger = mylogger

  def plot_xdec_dec(self, title_str="Offsets"):
    """
    plots the xdec and dec positions at which data were taken
    
    These data come from table 'raster'. These should agree with the grid
    points computed for the channel selected for the Z-plot realtime display.
    For the other channels, the data taken from table 'tlogs' may differ
    because of timing differences in the voltage-to-frequency converters.
    The table records 'az' and 'el' for each 'chan' which can be used to
    compute 'RA' and 'dec' and then 'center_map' can generate relative offsets.
    """
    if 'xdec_offset' in self.data and 'dec_offset' in self.data:
      PL.plot(self.data['xdec_offset'], self.data['dec_offset'], '-')
      PL.plot(self.data['xdec_offset'], self.data['dec_offset'], '.')
      PL.xlabel("Cross-declination")
      PL.ylabel("Declination")
      PL.title(title_str)
      PL.grid()
    else:
      self.logger.warning("use self.get_offsets()")      

  def plot_azel(self, title_str="Horizontal Positions"):
    """
    Plot elevation vs azimuth
    """
    if 'az' in self.data and 'el' in self.data:
      PL.plot(self.data['az'], self.data['el'], '-')
      PL.plot(self.data['az'], self.data['el'], '.')
      PL.xlabel("Azimuth")
      PL.ylabel("Elevation")
      PL.title(title_str)
      PL.grid()
    else:
      self.logger.warning("use self.azel_from_radec()")
      
  def plot_ra_dec(J2000=False, title_str=None):
    """
    Plot declination vs right ascension
    """
    if J2000:
      if 'ra2000' in self.data and 'dec2000' in self.data:
        ras,decs = self.data['ra2000'],self.data['dec2000']
        title_str = "J2000 Celestial"
      else:
        self.logger.warning("use self.radec2000_from_radec()")
        return
    else:
      if 'ra' in self.data and 'dec' in self.data:
        ras,decs = self.data['ra'],self.data['dec']
        title_str = "Apparent Celestial"
      else:
        self.logger.warning("use self.radec_from_azel() or self.radec_from_radec2000()")
        return
    PL.plot(ras,decs,'-')
    PL.plot(ras,decs,'.')
    PL.xlabel("Right Ascension")
    PL.ylabel("Declination")
    PL.title(title_str)
    PL.grid()

  def contours(self, channel, width=1.0, height=1.0, contours=None):
    """
    make contour maps; this calls 'regrid'
    """
    if 'grid_x' in self.data and 'grid_z' in self.data and  'grid_z' in self.data:
      pass
    else:
      x,y,z = self.regrid(width=width, height=height)
    fig, ax = PL.subplots(nrows=1, ncols=1)
    if contours:
      PL.contour( x, y, z[channel], contours, linewidths=0.5, colors='k')
      PL.contourf(x, y, z[channel], contours, cmap=PL.cm.jet)
    else:
      contourf(x, y, z[channel], cmap=PL.cm.jet)
    PL.grid()
    ax.set_aspect('equal')
    ax.set_xlim(-width/2., width/2.)
    ax.set_ylim(-height/2., height/2.)
    ax.set_PL.xlabel("Cross-declination offset")
    ax.set_ylabel("Declination offset")
    ax.set_title(self.name+" at " + str(self.rss_cfg[channel]['sky_freq']) + 
                 " MHz " + self.rss_cfg[channel]['pol'][0].upper())
    colorbar(shrink=0.6) # draw colorbar
    #hh = self.cfg['utc'].seconds/3600
    #mm = (self.cfg['utc'].seconds-3600*hh)/60
    timestruct = time.gmtime(self.start)
    yr = timestruct.tm_year
    dy = timestruct.tm_yday
    hh = timestruct.tm_hour
    mm = timestruct.tm_min
    fig.suptitle("%4d/%03d %02d:%02d UT" % (yr,dy,hh,mm))
