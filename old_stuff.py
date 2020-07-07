  def oldregrid(self, width=1.0, height=1.0, step=None):
    """
    converts a map from observed coordinates to map coordinates
    
    @param width : map width in deg
    @type  width : float
    
    @param height : map height in deg
    @type  height : float
    
    @param step : map step size in deg
    @type  step : float
    """
    if step == None:
      # use the original step size
      step = self.cfg['step']
    nx = int(round(width/step))
    ny = int(round(height/step))
    self.data['grid_x'] = NP.arange( -width/2,  width/2 + step, step/2)
    self.data['grid_y'] = NP.arange(-height/2, height/2 + step, step/2)
    self.data['grid_z'] = {}
    for channel in self.channels:
      cx = self.data['xdec_offset']
      cy = self.data['dec_offset']
      cz = self.data['VFC_counts'][channel]
      xi = self.data['grid_x']
      yi = self.data['grid_y']
      try:
        self.data['grid_z'][channel] = interp.griddata(cx,cy,cz, xi, yi, method='nearest')
      except ValueError as details:
        self.logger.error("regrid: gridding failed: %s", str(details))
        self.logger.debug("regrid: channel %d length of cx is %d", channel, len(cx))
        self.logger.debug("regrid: channel %d length of cy is %d", channel, len(cy))
        self.logger.debug("regrid: channel %d length of cz is %d", channel, len(cz))
        continue
    return self.data['grid_x'], self.data['grid_y'], \
           self.data['grid_z']

    
  def oldcenter_map(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates a map in coordinates relative to a source

    @param body : source at map center
    @type  body : ephem source instance

    @return: (dxdecs,ddecs) in degrees
    """
    try:
      list(self.data.keys())
    except AttributeError:
      self.maps_from_tlogs()
    if 'MPL_datenum' in self.data:
      pass
    else:
      # 'tlog' did not have good data
      return None
    if source.lower() == "sun":
      src = ephem.Sun()
    else:
      pass
    self.data['dec_offset'] = []
    self.data['xdec_offset'] = []
    for count in range(len(self.data['MPL_datenum'])):
      dt = MPL.num2date(self.data['MPL_datenum'][count])
      if source.lower() == "sun":
        src.compute(dt)
      else:
        pass
      ra_center = src.ra*12/math.pi    # hours
      dec_center = src.dec*180/math.pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['RA'][count] - ra_center)*15*math.cos(decrad) )
      self.data['dec_offset'].append(  dec_ofst + 
                        self.data['declination'][count] - dec_center)
    return self.data['xdec_offset'], self.data['dec_offset']

  def old_get_offsets(self, source="Sun", xdec_ofst=0., dec_ofst=0.):
    """
    Generates coordinates relative to a source
    
    If the source is the default, the position of the Sun will be computed for
    the time of each sample. IT SEEMS LIKE A GOOD IDEA TO DO THIS FOR PLANETS
    ALSO.
    
    This adds elements with keys 'xdec_offset' and 'dec_offset' to the
    attribute 'data'.  If this attribute does not exist then
    'get_data_from_tlogs' is called first

    @param source : source at map center
    @type  source : ephem source instance

    @param xdec_ofst : relative X-dec position of sample
    @type  xdec_ofst : float

    @param dec_ofst : relative dec position of sample
    @type  dec_ofst : float
    
    @return: (dxdecs,ddecs) in degrees
    """
    try:
      list(self.data.keys())
    except AttributeError:
      self.get_data_from_tlogs()
    if 'MPL_datenum' in self.data:
      pass
    else:
      # 'tlog' did not have good data
      return None
    if source.lower() == "sun":
      src = ephem.Sun()
    else:
      src = self.calibrator
    self.data['dec_offset'] = []
    self.data['xdec_offset'] = []
    for count in range(len(self.data['MPL_datenum'])):
      dt = MPL.num2date(self.data['MPL_datenum'][count])
      if type(src) == AE.Quasar:
        pass
      else:
        src.compute(dt)
      ra_center = src.ra*12/math.pi    # hours
      dec_center = src.dec*180/math.pi # degrees
      decrad = src.dec
      # right ascension increases to the left, cross-dec to the right
      self.data['xdec_offset'].append(xdec_ofst - 
                       (self.data['RA'][count] - ra_center)*15*math.cos(decrad) )
      self.data['dec_offset'].append(  dec_ofst + 
                        self.data['declination'][count] - dec_center)
    return self.data['xdec_offset'], self.data['dec_offset']
             
class OldMapPlotter(DR.Map):
  """
  """
  def __init__(self, parent=None, name=None, dss=None, date=None, project=None):
    """
    initialize a base Map plotter
    Args:
      parent (Session): optional, an observing session to which this belongs
      name (str):       an identifier, like a scan number; default if not given
      dss (int):        required; station where the data were taken
      date (str):       required; date of observation as "YEAR/DOY"
      project (str):    required; project for which this observation was made
    """
    mylogger = logging.getLogger(logger.name+".MapPlotter")
    DR.Map.__init__(self, parent, name=name, dss=dss, date=date, project=project)
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
      PL.figure()
      PL.plot(self.data['xdec_offset'], self.data['dec_offset'], '-')
      PL.plot(self.data['xdec_offset'], self.data['dec_offset'], '.')
      PL.xlabel("Cross-declination")
      PL.ylabel("Declination")
      PL.title(title_str)
      PL.grid()
      PL.show()
    else:
      self.logger.warning("use self.get_offsets()")      

  def plot_azel(self, title_str="Horizontal Positions"):
    """
    Plot elevation vs azimuth
    """
    if 'az' in self.data and 'el' in self.data:
      PL.figure()
      PL.plot(self.data['az'], self.data['el'], '-')
      PL.plot(self.data['az'], self.data['el'], '.')
      PL.xlabel("Azimuth")
      PL.ylabel("Elevation")
      PL.title(title_str)
      PL.grid()
      PL.show()
    else:
      self.logger.warning("use self.azel_from_radec()")
      
  def plot_ra_dec(self,J2000=False, title_str=None):
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
    PL.figure()
    PL.plot(ras,decs,'-')
    PL.plot(ras,decs,'.')
    PL.xlabel("Right Ascension")
    PL.ylabel("Declination")
    PL.title(title_str)
    PL.grid()
    PL.show()

  def contours(self, chnl, width=None, height=None, contours=None):
    """
    Make contour maps; this calls 'regrid' if necessary
    
    Args:
      chnl (str):                    required signal channel name
      width (float):                 optional, defaults to data cross-dec range
      height (float):                optional, defaults to data dec range
      contours (int or list of int): number of contours or contour levels
    """
    if width and height:
      pass
    else:
      width = (self.data['xdec_offset'].max()-self.data['xdec_offset'].min()).round(2)
      height = (self.data['dec_offset'].max()-self.data['dec_offset'].min()).round(2)
    if 'grid_x' in self.data and 'grid_z' in self.data and  'grid_z' in self.data:
      x,y,z = self.data['grid_x'], self.data['grid_y'], self.data['grid_z']
    else:
      x,y,z = self.regrid(width=width, height=height)
    fig, ax = PL.subplots(nrows=1, ncols=1)
    if contours:
      PL.contour( x, y, z[chnl], contours, linewidths=0.5, colors='k')
      PL.contourf(x, y, z[chnl], contours, cmap=PL.cm.jet)
    else:
      PL.contourf(x, y, z[chnl], cmap=PL.cm.jet)
    PL.grid()
    ax.set_aspect('equal')
    ax.set_xlim(-width/2., width/2.)
    ax.set_ylim(-height/2., height/2.)
    ax.set_xlabel("Cross-declination offset")
    ax.set_ylabel("Declination offset")
    title = self.name
    if self.channel[chnl]['freq']:
      title += " at " + str(self.channel[chnl]['freq']) + " MHz"
    if self.channel[chnl]['pol']:
      title += "  "+ self.channel[chnl]['pol'][0].upper()
    ax.set_title(title)
    PL.colorbar(shrink=0.6) # draw colorbar
    #hh = self.cfg['utc'].seconds/3600
    #mm = (self.cfg['utc'].seconds-3600*hh)/60
    timestruct = time.gmtime(self.start)
    yr = timestruct.tm_year
    dy = timestruct.tm_yday
    hh = timestruct.tm_hour
    mm = timestruct.tm_min
    fig.suptitle("%4d/%03d %02d:%02d UT" % (yr,dy,hh,mm))
    fig.show()
