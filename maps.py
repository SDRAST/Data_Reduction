# -*- coding: utf-8 -*-
from pylab import *
from matplotlib.mlab import griddata
import ephem
from Astronomy.solar import ra_dec_to_xy

radian = 180./pi

def center_data(date_nums, ras, decs, body, site):
  """
  Generates a map in relative coordinates relative to a source

  @param date_nums : datetime of data point
  @type  date_nums : list of mpl datenums

  @param ras : right ascensions (hrs) of data point
  @type  ras : list of floats

  @param decs : declinations (degs) of data points
  @type  decs : list of floats

  @param body : source at map center
  @type  body : ephem source instance

  @param site : observatory location
  @type  site : ephem site instance

  @return: (dxdecs,ddecs) in degrees
  """

  mjds = []
  ddecs = []
  dxdecs = []
  for count in range(len(date_nums)):
    dt = num2date(date_nums[count])
    body.compute(dt)
    ra_center = body.ra*12/pi    # hours
    dec_center = body.dec*180/pi # degrees
    decrad = body.dec
    # right ascension increases to the left, cross-dec to the right
    dxdecs.append(-(ras[count] - ra_center)*15*cos(decrad))
    ddecs.append(decs[count] - dec_center)
  return dxdecs,ddecs

def plot_ra_dec(ras,decs,title_str):
  """
  Plot declination vs right ascension
  """
  plot(ras,decs,'-')
  plot(ras,decs,'.')
  xlabel("Right Ascension")
  ylabel("Declination")
  title(title_str)
  grid()

def plot_xdec_dec(xdecs,decs,title_str):
  """
  plot declination vs cross-declination
  """
  plot(xdecs,decs,'-')
  plot(xdecs,decs,'.')
  xlabel("Cross-declination")
  ylabel("Declination")
  title(title_str)
  grid()

def plot_azel(azs,els,title_str):
  """
  Plot elevation vs azimuth
  """
  plot(azs,els,'-')
  plot(azs,els,'.')
  xlabel("Azimuth")
  ylabel("Elevation")
  title(title_str)
  grid()

def show_body_orientation(body, meantime, map_params, solar_data):
  """
  """
  xdec_off,dec_off,mapwidth,mapheight = map_params
  if body.capitalize() == "Sun":
    obj = ephem.Sun()
  else:
    return
  obj.compute(meantime)
  obj_ra = obj.ra*radian
  obj_dec = obj.dec*radian
  for ddec in arange(-1., 1.+0.1, 0.1):
    x = []
    y = []
    for dra in arange(-1., 1.+0.1, 0.1):
      ra = obj_ra + dra
      dec = obj_dec + ddec
      xy = ra_dec_to_xy(ra, dec, -solar_data["polar angle"],
                        obj_ra, obj_dec)
      x.append(xy[0])
      y.append(xy[1])
      plot(x,y,"y:")
      plot(-array(y),array(x),"y:")
  # draw the circle
  x = []
  y = []
  for theta in linspace(0,360/radian,180):
    x.append(solar_data["semidiameter"]*cos(theta)/60.)
    y.append(solar_data["semidiameter"]*sin(theta)/60.)
  plot(x,y,"y-")
  return

