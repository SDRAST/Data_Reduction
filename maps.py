# -*- coding: utf-8 -*-
"""
Obsolete module replaced by ``plotting`` and the ``MapPlotter`` class.
"""
import Astronomy.solar as solar
import math
import pylab as PL
import ephem
# from Astronomy.solar import ra_dec_to_xy

radian = 180./math.pi

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
    ra_center = body.ra*12/math.pi    # hours
    dec_center = body.dec*180/math.pi # degrees
    decrad = body.dec
    # right ascension increases to the left, cross-dec to the right
    dxdecs.append(-(ras[count] - ra_center)*15*math.cos(decrad))
    ddecs.append(decs[count] - dec_center)
  return dxdecs,ddecs

def plot_ra_dec(ras,decs,title_str):
  """
  Plot declination vs right ascension
  """
  PL.plot(ras,decs,'-')
  PL.plot(ras,decs,'.')
  PL.xlabel("Right Ascension")
  PL.ylabel("Declination")
  PL.title(title_str)
  PL.grid()

def plot_xdec_dec(xdecs,decs,title_str):
  """
  plot declination vs cross-declination
  """
  PL.plot(xdecs,decs,'-')
  PL.plot(xdecs,decs,'.')
  PL.xlabel("Cross-declination")
  PL.ylabel("Declination")
  PL.title(title_str)
  PL.grid()

def plot_azel(azs,els,title_str):
  """
  Plot elevation vs azimuth
  """
  PL.plot(azs,els,'-')
  PL.plot(azs,els,'.')
  PL.xlabel("Azimuth")
  PL.ylabel("Elevation")
  PL.title(title_str)
  PL.grid()

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
      xy = solar.ra_dec_to_xy(ra, dec, -solar_data["polar angle"],
                        obj_ra, obj_dec)
      x.append(xy[0])
      y.append(xy[1])
      PL.plot(x,y,"y:")
      PL.plot(-array(y),array(x),"y:")
  # draw the circle
  x = []
  y = []
  for theta in linspace(0,360/radian,180):
    x.append(solar_data["semidiameter"]*math.cos(theta)/60.)
    y.append(solar_data["semidiameter"]*math.sin(theta)/60.)
  PL.plot(x,y,"y-")
  return

