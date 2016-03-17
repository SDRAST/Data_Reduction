# -*- coding: utf-8 -*-
"""
Support for GAVRT mysql database

Uses local module Mysql. Example::
 

The databases and their schemas are described in
http://gsc.lewiscenter.org/data_info/dss28_spec.php
The server has these databases::
 'dss28_eac_v2'
 'dss28_spec'
 'gavrt_sources'.

Database 'dss28_eac_v2' has these tables::
 'five_point',
 'raster',
 'raster_cfg',
 'rss_cfg',
 'seti_cfg',
 'seti_frame',
 'weather',
 'xpwr',
 'xpwr_cfg',
 'xscan'

Database 'gavrt_sources' has these tables::
 'catalog',
 'class',
 'source'

'rss_cfg' columns::
  rss_cfg_id,
  year, doy, utc, chan, sky_freq, feed, pol, nd, if_mode, if_bw, bb_bw, fiber_chan
'weather' columns::
  weather_id,
  datetime, pressure, temp, humidity, wind_speed, wind_dir
'xpwr' columns::
  xpwr_id, xpwr_cfg_id,
  year, doy, utc, epoch, tsys, az, el, ha, dec, offset
'xpwr_cfg' columns::
  xpwr_cfg_id, source_id,
  axis, chan
'xscan' columns::
  xscan_id, xpwr_cfg_id,
  year, doy, utc, epoch, tsrc, stdev, bl_stdev, az, az_offset, el, el_offset,
  ha, dec, offset, bw, corr
'source' columns::
  source_id, catalog_id, class_id,
  name, RA, Dec, size_dec, size_xdec, reference, aka
"""
import MySQLdb
import Mysql
import os
import stat
import pickle
import numpy as np
import logging

mylogger = logging.getLogger("__main__."+__name__)

_host,_user,_pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p", "rb" ))

def _validate_login(host=None, user=None, pw=None):
  """
  """
  if host == None:
    host = "kraken.gavrt.org"
  if user == None:
    user = "ops"
  if pw == None:
    pw = 'v3H!cle'
  return host, user, pw
  
def create_GAVRT_mysql_auth(host=None, user=None, pw=None):
  """
  Create a fairly private authentication file for GAVRT mysql db

  The file goes into your home directory.  It will authenticate you
  automatically.

  @param host : database platform
  @type  host : str

  @param user : authorized user
  @type  user : str

  @param pw : password
  @type  pw : str

  @return: True if successful
  """
  homepath = os.environ['HOME']
  host,user,pw = _validate_login(host,user,pw)
  if host and user and pw:
    GAVRTdb_login = (host,user,pw)
    fd = open(os.environ['HOME']+"/.GAVRTlogin.p", "wb")
    pickle.dump(GAVRTdb_login, fd)
    os.fchmod(fd, int(stat.S_IREAD | stat.S_IWRITE))
    fd.close()
    return True
  else:
    return False

def open_db(database):
  """
  Establishes a connection to the DSS_28 database.

  Returns a connection object and a cursor.  The connection object will be
  used to close the connection.  The cursor object is used to make queries.

  @param database : name of the database (see module doc)
  @type  database : str

  @return: connection object, cursor object
  """
  host, user, pw = _validate_login()
  try:
    conn = MySQLdb.connect(host = '163.150.129.236',
                           port = 3306,
                           user = 'ops',
                           passwd = 'v3H!cle',
                           db = database,
                           compress = True)
    return conn, conn.cursor()
  except AttributeError,info:
    print "Connection failed because",info
    return None,None
  except:
    return None,None

def send_query(crsr,query):
  """
  Send a MySQL query

  @param crsr : cursor object for a database

  @param query : MySQL query
  @type  query : str

  @return: str with query response
  """
  query += ";"
  crsr.execute(query)
  response = crsr.fetchall()
  return response

def extract_boresight_data(year,doy):
  """
  Get the metadata for the boresights on the designated day.

  The boresights are extracted from table 'xscan'.  Missing 'el' data are
  obtained from table 'xpwr'.  The source, scan axis and channel are obtained
  from table 'xpwr_cfg'.  The receiver data are obtained from table 'rss_cfg'.

  Returns a dictionary like this:
  {'utc':        list of datetime.timedelta,
   'epoch':      list of float,
   'az':         list of float,
   'el':         list of value,
   'chan':       list of int,
   'tsrc':       list of float,
   'axis':       list of str,
   'source':     list of str,
   'xpwr_cfg_id: list of int',
   'xscan_id':   list of int,
   'source_id':  list of int,
   'rx':         list of dict}
  An 'rx' dict looks like this:
  { 2: {'if_bw':    float,
        'if_mode':  str,
        'pol':      str,
        'sky_freq': float,
        'utc':      datetime.timedelta},
    4: { ... },
    ....
    16: { ... }}

  @param year : year of observation
  @type  year : int

  @param doy : day of year
  @type  doy : int

  @return: dict
  """
  db = Mysql.BaseDB(_host, _user, _pw, 'dss28_eac_v2')

  # Get the boresight data
  columns = "utc, epoch, tsrc, az, el, xscan_id, xpwr_cfg_id"
  boresight_data = db.get_rows_by_date("xscan",columns,year,doy)

  # Get the missing elevation data
  times = boresight_data['utc']
  power_data = db.get_rows_by_time('xpwr',['utc','el','tsys'],
                                   year,doy,times)
  # Fix the missing elevation data
  boresight_data['el'] = power_data['el']

  # Get the source information
  columns = "source_id, axis, chan"
  for column in columns.split(','):
    boresight_data[column.strip()] = []
  for cfg_id in boresight_data['xpwr_cfg_id']:
    response = db.get("select "
                      + columns
                      + " from xpwr_cfg where xpwr_cfg_id="+str(cfg_id)+";")
    for key in response.keys():
      boresight_data[key].append(response[key][0])
  boresight_data['source'] = []
  for source_id in boresight_data['source_id']:
    response = db.get("select name from gavrt_sources.source where source_id="
                      +str(source_id)+";")
    boresight_data['source'].append(response['name'][0])

  # Get the receiver information
  columns = "utc,sky_freq,pol,if_mode,if_bw"
  boresight_data['rx'] = []
  for time in times:
    boresight_data['rx'].append(get_receiver_data(db,year,doy,time,columns))

  db.c.close()
  return boresight_data

def get_receiver_data(db,year,doy,time,columns):
  """
  Get the receiver state at a given time

  This creates a dictionary keyed with channel number and returns a dictionary
  of the receiver configuration, keyed with specified in the columns, that was
  in effect at the given time.

  Notes
  =====
  Logic
  -----
  The challenge here is to get the latest configuration data for each channel
  at or prior to the specified time.  That channel may have been configured on
  the same day or a prior day. The method we'll use is to find the ID of last
  configuration change and assume that the IDs are sequential in date/time.

  @param db : database
  @type  db : Mysql.BaseDB instance

  @param year : year of observation
  @type  year : int

  @param doy : day of year
  @type  doy : int

  @param time : UTC for the requested receiver state
  @type  time : datetime.timedelta

  @param columns : data items to be returned
  @type  columns : list of str

  @return: dict
  """
  columns = columns.replace(" ","")
  column_keys = columns.split(',')
  receiver_state = {}
  latest_data = db.get("select rss_cfg_id,year,doy,utc from rss_cfg"
                       +" where year <= "+str(year)
                       +" and doy <= "+str(doy)
                       +" and utc <= '"+str(time)
                       +"' order by year desc, doy desc, utc desc limit 1;")
  mylogger.debug("get_receiver_data query response: %s",latest_data)
  cfg_ID = latest_data['rss_cfg_id'][0]
  for chan in [2,4,6,8,10,12,14,16]:
    rx_data = db.get("select "+columns
                     +" from rss_cfg where rss_cfg_id <= "+str(cfg_ID)
                     +" and chan = "+str(chan)
                     +" order by rss_cfg_id desc limit 1;")
    receiver_state[chan] = {}
    for key in column_keys:
      receiver_state[chan][key] = rx_data[key][-1]
  return receiver_state
