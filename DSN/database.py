"""
Functions for extracting complicated information from the DSN database.

The methods obtain data which are stored in different tables according to various
criteria.  They may just be a front end for a MySql query.
"""

import sys
import time as T

import Data_Reduction.GAVRT.mysql as sql
from numpy import isnan

diag = False

beam_keys = ['beamaoff','beamxoff','beameoff','beamhoff','beamcoff','beamdoff']

def get_scan_times(db,dss,start_time,end_time):
  """
  Gets the times when antenna pointing was stable and the data valid.
  
  This reads the database 'dssXX' tables 'pointing_status' to find the
  times when the antenna pointing is stable and 'pointing_cmd' to get the
  pointing offsets.

  @param db : a database connection object obtained from Mysql.open_db()
  @type  db : database object
    
  @param dss : Station number
  @type  dss : int
    
  @param start_time : Unix timestamp
  @type  start_time : float
   
  @param end_time : Unix timestamp
  @type  end_time : float

  @return: (dictionary) Selected records from dssXX.pointing_status as a dict.
  """
  c = db.cursor()
  query = "SELECT * FROM dss" + str(dss) + \
    ".pointing_status WHERE UnixTime >= " + str(start_time) + \
    " AND UnixTime <= " + str(end_time) +" ORDER BY UnixTime;"
  if diag:
    print("get_scan_times query:\n",query)
  result = sql.get_as_dict(db,query)
  if diag:
    print(result)
  return result

def get_scans(db,dss,scan_times_dict):
  """
  Processes a scan times dictionary and returns data about each scan.

  By following foreign keys, this expands on the minimal data returned by
  get_scan-times() and creates a dictionary with lists for each key.

  Notes
  =====
  We begin in a state 'onpoint' == False and scan the table
  'pointing_status until status = 'on point' is detected.  This
  is the start of a scan. We use 'pointing_id' from that row to
  go to the table 'pointing', where we get the aximuth, elevation
  and pointing_cmd_id. If pointing_cmd_id has a valid value then
  we get all the data for that row, which includes source and
  offsets. We then set 'onpoint' to True
  
  If state 'onpoint' == True and status = 'moving' is
  encountered we compute the scan duration, average azimuth and
  elevation, and system temperature, and put it all in 'scan_data'.
  We the set 'onpoint' to False.

  @param db : Obtained from Mysql.open_db()
  @type  db : database connection object

  @param dss : Station ID
  @type  dss : int

  @param scan_times_dict : Result from get_scan-times()
  @type  scan_times_dict : dict

  @return: (dict) Keys of the dictionary are::
      - source:    the source observed,
      - scan_time: the scan start time,
      - exposure:  the scan length,
      - beamaoff:  the beam offset in azimuth,
      - beamxoff:  the beam offset in cross-elevation,
      - beameoff:  the beam offset in elevation,
      - beamhoff:  the beam offset in hour angle,
      - beamcoff:  the beam offset in cross-declination,
      - beamdoff:  the beam offset in declination,
      - azimuth:   the antenna azimuth, and
      - elevation: the antenna elevation
  """
  num_scans = len(scan_times_dict['ID'])
  onpoint = False
  scan_data = {}
  scan_data['source'] = []
  scan_data['scan_time'] = []
  scan_data['exposure'] = []
  scan_data['beamaoff'] = []
  scan_data['beamxoff'] = []
  scan_data['beameoff'] = []
  scan_data['beamhoff'] = []
  scan_data['beamcoff'] = []
  scan_data['beamdoff'] = []
  scan_data['azimuth'] = []
  scan_data['elevation'] = []
  scan_data['tsys'] = []
  for i in range(num_scans):
    status = scan_times_dict['status'][i]
    pointing_id = scan_times_dict['pointing_id'][i]
    # The logic here is that if the antenna status is 'on' then the beam
    # offsets are correct. The next status 'moving' is the end of the scan
    if status == 'on point' and onpoint == False:
      # This starts a scan
      start_scan = scan_times_dict['UnixTime'][i]
      # get the data from table 'pointing'
      query = "SELECT azimuth, elevation, pointing_cmd_id FROM dss" + \
          str(dss) + ".pointing WHERE ID = " + str(pointing_id) + ";"
      if diag:
        print("get_scans query (on point):\n",query)
      pointing_dict = sql.get_as_dict(db, query)
      if diag:
        print(pointing_dict)
      az_start = pointing_dict['azimuth'][0]
      el_start = pointing_dict['elevation'][0]
      pointing_cmd_id = pointing_dict['pointing_cmd_id'][0]
      # get the data from the table 'pointing_cmd'
      if not isnan(pointing_cmd_id):
        # NaN can happen near the beginning of the database when
        # no prior pointing command has been found
        query = "SELECT * FROM dss" + str(dss) + \
          ".pointing_cmd WHERE ID = " + str(pointing_cmd_id) + ";"
        if diag:
          print("get_scans_query:\n",query)
        result_dict = sql.get_as_dict(db,query)
        if diag:
          print(result_dict)
        source_name_id = result_dict['source_name_id'][0]
        beamaoff = result_dict['beamaoff'][0] # azimuth offset
        beamxoff = result_dict['beamxoff'][0] # cross-elevation offset
        beameoff = result_dict['beameoff'][0] # elevation offset
        beamhoff = result_dict['beamhoff'][0] # hour angle offset
        beamcoff = result_dict['beamcoff'][0] # cross-declination offset
        beamdoff = result_dict['beamdoff'][0] # declination offset
        # get the source names
        if diag:
          print("Source name ID:",source_name_id)
        query = "SELECT source_name FROM dsn.source_name WHERE ID = " + \
          str(source_name_id) + ";"
        if diag:
          print("Source name query:",query)
        source_name = sql.ask_db(db,query)[0][0]
        if diag:
         print("Result:",source_name)
        onpoint = True
    if status == 'moving' and onpoint == True:
      # This ends a scan
      end_scan = scan_times_dict['UnixTime'][i]
      query = "SELECT azimuth, elevation, pointing_cmd_id FROM dss" + \
          str(dss) + ".pointing WHERE ID = " + str(pointing_id) + ";"
      if diag:
        print("get_scans query (moving):",query)
      pointing_dict = sql.get_as_dict(db, query)
      if diag:
        print(pointing_dict)
      az_end = pointing_dict['azimuth'][0]
      el_end = pointing_dict['elevation'][0]
      query = "SELECT freq,pol,tsys FROM dss" + str(dss) + \
        ".tsys WHERE UnixTime >= " + str(start_scan) + \
        " AND UnixTime <= " + str(end_scan) + ";"
      if diag:
        print("get_scans query:",query)
      tsys_data = sql.get_as_dict(db,query)
      tsys = 0.0
      num_tsys = len(tsys_data['tsys'])
      for i in range(num_tsys):
        tsys += tsys_data['tsys'][i]
      tsys /= num_tsys
      # Now add to the scan_data
      scan_data['source'].append(source_name)
      scan_data['scan_time'].append(start_scan)
      scan_data['exposure'].append(end_scan - start_scan)
      scan_data['beamaoff'].append(beamaoff)
      scan_data['beamxoff'].append(beamxoff)
      scan_data['beameoff'].append(beameoff)
      scan_data['beamhoff'].append(beamhoff)
      scan_data['beamcoff'].append(beamcoff)
      scan_data['beamdoff'].append(beamdoff)
      scan_data['azimuth'].append((az_start+az_end)/2.)
      scan_data['elevation'].append((el_start+el_end)/2.)
      scan_data['tsys'].append(tsys)
      onpoint = False
  return scan_data

def report_scans(scan_data):
  """
  Pretty print the scan data obtained from get_scans()

  @param scan_data : dictionary
  """

  for i in range(len(scan_data['scan_time'])):
    print(T.ctime(scan_data['scan_time'][i]), \
          ("%6s" % scan_data['source'][i]), \
          ("%4d" % scan_data['exposure'][i]), \
          ("%5.1f" % scan_data['azimuth'][i]), \
          ("%5.1f" % scan_data['elevation'][i]), \
          ("%5.1f" % scan_data['tsys'][i]), end=' ')
    for k in beam_keys:
      if scan_data[k][i] != 0.0:
        print(k,'=',scan_data[k][i], end=' ')
    print("")

def test():
  """
  Test this module
  """
  host = 'localhost'
  dss = 13
  dbname = "dss"+str(dss)
  user = 'kuiper'
  passwd = 'dsnra'
  ion()
  if check_database(host,user,passwd,dbname) == False:
    print("Could not access database",dbname)
    sys.exit(0)
  else:
    db = sql.open_db(dbname,host,user,passwd)
    scan_times_dict = get_scan_times(db,1275017199.0,1275023207.0)
    if diag:
      print(len(scan_times_dict),"scans")
      print("Keys:",list(scan_times_dict.keys()))
    scan_data = get_scans(scan_times_dict)
    report_scans(scan_data)

if __name__ == "__main__":
  test()
