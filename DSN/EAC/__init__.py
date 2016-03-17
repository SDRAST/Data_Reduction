"""Support for obtaining EAC data from the antenna database"""

import Mysql as sql

def check_for_pointing_cmd(db,dss,UnixTime):
  # Is there a pointing_cmd record for this time?
  query = "SELECT id FROM dss" + str(dss) + \
          ".pointing_cmd WHERE UnixTime = " + str(UnixTime) + ";"
  result = sql.ask_db(db,query)
  return result
