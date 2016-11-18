"""
Methods to generate data for the tables in the DSN databases.

The order of events in filling antenna data tables should be:
 1. get the antenna commands from the EAC proc_macro log. If this log
  - doesn't exist, the information can be gleaned with difficulty from
  - the antenna status table.  This will put data in DSSxx.POINTING_CMD
  - and if so modifies entries in DSSxx.POINTING.
 2. fill tables DSSxx.FIVE_PT and DSSxx.BORESIGHT
  - This also puts data into the DSSxx.SUBREFLECTOR, DSSxx.POINTING_CMD
  - and DSSxx.POINTING TABLES.
 3. fill table DSSxx.MINICAL
 4. fill the table DSSxx.TSYS
The order is important because primary keys must be generated before foreign
keys can be created to refer to them.

NEEDS MAJOR UPDATES FROM MODULE Observatory TO MonitorControl
"""
from os.path import basename
import Astronomy as A
import DatesTimes as DT
import time as T
import Data_Reduction.DSN.GAVRT.Mysql as sql
from Observatory.mysql import *
from Data_Reduction.DSN.EAC import *
from numpy import *

diag = False

def find_prior_pointing_command(db,dss,UnixTime):
  """Returns the ID for most recent time prior to 'UnixTime' that
  a record exists in the table 'pointing_cmd'."""
  if diag:
    print "find_prior_pointing_command: UnixTime =",UnixTime
  query = "SELECT  MAX(UnixTime) FROM dss" + \
          str(dss) + ".pointing_cmd WHERE UnixTime < " + \
          str(UnixTime) + ";"
  result = sql.ask_db(db,query)
  if diag:
    print "find_prior_pointing_command query:\n",query
    print "find_prior_pointing_command result:",result
  max_time = result[0][0]
  if max_time != None:
    if diag:
      print "find_prior_pointing_command: max_time =",max_time
    query = "SELECT id FROM dss" + str(dss) + \
        ".pointing_cmd WHERE UnixTime = " + str(max_time) + ";"
    result = sql.ask_db(db,query)
    if diag:
      print "find_prior_pointing_command query:\n",query
      print "find_prior_pointing_command result: ", result
    prior_record_id = result[0][0]
    if diag:
      print "find_prior_pointing_command: prior record", \
          prior_record_id,  "at", T.ctime(UnixTime)
    return prior_record_id
  else:
    return -1

def pointing_cmd_event(db,year,timestr,src_id,src_name,src_name_id,dss,
  beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff):
  """Generates a dictionary with offset event particulars to be added
  as another row in table pointing_cmd"""
  UnixTime = DT.macro_log_time_to_UnixTime(year,timestr)
  return pointing_cmd_event2(db,UnixTime,src_id,src_name,src_name_id,dss,
  beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff)

def pointing_cmd_event2(db,UnixTime,src_id,src_name,src_name_id,dss,
  beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff):
  """Generates a dictionary with offset event particulars. It is added
  as another row in table pointing_cmd, or if a row already exists for
  that time, it updates the row."""
  pointing_event = {}
  pointing_event['UnixTime'] = UnixTime
  pointing_event['azimuth'], \
  pointing_event['elevation'] = A.mysql.get_pointing_data(db,
                                                          UnixTime,
                                                          src_name,
                                                          dss)
  offset_event = {}
  offset_event['UnixTime'] = UnixTime
  offset_event['source_name_id'] = src_name_id
  offset_event['beamaoff'] = beamaoff
  offset_event['beamxoff'] = beamxoff
  offset_event['beameoff'] = beameoff
  offset_event['beamhoff'] = beamhoff
  offset_event['beamcoff'] = beamcoff
  offset_event['beamdoff'] = beamdoff
  # Try to insert the new event
  pointing_cmd_id,msg = sql.insert_record(db,'dss13.pointing_cmd',offset_event)
  if pointing_cmd_id == -2:
    # A record exists for this time; update it
    print "pointing_cmd_event2: Duplicate pointing command for", \
      T.ctime(offset_event['UnixTime'])
    query = "SELECT ID FROM dss"+str(dss)+".pointing_cmd WHERE UnixTime = "+ \
      str(UnixTime)+";"
    result = sql.ask_db(db,query)
    pointing_cmd_id = result[0][0]
    condition = ("ID",pointing_cmd_id)
    sql.update_record(db,"dss"+str(dss)+".pointing_cmd",offset_event,condition)
  elif pointing_cmd_id < 0:
    print "pointing_cmd_event2: Pointing command for", \
      T.ctime(offset_event['UnixTime']), \
      "could not be stored. Error",pointing_cmd_id
    print msg
  if pointing_cmd_id > 0:
    # update the foreign key in the associated record in 'pointing' 
    pointing_event['pointing_cmd_id'] = pointing_cmd_id
    pointing_id,msg = sql.insert_record(db,'dss13.pointing',pointing_event)
    if pointing_id == -2:
      print "pointing_cmd_event2: Duplicate pointing data for",T.ctime(UnixTime)
      msg_parts = msg.split("'")
      print "conflicts with record for",T.ctime(float(msg_parts[1]))
      query = "SELECT ID FROM dss"+str(dss)+".pointing WHERE UnixTime = "+ \
        str(UnixTime)+";"
      result = sql.ask_db(db,query)
      pointing_id = result[0][0]
      condition = ("ID",pointing_id)
      sql.update_record(db,"dss"+str(dss)+".pointing",pointing_event,condition)
    elif pointing_id < 0:
      print "pointing_cmd_event2: Pointing data for",T.ctime(UnixTime), \
        "could not be stored. Error",pointing_id
      print msg
  return pointing_id,pointing_cmd_id

def get_antenna_status(db,dss,fname):
  """This parses a pfYYYDDD.txt file which also has antenna status
  information.  The file is in the form of a space-separated table and
  has these columns::
  SOURCE - source name
  STATUS - see below
  DOY    - UTC day of year
  UTC    - HH:MM:SS
  AZ     - decimal degrees
  EL     - decimal degrees
  The status column has these possible values::
  cmd_off      - position offset (including clear offset) command issued;
                 implicitly, the antenna is not on point
  subr_off     - subreflector error; ignore this
  tracking_off - not tracking and not in maintenance mode
  on           - antenna is at commanded position (including offset)
  off          - antenna is off the commanded position
  A 'cmd_off' is usually followed by an 'on' about 10(!) seconds later.
  'off's and 'on's may also happen in pairs, separated by a few seconds, if
  the antenna is momentarily outside the tolerance circle.

  In trying to parse this information for implicit antenna commands, bear
  in mind that when a 'cmd_off' occurs, the antenna is at the position
  where it had previously been commanded to be.  The subsequent 'on' shows
  that the new position has been achieved.  Other 'off' and 'on' pairs can
  be ignored.

  Most of the 10 sec delay is probably due to settling at the commanded
  position.  We can probably assume that the antenna is actually at the
  commanded position TBD seconds before the 'on' status is reported.  The
  time between 'cmd_off' and 'on'-TBD will be considered invalid data time.

  The antenna pointing status status file has time, azimuth and elevation
  so that data can be entered into the table 'pointing'.  However, no
  pointing commands have yet been processed so that column will have only
  nulls.
  """
  # parse file name for year
  year = int(basename(fname)[2:6])
  if diag:
    print "get_antenna_status: year =",year
  cmd_file = open(fname,"r")
  lines = cmd_file.readlines()
  cmd_file.close()
  data = [] # a list of all the events contained in dictionaries
  keys = [] # the table keys (or column headers)
  # parse the file into list 'data' of dictionaries for each line
  print "get_antenna_status:",len(lines),"read"
  for l in lines:
    linedict = {}
    line = l.rstrip().split(' ')
    if line[0] == 'SOURCE':
      # This is a line with column headers; save them as keys
      if keys == []:
        keys = line
    else:
      # parse the line data into the dictionary 'linedict'
      for k in range(len(keys)):
        if keys[k] == 'SOURCE' or keys[k] == 'STATUS':
          # string
          linedict[keys[k]] = line[k]
        elif keys[k] == 'DOY':
          DOY = int(line[k])
        elif keys[k] == 'UTC':
          start_sec = DT.HHMMSS_to_seconds(line[k])
        else:
          linedict[keys[k]] = float(line[k])
      UnixTime = DT.VSR_tuple_to_timestamp(year, DOY, start_sec)
      linedict['UnixTime'] = UnixTime # not needed
      # Now the data on the line have been converted
      # prepare the data for insertion into the tables
      pointing_data = {}
      status_data = {}
      pointing_data['UnixTime'] = UnixTime
      pointing_data['azimuth'] = linedict['AZ']
      pointing_data['elevation'] = linedict['EL']
      prior_pointing_command_id = find_prior_pointing_command(db,dss,UnixTime)
      if prior_pointing_command_id > 0:
        pointing_data['pointing_cmd_id'] = prior_pointing_command_id
      # At this point we don't know if there is an associated pointing
      # command.  Try to insert a record.
      if diag:
        print "get_antenna_status: inserting\n",pointing_data
      pointing_id,msg = sql.insert_record(db,'dss'+str(dss)+".pointing",
                                     pointing_data)
      if pointing_id == -2:
        # duplicate record
        if diag:
          parts = msg.split("'")
          UxTm = float(parts[1])
          print "Duplicate record for",T.ctime(UxTm)
        # what is the existing record's ID?
        result = sql.ask_db(db,
          "SELECT ID FROM dss"+str(dss)+".pointing WHERE UnixTime = " +
          str(UnixTime) + ";")
        if result == None or result == ():
          print 'Pointing data for',T.ctime(UnixTime),'is not in the database'
        else:
          pointing_id = result[0][0]
      elif pointing_id < 0:
        # other error
        print "get_antenna_status: pointing table query for",T.ctime(UnixTime), \
            "returned",pointing_id
        print msg
      if pointing_id > 0:
        # good insertion
        status_data['UnixTime'] = UnixTime
        status_data['pointing_id'] = pointing_id
        if linedict['STATUS'] == 'on':
          status_data['status'] = "on point"
        elif linedict['STATUS'] == "cmd_off":
          status_data['status'] = "moving"
        else:
          status_data['status'] = "off point"
        (src_name_id,src_id) = A.mysql.get_source_IDs(db,linedict['SOURCE'])
        status_data['source_id'] = src_id
        status_id,msg = sql.insert_record(db,'dss'+str(dss)+".pointing_status",
                                 status_data)
      data.append(linedict)
  # Now figure out what each line meant
  print len(data),"events found"
  return data
  
def get_EAC_offsets(db,fname,dss):
  """Parses an EAC proc_macro log for commands which change the antenna
  pointing, including new sources and changes in manual pointing offsets
  (not boresight offsets).

  The EAC macro log is a narrative log whose form evolved during the first
  half of 2010, mainly by adding date-time codes at the beginning of most
  (but not yet all) lines.

  Figuring out what commands were issued when gets to be pretty tricky if
  the macro log was not kept.  Be sure to pipe the macro processor output to
  stdout (console) into a log file.

  The required information for the table pointing_cmd is the time, the
  name of the source being tracked, and the six offsets.

  This function assumes no offsets at the start, although it would be nice
  to have a 'CLR PO ALL' early in the macro log to confirm this. When such
  a command is detected, it is recorded as an actual 'pointing_cmd' event at
  the next 'clr po ... completed'.

  For 'CLR PO {something}' the appropriate offsets in table 'pointing_cmd'
  are set to zero.

  A 'poffset' command contains the offset values which are used to set
  the appropriate parameters.  The offsets are stored as a 'pointing_cmd'
  event at the next 'poffset: completed' log entry.

  Command 'source' causes this function to remember the source name, its
  key ID in the table 'source_name' and in table 'source'.  All subsequent
  'pointing_cmd' entries reference that source.  A pointing_cmd event is
  recorded at the subsequent 'onpoint: Tracking'.  The offset values for that
  event are whatever was last determined.

  During a boresight, indicated by a 'pscan xdecdec' or something like
  that followed by a 'pscan: completed', the tables 'pointing_status' and
  'pointing' are processed for the boresight details. Every 'cmd_off' is
  the start of a move and the subsequent 'on' is its completion.  There
  are ten moves, five in elevation and five in cross-elevation.  They are::
    - baseline-offset
    - hpbw/2 offset
    0 offset (center)
    + hpbw/2 offset
    + baseline offset
  The next 'on' after the last offset is assumed to be the boresight
  position.

  When all the commands have been processed, a second pass is needed over
  table 'pointing' to replace all the NULL foreign key references to
  table 'pointing_cmd' with the previous value in time.
  """
  Long,Lat,elev,x,y,x = get_obs_coords(db,dss)
  offset_changed = False
  beamaoff = 0.0
  beamxoff = 0.0
  beameoff = 0.0
  beamhoff = 0.0
  beamcoff = 0.0
  beamdoff = 0.0
  source_name_id = 2 # defaults to no name
  src_id = 2
  in_boresight = False
  fd = open(fname,'r')
  lines = fd.readlines()
  # get the year from the filename
  year = int(basename(fname).split('-')[2])
  for l in lines:
    offset_event = {}
    line = l.strip().split()
    if line[1] == 'clr' and line[2] == 'PO':
      # clear offsets command; these values do not take effect until
      # 'clr po completed' has been received
      timestr = line[0]
      offset_changed = True
      if line[3] == 'ALL':
        beamaoff = 0.0
        beamxoff = 0.0
        beameoff = 0.0
        beamhoff = 0.0
        beamcoff = 0.0
        beamdoff = 0.0
      elif line[3] == 'AZ':
        beamaoff = 0.0
      elif line[3] == 'EL':
        beameoff = 0.0
      elif line[3] == 'XEL':
        beamxoff = 0.0
      elif line[3] == 'HA':
        beamhoff = 0.0
      elif line[3] == 'DEC':
        beamdoff = 0.0
      elif line[3] == 'XDEC':
        beamcoff = 0.0
    if line[1] == 'source':
      # New source command
      # We need the source to get the coordinates
      # lines with 'source:' (note colon) are ignored
      src_name = line[2]
      (src_name_id,src_id) = A.mysql.get_source_IDs(db,src_name)
    elif line[1] == 'onpoint:' and line[2] == 'Tracking' and src_id > 2:
      # It's hard to imagine a 'Tracking' status without a source having
      # been specified, but we check just in case
      pointing_id, pointing_cmd_id  = \
        pointing_cmd_event(db,year,line[0],src_id,src_name,src_name_id,
        dss,beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff)
    elif line[1] == 'poffset':
      # position offset command; not an event until it is completed
      if line[2] == 'AZ':
        beamaoff = float(line[3])
      elif line[2] == 'XEL':
        beamxoff = float(line[3])
      elif line[2] == 'EL':
        beameoff = float(line[3])
      elif  line[2] == 'HA':
        beamhoff = float(line[3])
      elif  line[2] == 'XDEC':
        beamcoff = float(line[3])
      elif  line[2] == 'DEC':
        beamdoff = float(line[3])
      offset_changed = True
    elif line[1] == 'clr' and line[2] == 'po':
      if len(line) > 3 and line[4] == 'completed':
        # Confirmation of 'clear offset'; have new offset
        if src_id > 2 and offset_changed:
          timestr = line[0]
          event = pointing_cmd_event(db,year,line[0],src_id,src_name,
            src_name_id,dss,beamaoff,beamxoff,beameoff,beamhoff,beamcoff,
            beamdoff)
          offset_changed = False
    elif line[1] == 'poffset:' and line[2] == 'completed':
      # confirmation of offset command completed
        timestr = line[0]
        if src_id > 2 and offset_changed:
          event = pointing_cmd_event(db,year,line[0],src_id,src_name,
            src_name_id,dss,beamaoff,beamxoff,beameoff,beamhoff,beamcoff,
            beamdoff)
          offset_changed = False
    elif line[1] == 'pscan' and line[2] != 'completed':
      # This is the beginning of a boresight
      in_boresight = True
      boresight_start = DT.macro_log_time_to_UnixTime(year,line[0])
      boresight_mode = line[2]
    elif line[1] == 'pscan:' and line[2] == 'completed':
      # end of the boresight
      boresight_end = DT.macro_log_time_to_UnixTime(year,line[0])
      # Now we can try to make sense of what is in the pointing and
      # status tables.  Grab all the pointing data for this boresight
      query = "SELECT UnixTime,status,pointing_id FROM dss" + str(dss) + \
        ".pointing_status WHERE UnixTime >= "+ \
        str(boresight_start)+" AND UnixTime <= "+str(boresight_end)+ \
        " ORDER BY UnixTime;"
      if diag:
        print "get_EAC_offsets query:",query
      result = sql.ask_db(db,query)
      count = 0
      ax_offs = []
      el_offs = []
      bf_time = []
      looking_for_onpoint = False
      for item in result:
        UnixTime, status, pointing_id = item
        query = "SELECT azimuth,elevation FROM dss"+str(dss)+ \
          ".pointing WHERE dss"+str(dss)+".pointing.UnixTime = "+ \
          str(UnixTime)+";"
        result = sql.ask_db(db,query)
        if diag:
          print "get_EAC_offsets:", T.ctime(UnixTime),status, \
            ("%8.4f" % result[0][0]), ("%8.4f" % result[0][1]), \
            looking_for_onpoint
        # You can compare floats but the MySQL data type Decimal, converted
        # to string, will match if equal.
        if str(UnixTime) != ("%12.1f" %boresight_start):
          if status == "moving":
            looking_for_onpoint = True
          if looking_for_onpoint and status == "on point":
            Az,El = result[0]
            # empirically, if there is already a record in pointing_cmd then
            # the next item is the start of the boresight, so skip it
            if check_for_pointing_cmd(db,dss,UnixTime) == ():
              # OK. No pointing command for this pointing record
              # calculate az and el
              az, el = A.mysql.calculated_AzEl(db,src_name,UnixTime,Lat)
              total_ax_offset = (Az - az)*cos(El*pi/180.)
              ax_offs.append(total_ax_offset)
              total_el_offset = El - el
              el_offs.append(total_el_offset)
              bf_time.append(float(UnixTime)) # UnixTime is type 'Decimal'
              count += 1
              print "get_EAC_offsets:", count,T.ctime(UnixTime), \
                total_ax_offset,total_el_offset,status
              looking_for_onpoint = False
            else:
              # This should not be possible for an initially empty table
              print "get_EAC_offsets: Pointing cmd exists for", \
                T.ctime(UnixTime),UnixTime
        else:
          print str(UnixTime),"=",("%12.1f" %boresight_start)
      if count < 10:
        # get the next on point
        query = "SELECT min(UnixTime) FROM dss" + str(dss)+ \
          ".pointing_status WHERE UnixTime > "+str(boresight_end) + \
          " AND status = 'on point';"
        result = sql.ask_db(db,query)
        onpoint_time = result[0][0]
        query = "SELECT azimuth, elevation FROM dss" + str(dss)+ \
          ".pointing WHERE UnixTime = " + str(onpoint_time) + ";"
        result = sql.ask_db(db,query)
        if diag:
          print "get_EAC_offsets query:", query
          print "get_EAC_offsets result:", result
        Az,El = result[0]
        az, el = A.mysql.calculated_AzEl(db,src_name,onpoint_time,Lat)
        total_ax_offset = (Az - az)*cos(El*pi/180.)
        ax_offs.append(total_ax_offset)
        total_el_offset = El - el
        el_offs.append(total_el_offset)
        bf_time.append(float(UnixTime)) # UnixTime is type 'Decimal'
        count += 1
      if diag:
        print "get_EAC_offsets:", count,"boresight measurements"
      # We've got all the data.  Fit the non-scan direction positions
      #     elevation data
      Xe_off  = array(ax_offs[0:5])
      Ee_off  = array(el_offs[0:5])
      Te = array(bf_time[0:5])
      coef_Xe_off = polyfit(Te, Xe_off,1)
      linear_X_off = poly1d(coef_Xe_off)
      center_X_off = linear_X_off(Te[2])
      #     azimuth data
      Xa_off  = array(ax_offs[5:10])
      Ea_off  = array(el_offs[5:10])
      Ta = array(bf_time[5:10])
      coef_Ea_off = polyfit(Ta,Ea_off,1)
      linear_E_off = poly1d(coef_Ea_off)
      center_E_off = linear_E_off(Ta[2])
      # Now estimate the commanded boresight offsets and correct the
      # pointing_cmd table.  In the non-scan direction we can remove a linear
      # fit.  In the scan direction we can only remove the mean determined
      # when it was the non-scan axis.
      #   elevation scan
      if boresight_mode == 'xdecdec':
        beamaoff = 0.0
        beamhoff = 0.0
        beamcoff = 0.0
        beamdoff = 0.0
        for i in range(5):
          fields = {}
          # linear correction to cross-elevation
          beamxoff = (Xe_off[i] - linear_X_off(Te[i]))
          # mean correction to elevation
          beameoff = Ee_off[i] - center_E_off
          pointing_id, pointing_cmd_id = \
            pointing_cmd_event2(db,Te[i],src_id,src_name, src_name_id,dss,
            beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff)
          if diag:
            print \
              "get_EAC_offsets: calculated azimuth and elevation stored in", \
              pointing_id
            print "get_EAC_offsets:", beamxoff, beameoff, "stored in", \
              pointing_cmd_id
        for i in range(5):
          # mean correction to cross-elevation
          beamxoff = (Xa_off[i] - center_X_off)
          # linear correction to elevation
          beameoff = Ea_off[i] - linear_E_off(Ta[i])
          pointing_id, pointing_cmd_id = \
            pointing_cmd_event2(db,Ta[i],src_id,src_name, src_name_id,dss,
            beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff)
          if diag:
            print \
              "get_EAC_offsets: calculated azimuth and elevation stored in", \
              pointing_id
            print "get_EAC_offsets:", beamxoff,beameoff, "stored in", \
              pointing_cmd_id
        # The next status event should be an 'on' and it should be at the
        # boresighted position.  Get the onpoint time.
        query = 'select min(UnixTime) from dss' + str(dss) + \
          '.pointing_status where UnixTime > ' + str(Ta[4]) + \
          ' and status = "on point";'
        result = sql.ask_db(db,query)
        onpoint_time = result[0][0]
        # get the onpoint encoder azimuth and elevation
        query = "select * from dss"+str(dss)+".pointing where UnixTime = "+ \
          str(onpoint_time)+";"
        if diag:
          print "get_EAC_offsets query:", query
        result_dict = sql.get_as_dict(db,query)
        if diag:
          print "get_EAC_offsets result:", result_dict
        # get the onpoint calculated azimuth and elevation
        az, el = A.mysql.calculated_AzEl(db,src_name,onpoint_time,Lat)
        # correction to cross-elevation
        delta_az = result_dict['azimuth'] - az
        delta_el = result_dict['elevation'] - el
        beamxoff = (delta_az - center_X_off)*cos(el*pi/180.)
        beameoff = (delta_el - center_E_off)
        pointing_id, pointing_cmd_id = \
            pointing_cmd_event2(db,onpoint_time,src_id,src_name,
            src_name_id,dss,
            beamaoff,beamxoff,beameoff,beamhoff,beamcoff,beamdoff)
      else:
        print "get_EAC_offsets: Boresight mode", boresight_mode, \
          "not yet coded"
  # Now replace all the NULL foreign key values in table 'pointing' with the
  # previous reference.
  # find out where to start
  query = " SELECT min(ID) FROM dss" + str(dss) + \
    ".pointing WHERE pointing_cmd_id IS NOT NULL;"
  start_id = int(sql.ask_db(db,query)[0][0])
  if diag:
    print "Will start processing pointing from ID=",start_id
  query = "SELECT ID, pointing_cmd_id FROM dss" + str(dss) \
    + ".pointing ORDER by UnixTime;"
  result_dict = sql.get_as_dict(db,query)
  if diag:
    print "get_EAC_offsets query:", query
    print "get_EAC_offsets result:", result_dict
  pointing_cmd_id = None
  for i in range(start_id,len(result_dict['ID'])):
    print i, "Found:", \
        result_dict['ID'][i], \
        pointing_cmd_id, \
        result_dict['pointing_cmd_id'][i]
    if pointing_cmd_id == None and result_dict['pointing_cmd_id'][i] > 0:
      # This should cause the first foreign key to be remembered
      pointing_cmd_id = result_dict['pointing_cmd_id'][i]
      if diag:
        print "First foreign key value:",pointing_cmd_id
    elif pointing_cmd_id != None and isnan(result_dict['pointing_cmd_id'][i]):
      if diag:
        print "Updating record",i,"with",pointing_cmd_id
      new_dict = {'pointing_cmd_id': pointing_cmd_id}
      condition = ("ID",result_dict['ID'][i])
      sql.update_record(db,"dss"+str(dss)+".pointing",new_dict,condition)
    elif not isnan(result_dict['pointing_cmd_id'][i]):
      pointing_cmd_id = result_dict['pointing_cmd_id'][i]
      if diag:
        print "Found new foreign key value:",pointing_cmd_id
    else:
      print "get_EAC_offsets query: impossible if case"
    print i, \
        result_dict['ID'][i], \
        pointing_cmd_id, \
        result_dict['pointing_cmd_id'][i]
  return l

def upload_5p_log(database,logpath,dss):
  """This parses an EAC boresight log for boresight results::
    logpath - the full path and log file name, which must be 5pYYYYDDD.pri
    dss -     the DSN station for this log
  The data are stored in these tables in the order shown because foreign key
  values must be determined for some tables::
    subreflector
    pointing_cmd (needed by pointing)
    pointing (needed by five_pt and boresight)
    five_pt (needed by boresight)
    boresight.
  
  The fields in 'pointing' are::
    UnixTime,
    antenna_id,
    source_name_id, azimuth, elevation,
    beamaoff, beamxoff, beameoff, beamhoff, beamcoff, beamdoff.

  The fields in 'subreflector' are::
    UnixTime, ypos, zpos, yoff, zoff.
  
  The fields in 'five_pt' to be provided are::
    pointing_id,
    freq, axhxpos, eldepos,
    deltaXp, deltaX0, deltaXm, baseXp, baseXm, X_bw, X_err, Xpos,
    deltaYp, deltaY0, deltaYm, baseYp, baseYm, Y_bw, Y_err, Ypos,
    subreflector_id.

  The fields in 'boresight' are::
    pointing_id,
    xeloff, eloff, xdecoff, decoff,
    five_pt_id
  """
  # These are the entries on each row in the order in which they appear
  fivept_keys = ['DOY', 'UTC', 'SOURCE', 'FREQ', 'POL', 'HA', 'DEC',
                 'AZ', 'EL', 'AXHXPOS', 'ELDEPOS',
                 'DELTAXP', 'DELTAX0', 'DELTAXM', 'NUMPTSX', 'RMSX',
                 'BASEXP', 'ONX1', 'ONX2', 'BASEXM', 'X_BW', 'X_ERR', 'XPOS',
                 'DELTAYP', 'DELTAY0', 'DELTAYM', 'NUMPTSY', 'RMSY',
                 'BASEYP', 'ONY1', 'ONY2', 'BASEYM', 'Y_BW', 'Y_ERR', 'YPOS',
                 'XELOFF', 'ELOFF', 'XDECOFF', 'DECOFF',
                 'YPOS', 'ZPOS', 'YOFF', 'ZOFF']
  numkeys = len(fivept_keys)
  logobj = open(logpath,'r')
  logdata = logobj.readlines()
  
  # get a cursor for the database
  try:
    cursor = database.cursor()
  except Exception, details:
    print "Could not get a cursor to database"
    print details[1]
    return False
  
  # check that the required tables exist and create data lists
  if sql.check_table(database,"pointing"):
    pointing_data = {}
  else:
    print "Could not access table 'pointing'"
    return False
  
  if sql.check_table(database,"pointing_cmd"):
    pointing_cmd_data = {}
  else:
    print "Could not access table 'pointing_cmd'"
    return False
  
  if sql.check_table(database,"subreflector"):
    subrefl_data = {}
  else:
    print "Could not access table 'subreflector'"
    return False
  
  if sql.check_table(database,"five_pt"):
    five_pt_data = {}
  else:
    print "Could not access table 'five_pt'"
    return False
  
  if sql.check_table(database,"boresight"):
    boredata = {}
  else:
    print "Could not access table 'boresight'"
    return False
  
  # Parse the file name for the year
  year = int(basename(logpath)[2:6])

  # Open the database
  #db = open_db('dss'+str(dss),host,user,passwd)
  
  # Process the rows in the file
  for line in logdata:
    linedata = line.strip().split()
    if linedata[0] == 'DOY':
      # get the keys
      pass
    elif len(linedata) != numkeys:
      #there's something wrong with this line
      print "Oops:",line
    else:
      linedict = {} # dictionary containing the data on the line
      for index in range(numkeys):
        linedict[fivept_keys[index]] = linedata[index]
      # compute the UNIX timestamp
      start_sec = DT.HHMMSS_to_seconds(linedict['UTC'])
      UnixTime = DT.VSR_tuple_to_timestamp(year,
                                           int(linedict['DOY']),
                                           start_sec)

      # First do the 'subreflector' table and remember the key
      subrefl_data['UnixTime'] = UnixTime
      subrefl_data['ypos'] = linedict['YPOS']
      subrefl_data['zpos'] = linedict['ZPOS']
      subrefl_data['yoff'] = linedict['YOFF']
      subrefl_data['zoff'] = linedict['ZOFF']
      subrefl_id,msg = sql.insert_record(database,'subreflector',subrefl_data)
      if subrefl_id < 0:
        if subrefl_id == -2:
          print "Duplicate entry:",subrefl_data
        else:
          print "Writing subreflector data failed;",
          print report_insert_error(subrefl_id,msg)
          return False
        
      # Now do the 'pointing' table and remember it's key
      pointing_cmd_data['UnixTime'] = UnixTime
      pointing_cmd_data['source_name_id'],src_id = \
        A.mysql.get_source_IDs(database,linedict['SOURCE'])
      # One never does pointing with a manual offset (I hope)
      pointing_cmd_data['beamaoff'] = 0.0
      pointing_cmd_data['beamxoff'] = 0.0
      pointing_cmd_data['beameoff'] = 0.0
      pointing_cmd_data['beamhoff'] = 0.0
      pointing_cmd_data['beamcoff'] = 0.0
      pointing_cmd_data['beamdoff'] = 0.0
      pointing_cmd_id,msg = sql.insert_record(database,'pointing_cmd',
                                          pointing_cmd_data)
      if pointing_cmd_id < 0:
        if pointing_cmd_id != -2:
          print "Writing pointing command data failed;",
          print sql.report_insert_error(pointing_cmd_id,msg)
          return False
        
      pointing_data['UnixTime'] = UnixTime
      pointing_data['azimuth'] = float(linedict['AZ'])
      pointing_data['elevation'] = float(linedict['EL'])
      pointing_data['pointing_cmd_id'] = pointing_cmd_id
      pointing_id,msg = sql.insert_record(database,'pointing',pointing_data)
      if pointing_id < 0:
        if pointing_id != -2:
          print "Writing pointing data failed;",
          print sql.report_insert_error(pointing_id,msg)
          return False
        else:
          UnixTime_str = "%12.1f" % UnixTime
          if diag:
            print "Duplicate entry for",UnixTime_str, \
              "while writing pointing command data"
          query = "SELECT id FROM pointing WHERE UnixTime = "+UnixTime_str+";"
          result = sql.ask_db(database,query)
          if diag:
            print query
            print "Obtained pointing ID",result
          if result == ():
            print "Skipped data for",UnixTime_str,T.ctime(UnixTime)
      else:
        # Now we can do the 'five_pt' table
        pointing_id = result[0][0]
        five_pt_data['pointing_id'] = pointing_id
        five_pt_data['freq'] = linedict['FREQ']
        five_pt_data['axhxpos'] = linedict['AXHXPOS']
        five_pt_data['eldepos'] = linedict['ELDEPOS']
        five_pt_data['deltaXp'] = linedict['DELTAXP']
        five_pt_data['deltaX0'] = linedict['DELTAX0']
        five_pt_data['deltaXm'] = linedict['DELTAXM']
        five_pt_data['baseXp'] = linedict['BASEXP']
        five_pt_data['baseXm'] = linedict['BASEXM']
        five_pt_data['X_bw'] = linedict['X_BW']
        five_pt_data['X_err'] = linedict['X_ERR']
        five_pt_data['Xpos'] = linedict['XPOS']
        five_pt_data['deltaYp'] = linedict['DELTAYP']
        five_pt_data['deltaY0'] = linedict['DELTAY0']
        five_pt_data['deltaYm'] = linedict['DELTAYM']
        five_pt_data['baseYp'] = linedict['BASEYP']
        five_pt_data['baseYm'] = linedict['BASEYM']
        five_pt_data['Y_bw'] = linedict['Y_BW']
        five_pt_data['Y_err'] = linedict['Y_ERR']
        five_pt_data['Ypos'] = linedict['YPOS']
        five_pt_data['subreflector_id'] = subrefl_id
        five_pt_id,msg = sql.insert_record(database,'five_pt',five_pt_data)
        if five_pt_id < 0:
          if five_pt_id == -2:
            print "Duplicate entry in table 'five_pt':",five_pt_data
          else:
            print "Writing five_pt data failed;",
            print report_insert_error(five_pt_id,msg)
            return False
        else:
          if diag:
            print "Wrote five_pt record",five_pt_id

        # Now we can do the 'boresight' table
        boredata['pointing_id'] = pointing_id
        boredata['xeloff'] = linedict['XELOFF']
        boredata['eloff'] =linedict['ELOFF']
        boredata['xdecoff'] = linedict['XDECOFF']
        boredata['decoff'] = linedict['DECOFF']
        boredata['five_pt_id'] = five_pt_id
        boresight_id,msg = sql.insert_record(db,'boresight',boredata)
        if boresight_id < 0:
          if boresight_id == -2:
            print "Duplicate entry:",boredata
          else:
            print "Writing boresight data failed;",
            print report_insert_error(boresight_id,msg)
            return False
  return True

def process_minical_file(db,fname,dss):
  """Process EAC minical file"""
  keys = []
  fd = open(fname,'r')
  lines = fd.readlines()
  fd.close()
  # First process the header
  newpars = False # parameters not finished until the column headers found
  for l in lines:
    line = l.strip()
    if len(line) == 0:
      if diag:
        print "Blank line skipped"
    else:
      parts = line.split()
      if parts[0] == 'Year':
        year = int(parts[1])
      elif parts[0] == 'Dss':
        dss = int(parts[1])
      elif parts[0] == 'Frequency':
        freq = float(parts[1])
      elif parts[0] == 'Bandwidth':
        bw = float(parts[1])
      elif parts[0] == 'Sample_time(sec)':
        samp_time = float(parts[1])
      elif parts[0] == 'Samples/reading':
        num_samps = int(parts[1])
      elif parts[0] == 'Reads/point':
        num_reads = int(parts[1])
      elif parts[0] == 'DOY':
        newpars = True # triggers a new record for the pars file
        for key in parts:
          keys.append(key)
      else:
        # must be a data line; parse the data
        linedict = {}
        for index in range(len(keys)):
          linedict[keys[index]] = parts[index]
        start_sec = float(linedict['UT'])*3600
        UnixTime = DT.VSR_tuple_to_timestamp(year,
                                             int(linedict['DOY']),
                                             start_sec)
        if newpars == True:
          # We can write a record to 'minical_pars'
          fields = {}
          fields['UnixTime'] = UnixTime
          fields['Freq'] = freq
          fields['Bandwidth'] = bw
          fields['sample_time'] = samp_time
          fields['samples_per_reading'] = num_samps
          fields['readings_per_point'] = num_reads
          minical_pars_id,msg = sql.insert_record(db,'minical_pars',fields)
          if minical_pars_id < 0:
            if minical_pars_id == -2:
              print "Duplicate entry in minical_pars:",fields
            else:
              print "Writing minical pars data failed;",
              print report_insert_error(minical_pars_id,msg)
              return False
          # We don't need to process a header until we find another one
          newpars = False
        # Now we can process the data file's lines
        #   First save the weather data
        fields = {}
        fields['UnixTime'] = UnixTime
        fields['Press'] = float(linedict['Press'])
        fields['WxTemp'] = float(linedict['WxTemp'])
        fields['Hum'] = float(linedict['Hum'])
        fields['Wind_vel'] = float(linedict['Wind_vel'])
        fields['Wind_dir'] = float(linedict['Wind_dir'])
        weather_id,msg = sql.insert_record(db,'venus.weather',fields)
        if weather_id < 0:
          if minical_pars_id == -2:
            print "Duplicate entry in weather:",fields
          else:
            print "Writing minical pars data failed;",
            print report_insert_error(weather_id,msg)
            return False
        #   Save the minical data
        fields = {}
        fields['UnixTime'] = UnixTime
        fields['minical_pars_id'] = minical_pars_id
        fields['Pol'] = linedict['Pol']
        fields['Gain'] = float(linedict['Gain'])
        fields['Linearity'] = float(linedict['Lin'])
        fields['El'] = float(linedict['El'])
        fields['T_lna'] = float(linedict['T(lna)'])
        fields['T_follow'] = float(linedict['T(follow)'])
        fields['R1'] = float(linedict['R1'])
        fields['R2'] = float(linedict['R2'])
        fields['R3'] = float(linedict['R3'])
        fields['R4'] = float(linedict['R4'])
        fields['R5'] = float(linedict['R5'])
        fields['T_load'] = float(linedict['T(load)'])
        fields['Tsys'] = float(linedict['Tsys'])
        fields['Tdiode'] = float(linedict['Tdiode'])
        fields['weather_id'] = weather_id
        minical_id,msg = sql.insert_record(db,'minical',fields)
        if minical_id < 0:
          if minical_id == -2:
            print "Duplicate entry:",fields
          else:
            print "Writing minical data failed:",
            print report_insert_error(minical_id,msg)
            return False
  return True

def process_minical_file(db,fname,dss):
  """Process EAC minical file"""
  keys = []
  fd = open(fname,'r')
  lines = fd.readlines()
  fd.close()
  # First process the header
  newpars = False # parameters not finished until the column headers found
  for l in lines:
    line = l.strip()
    if len(line) == 0:
      if diag:
        print "Blank line skipped"
    else:
      parts = line.split()
      if parts[0] == 'Year':
        year = int(parts[1])
      elif parts[0] == 'Dss':
        dss = int(parts[1])
      elif parts[0] == 'Frequency':
        freq = float(parts[1])
      elif parts[0] == 'Bandwidth':
        bw = float(parts[1])
      elif parts[0] == 'Sample_time(sec)':
        samp_time = float(parts[1])
      elif parts[0] == 'Samples/reading':
        num_samps = int(parts[1])
      elif parts[0] == 'Reads/point':
        num_reads = int(parts[1])
      elif parts[0] == 'DOY':
        newpars = True # triggers a new record for the pars file
        for key in parts:
          keys.append(key)
      else:
        # must be a data line; parse the data
        linedict = {}
        for index in range(len(keys)):
          linedict[keys[index]] = parts[index]
        start_sec = float(linedict['UT'])*3600
        UnixTime = DT.VSR_tuple_to_timestamp(year,
                                             int(linedict['DOY']),
                                             start_sec)
        if newpars == True:
          # We can write a record to 'minical_pars'
          fields = {}
          fields['UnixTime'] = UnixTime
          fields['Freq'] = freq
          fields['Bandwidth'] = bw
          fields['sample_time'] = samp_time
          fields['samples_per_reading'] = num_samps
          fields['readings_per_point'] = num_reads
          minical_pars_id,msg = sql.insert_record(db,'minical_pars',fields)
          if minical_pars_id < 0:
            if minical_pars_id == -2:
              print "Duplicate entry in minical_pars:",fields
            else:
              print "Writing minical pars data failed;",
              print report_insert_error(minical_pars_id,msg)
              return False
          # We don't need to process a header until we find another one
          newpars = False
        # Now we can process the data file's lines
        #   First save the weather data
        fields = {}
        fields['UnixTime'] = UnixTime
        fields['Press'] = float(linedict['Press'])
        fields['WxTemp'] = float(linedict['WxTemp'])
        fields['Hum'] = float(linedict['Hum'])
        fields['Wind_vel'] = float(linedict['Wind_vel'])
        fields['Wind_dir'] = float(linedict['Wind_dir'])
        weather_id,msg = sql.insert_record(db,'venus.weather',fields)
        if weather_id < 0:
          if minical_pars_id == -2:
            print "Duplicate entry in weather:",fields
          else:
            print "Writing minical pars data failed;",
            print report_insert_error(weather_id,msg)
            return False
        #   Save the minical data
        fields = {}
        fields['UnixTime'] = UnixTime
        fields['minical_pars_id'] = minical_pars_id
        fields['Pol'] = linedict['Pol']
        fields['Gain'] = float(linedict['Gain'])
        fields['Linearity'] = float(linedict['Lin'])
        fields['El'] = float(linedict['El'])
        fields['T_lna'] = float(linedict['T(lna)'])
        fields['T_follow'] = float(linedict['T(follow)'])
        fields['R1'] = float(linedict['R1'])
        fields['R2'] = float(linedict['R2'])
        fields['R3'] = float(linedict['R3'])
        fields['R4'] = float(linedict['R4'])
        fields['R5'] = float(linedict['R5'])
        fields['T_load'] = float(linedict['T(load)'])
        fields['Tsys'] = float(linedict['Tsys'])
        fields['Tdiode'] = float(linedict['Tdiode'])
        fields['weather_id'] = weather_id
        minical_id,msg = sql.insert_record(db,'minical',fields)
        if minical_id < 0:
          if minical_id == -2:
            print "Duplicate entry:",fields
          else:
            print "Writing minical data failed:",
            print report_insert_error(minical_id,msg)
            return False
  return True

def get_Tsys_data(db,fname,dss):
  """This parses a system temperature file.  It uses takes the az/el
  data and adds it to the pointing table."""
  fd = open(fname,'r')
  lines = fd.readlines()
  if diag:
    print len(lines),"read"
  fd.close()
  # get the year from the file name, which has the form tYYYYDDD.X
  # where X is A or B.
  year = int(basename(fname)[1:5])
  # db = sql.open_db('dss'+str(dss),host,user,passwd)
  num_skipped = 0
  for l in lines:
    line = l.strip().split()
    if line[0] != 'DOY':
      fields = {}
      DOY = int(line[0])
      timestr = line[1]
      UnixTime = DT.VSR_script_time_to_timestamp(year,line[0]+"/"+line[1])
      if diag:
        print "get_Tsys_data: UnixTime =", UnixTime
      # Get a pointer to the prior command record in the pointing_cmd table
      prior_record_id = find_prior_pointing_command(db,dss,UnixTime)
      if diag:
        print "get_Tsys_data: prior pointing command ID =", prior_record_id
      # Now modify the record
      fields['UnixTime'] = UnixTime
      fields['azimuth'] = float(line[7])
      fields['elevation'] = float(line[8])
      if prior_record_id > 0:
        fields['pointing_cmd_id'] = prior_record_id
      pointing_id,msg = sql.insert_record(db,'dss'+str(dss)+'.pointing',fields)
      if pointing_id == -2:
        # Duplicate entry
        if diag:
          print "get_Tsys_data: Duplicate pointing record for", \
            UnixTime, T.ctime(UnixTime)
        # Get the ID of the already stored record
        response = sql.ask_db(db,
          "select id from dss"+str(dss)+".pointing where UnixTime = " +
          str(UnixTime) + ";")
        if response == ():
          print "get_Tsys_data: Empty response from 'select id from pointing'"
          print "Record for",str(UnixTime),T.ctime(UnixTime),"skipped"
          num_skipped += 1
        else:
          pointing_id = response[0][0]
          if diag:
            print "exists as record",pointing_id
      elif pointing_id < 0:
        print "Problem storing pointing record; error", pointing_id
        print msg
      else:
        if diag:
          print "Pointing data for",T.ctime(UnixTime),"stored in record", \
            pointing_id
      # Since pointing_id may have changed (if originally -2), test again
      if pointing_id > 0:
        # Now the Tsys table
        fields = {}
        fields['UnixTime'] = UnixTime
        fields['freq'] = float(line[4])
        fields['pol'] = "'"+line[5]+"'"
        fields['tsys'] = line[2]
        fields['pointing_id'] = pointing_id
        tsys_id,msg =  sql.insert_record(db,'dss'+str(dss)+'.tsys',fields)
        if tsys_id == -2:
          # duplicate OK
          if diag:
            print "Duplicate Tsys record for",T.ctime(UnixTime)
            msg_parts = msg.split("'")
            print "conflicts with record for",T.ctime(float(msg_parts[1]))
        elif tsys_id < 0:
          print "Problem inserting into Tsys table: error", tsys_id
          print msg
        elif diag:
          print "wrote tsys record",tsys_id
  if diag:
    print num_skipped,"records skipped"
  return True
