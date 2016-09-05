# -*- coding: utf-8 -*-
"""
Mysql - stuff to make using MySQLdb a little easier to use

This module has two sections.

The module's functions in the first section, mostly written by Glenn
Jones with some mods by Tom Kuiper, can be used to access a database
without OO programming, which some people may like.  You'll need 'open_db'
to get a connection object. The functions defined here are::
  insert_record(db,table,fields):
  update_record(db,table,fields,condition):
  get_as_dict(db,*args,**kwargs):
  get_last_id(db,table):
  get_last_record(db,table):
  get_record_by_id(db,table,id)
  rt_error(error,msg)

The second section defines the class BaseDB.  It can be used as is if
you give it all the connection parameters.  It may also serve as a super
class for subclasses which are specific to certain servers, databases
and users.
"""
import MySQLdb  #  This causes a warning; see above
import numpy as np
import re
import logging

reserved = ['ADD', 'ALL', 'ALTER', 'ANALYZE', 'AND', 'AS', 'ASC',
            'ASENSITIVE', 'BEFORE', 'BETWEEN', 'BIGINT', 'BINARY',
            'BLOB', 'BOTH', 'BY', 'CALL', 'CASCADE', 'CASE', 'CHANGE',
            'CHAR', 'CHARACTER', 'CHECK', 'COLLATE', 'COLUMN',
            'CONDITION', 'CONSTRAINT', 'CONTINUE', 'CONVERT', 'CREATE',
            'CROSS', 'CURRENT_DATE', 'CURRENT_TIME', 'CURRENT_TIMESTAMP',
            'CURRENT_USER', 'CURSOR', 'DATABASE', 'DATABASES', 'DAY_HOUR',
            'DAY_MICROSECOND', 'DAY_MINUTE', 'DAY_SECOND', 'DEC',
            'DECIMAL', 'DECLARE', 'DEFAULT', 'DELAYED', 'DELETE', 'DESC',
            'DESCRIBE', 'DETERMINISTIC', 'DISTINCT','DISTINCTROW', 'DIV',
            'DOUBLE', 'DROP', 'DUAL', 'EACH', 'ELSE', 'ELSEIF', 'ENCLOSED',
            'ESCAPED', 'EXISTS', 'EXIT', 'EXPLAIN', 'FALSE', 'FETCH',
            'FLOAT', 'FLOAT4', 'FLOAT8', 'FOR', 'FORCE', 'FOREIGN',
            'FROM', 'FULLTEXT', 'GRANT', 'GROUP', 'HAVING', 'HIGH_PRIORITY',
            'HOUR_MICROSECOND', 'HOUR_MINUTE', 'HOUR_SECOND', 'IF',
            'IGNORE', 'IN', 'INDEX', 'INFILE', 'INNER', 'INOUT',
            'INSENSITIVE', 'INSERT', 'INT', 'INT1', 'INT2', 'INT3', 'INT4',
            'INT8', 'INTEGER', 'INTERVAL', 'INTO', 'IS', 'ITERATE', 'JOIN',
            'KEY', 'KEYS', 'KILL', 'LEADING', 'LEAVE', 'LEFT', 'LIKE',
            'LIMIT', 'LINES', 'LOAD', 'LOCALTIME', 'LOCALTIMESTAMP',
            'LOCK', 'LONG', 'LONGBLOB', 'LONGTEXT', 'LOOP', 'LOW_PRIORITY',
            'MATCH', 'MEDIUMBLOB', 'MEDIUMINT', 'MEDIUMTEXT', 'MIDDLEINT',
            'MINUTE_MICROSECOND', 'MINUTE_SECOND', 'MOD', 'MODIFIES',
            'NATURAL', 'NOT', 'NO_WRITE_TO_BINLOG', 'NULL', 'NUMERIC', 'ON',
            'OPTIMIZE', 'OPTION', 'OPTIONALLY', 'OR', 'ORDER', 'OUT',
            'OUTER', 'OUTFILE', 'PRECISION', 'PRIMARY', 'PROCEDURE',
            'PURGE', 'READ', 'READS', 'REAL', 'REFERENCES', 'REGEXP',
            'RELEASE', 'RENAME', 'REPEAT', 'REPLACE', 'REQUIRE', 'RESTRICT',
            'RETURN', 'REVOKE', 'RIGHT', 'RLIKE', 'SCHEMA', 'SCHEMAS',
            'SECOND_MICROSECOND', 'SELECT', 'SENSITIVE', 'SEPARATOR',
            'SET', 'SHOW', 'SMALLINT', 'SONAME', 'SPATIAL', 'SPECIFIC',
            'SQL', 'SQLEXCEPTION', 'SQLSTATE', 'SQLWARNING',
            'SQL_BIG_RESULT', 'SQL_CALC_FOUND_ROWS', 'SQL_SMALL_RESULT',
            'SSL', 'STARTING', 'STRAIGHT_JOIN', 'TABLE', 'TERMINATED',
            'THEN', 'TINYBLOB', 'TINYINT', 'TINYTEXT', 'TO', 'TRAILING',
            'TRIGGER', 'TRUE', 'UNDO', 'UNION', 'UNIQUE', 'UNLOCK',
            'UNSIGNED', 'UPDATE', 'USAGE', 'USE', 'USING', 'UTC_DATE',
            'UTC_TIME', 'UTC_TIMESTAMP', 'VALUES', 'VARBINARY', 'VARCHAR',
            'VARCHARACTER', 'VARYING', 'WHEN', 'WHERE', 'WHILE', 'WITH',
            'WRITE', 'XOR', 'YEAR_MONTH', 'ZEROFILL']

class MysqlException(Exception):
  """
  Handles exceptions in local module Mysql
  """
  def __init__(self,message,*args):
    self.message = message
    self.args = args
  def __str__(self):
    return (self.message % self.args)

############################ Global Record Functions ##########################

def insert_record(db,table,fields):
  """
  Insert a record into data base
  
  @param db : return from a MySQLdb connect()
  @type  db : database connection object
  
  @param table : table into which the record will be inserted
  @type  table : str

  @param fields : dictionary with column names and the values for the record
  @type  fields : dict

  @return: unique ID of the record, message which is blank if successful.
  """
  keystr = ''
  valstr = ''
  values = []
  try:
    del fields['ID']
  except:
    pass
  for k,v in fields.items():
    keystr = keystr + k + ', '
    valstr = valstr + '%s ,'
    values.append(v)
  keystr = keystr[:-2]    #strip final ,
  valstr = valstr[:-2]
  c = db.cursor()
  query = "INSERT INTO "+ table + " (" + keystr + ") VALUES (" + valstr + ");"
  try:
    c.executemany(query,[tuple(values)])
    c.close()
    new_ids = int(db.insert_id())
    db.commit()
    return (new_ids,"")
  except Exception, details:
    c.close()
    if details[0] == 1062:
      # Duplicate entry
      print 'insert_record: Duplicate entry for'
      print query % tuple(values)
      return (-2,details[1])
    elif details[0] == 1025:
      # Lock timeout exceeded; retry?
      print 'insert_record: Timeout exceeded'
      return (-3,details[1])
    elif details[0] == 1146:
      print "\ninsert_record: table does not exist"
      return (-4,details[1])
    else:
      print "\ninsert_record: could not execute:"
      print query % tuple(values)
      print "Returned error",details[0]
      return (-1,details[1])

def update_record(db,table,fields,condition):
  """
  Updates a record in 'table' of data base 'db'
  
  @param db : database connection object returned from a MySQLdb connect()
  
  @param table : name of the table in which the record will be updated
  
  @param fields : dictionary with column names and the values to be updated
  
  @param condition :  for the condition which selects the record(s)
  @type  condition : (column name, value tuple)
  
  @return: unique ID of the record, message which is blank if successful.
  """
  keystr = ''
  valstr = ''
  values = []
  try:
    del fields['ID']
  except:
    pass
  query = "UPDATE "+ table + " SET "
  for k,v in fields.items():
    query += k + " = " + "%s" +", "
    values.append(v)
  query = query[:-2]    #strip final ,
  query += " WHERE " + condition[0] + " = " + str(condition[1]) + ";"
  try:
    c = db.cursor()
  except Exception,details:
    logging.error("Could not get cursor")
    return (-1,details)
  try:
    c.executemany(query,[tuple(values)])
    # This should return ()
    result = c.fetchall()
  except Exception, details:
    logging.error("insert_record: could not execute: "+query,tuple(values))
    result = (-2,details)
  logging.debug("update_record: "+result)
  c.close()
  db.commit()
  return result

def get_as_dict(db, *args, **kwargs):
  """
  Executes a query of the database

  At present, only keyword {'asfloat': True} is recognized and is the default
  if not given.  It will convert to float any values for which it is possible.

  If the query returns multiple rows, each value associated with a keyword will
  be a list. If nothing was found, an empty dictionary is returned.
  
  @param db : database connection object
  
  @param args: query to be executed
  
  @param kwargs : a dictionary with keyword arguments.
  
  @eturn: the record as a dictionary.
  """
  try:
    asfloat = kwargs['asfloat']
  except:
    asfloat = True
  c = db.cursor()
  c.execute(*args)
  res = np.array(c.fetchall())
  logging.debug("get_as_dict result shape: "+str(res.shape))
  logging.debug("get_as_dict result: "+str(res))
  # an empty result will have a shape of (1,)
  if len(res.shape) > 1:
    descr = [x[0] for x in c.description]
    if asfloat:
      rd = {}
      for x in range(res.shape[1]):
        try:
          r = res[:,x].astype('float')
        except:
          r = res[:,x]
        rd[descr[x]] = r
    else:
      rd = dict([(descr[x],res[:,x]) for x in range(res.shape[1])])
  else:
    rd = {}
  return rd

def get_last_id(db,table):
  """
  Returns the ID (integer) of the last record
  
  @param db : database connection object
  
  @param table : name of the table
  @type  table : str

  @return: int
  """
  c = db.cursor()
  c.execute("SELECT ID FROM " + table + " ORDER BY ID DESC LIMIT 1;")
  return int(c.fetchone()[0])

def get_last_record(db,table):
  """
  Returns the last record as a dictionary
  
  @param db : database connection object
  
  @param table : name of the table
  @type  table : str

  @return: record (dict)
  """
  c = db.cursor()
  # This returns the last row
  c.execute("SELECT * FROM " + table + " ORDER BY ID DESC LIMIT 1;")
  res = c.fetchone()
  # This returns the column names
  descr = [x[0] for x in c.description]
  # This returns the row as a dictionary
  return dict(zip(descr,res))

def get_record_by_id(db,table,id):
  """
  Returns the record with the given ID
  
  @param db : database connection object
  
  @param table : table name
  @type  table : str

  @param id : row ID
  @type  id : int

  @return: dict
  """
  c = db.cursor()
  c.execute("SELECT * FROM " + table + " WHERE ID = %s;",(id,))
  res = c.fetchone()
  descr = [x[0] for x in c.description]
  return dict(zip(descr,res))

def report_insert_error(error,msg):
  """
  Returns text string with cause of the insert error
  """
  if error == -3:
    return "Lock timeout exceeded"
  elif error == -2:
    return "Duplicate entry:"+msg
  elif error == -1:
    return "method 'report_insert' error:"+msg
  else:
    return "Unknown error:"+error+"; "+msg

############################## Class BaseDB ############################

class BaseDB():
  """
  This is a database superclass.
  
  Public attributes::
    db -   the database connection object
    c  -   the database cursor
    host - database host
    port - port used to connect to the database
    pw -   user's password
    user - authorized db user
  Methods::
    connect - returns a connection to a database.
    check_db - reconnects to the database if a connection has been lost
    cursor   - """
  def __init__(self,host,user,pw,db,port=3306):
    """
    Initializes a BaseDB instance by connecting to the database
    
    @param host : the IP address as a string, or fully qualified host name
    @type  host : str
    
    @param user : a valid user on the mysql server at host
    
    @param pw : user's password or "" if not required
    
    @param db : database name (string)
    
    Generates a cursor object BaseDB.c.
    """
    self.host = host
    self.port = port
    self.user = user
    self.pw = pw
    self.db = db
    self.connect()

  def connect(self):
    """
    Make a connection to the database. Creates a cursor object.
    
    Automatically invoked when an instance is created; can be called
    again if the connection is closed but the database object persists.
    """
    self.db = MySQLdb.connect(host = self.host,
                              port=self.port,
                              user=self.user,
                              passwd=self.pw,
                              db=self.db,
                              compress=True)
    self.c = self.db.cursor()

  def close(self):
    """
    Close a connection
    """
    self.c.close()
    
  def checkDB(self):
    """
    Reconnects to the database if the connection has been lost.
    """
    try:
      self.db.commit()
    except:
      self.connect()
    self.db.commit()

  def cursor(self):
    """
    Creates a database cursor object; same as BaseDB.c but this
    is better because it handles disconnected a database
    """
    self.checkDB()
    return self.db.cursor()

  def commit(self):
    """
    Commits the most recent database transaction
    """
    return self.db.commit()
        
  def insertRecord(self,table,rec):
    """
    Inserts a record into the database; handles a disconnected database

    @param table : table name
    @type  table : str

    @param rec : a dictionary with column names as keys.

    @return: record ID (int)
    """
    self.checkDB()
    return insert_record(self.db, table, rec)
  
  def getLastId(self,table):
    """
    ID of the last record
    
    @param table : the name of the table
    @type  table : str

    @return: ID (int)
    """
    self.checkDB()
    return get_last_id(self.db,table)
  
  def getLastRecord(self,table):
    """
    Returns the last record as a dictionary
    
    @param table : name of the table (string)

    @return: dict
    """
    self.checkDB()
    return get_last_record(self.db,table)
  
  def getRecordById(self,table,id):
    """
    Get the record with the given ID
    
    @param table : table name
    @type  table : str
    
    @param id : row ID
    @type  id : int

    @return: dict
    """
    self.checkDB()
    return get_record_by_id(self.db,table,id)
    
  def get(self,*args):
    """
    Executes a query of the database
    
    @param args : query to be executed
    
    @return: record (dict)
    """
    self.checkDB()
    return get_as_dict(self.db,*args,**dict(asfloat=True))
        
  def updateValues(self,vald,table):
    """
    Add row with updated values
    
    Add a new row to table with same values as previous row,
    except for keys in vald, which are updated with provided
    values.  This is useful for keeping logs.

    @param vald : updated values
    @type  vald : dict

    @param table : table name
    @type  table : str
    """
    lastrec = self.getLastRecord(table)
    lastrec.update(vald)
    self.insertRecord(table, lastrec)

  def get_public_tables(self):
    """
    List the table names.in the database.

    @return: tuple of tuples of str
    """
    try:
      self.c.execute("""SHOW TABLES;""")
      result = self.c.fetchall()
    except MySQLdb.Error, e:
      print "MySQLdb error: Cannot connect to server"
      print "Error code:",e.args[0]
      print "Error message:",e.args[1]
      result = None
    return result

  def report_table(self,table,columns):
    """
    Reports on the columns in a table

    @param table : table name
    @type  table : str

    @param columns : list of column names
    @type  columns : list of str

    @return: result of query
    """
    print "\nShowing name and type of columns for "+table
    response = self.get("show columns from "+table+";")
    for key in columns:
      print response[key]
    return response

  def get_rows_by_date(self,table,columns,year,doy,utcs=None):
    """
    Gets data from xscan table

    @param table : table name
    @type  table : str

    @param columns : list of columns to be selected
    @type  columns : list of str

    @type year : int

    @param doy : day of year
    @type  doy : int

    @param utcs : not used

    @return: dict of numpy arrays keyed on column name
    """
    columnstr = str(columns).lstrip('[').rstrip(']').replace("'","")
    try:
      response = self.get("select " + columnstr
                        + " from "+table+" where year=%s and doy=%s",
                        (year,doy))
    except Mysql.MySQLdb.OperationalError, details:
      print "MySQLdb OperationalError:",details
    else:
      return response

  def get_rows_by_time(self,table,columns,year,doy,utcs):
    """
    Queries a table for quantities in columns at designated times

    It takes the first row it finds matching the date and time.  So it has an
    effective resolution of one second.

    @param table : table name
    @type  table : str

    @param columns : list of columns to be selected
    @type  columns : list of str

    @type year : int

    @param doy : day of year
    @type  doy : int

    @param utcs : times to be selected; first occurrence is used
    @type  utcs : list of unixtimes (seconds since the epoch)

    @return: dict of numpy arrays keyed on column name
    """
    data = {}
    for col in columns:
      data[col] = []
    columnstr = str(columns).lstrip('[').rstrip(']').replace("'","")
    for utc in utcs:
      fmt = "select "+columnstr+" from "+table+" where year=%s and doy=%s and utc=%s;"
      result = self.get(fmt,(year,doy,utc))
      for col in columns:
        data[col].append(result[col][0])
    for col in columns:
      data[col] = np.array(data[col])
    return data
