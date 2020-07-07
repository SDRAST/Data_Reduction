# -*- coding: utf-8 -*-
"""
Classes for GAVRT data using MySQL tables

This module defines the classes BaseDB and MysqlException.  It also provides
subclasses for reducing data stored in the DSS-28 database.
"""
import logging
import MySQLdb
import numpy as np
import os
import pickle

logger = logging.getLogger(__name__)


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

############################## Classes #############################

class MysqlException(Exception):
  """
  Handles exceptions in local module Mysql
  """
  def __init__(self,message,*args):
    self.message = message
    self.args = args
  def __str__(self):
    return (self.message % self.args)

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
    cursor   - 
  """
  def __init__(self, host=None, user=None , pw=None, name=None, port=3306):
    """
    Initializes a BaseDB instance by connecting to the database
    
    This defaults to the GAVRT DSS-28 EAC database on 'kracken'.  However, use
    of sub-class DSS28db is recommended to avoid writing acciden
    
    @param host : the IP address as a string, or fully qualified host name
    @type  host : str
    
    @param user : a valid user on the mysql server at host
    
    @param pw : user's password or "" if not required
    
    @param db : database name (string)
    
    Generates a cursor object BaseDB.c.
    """
    self.logger = logging.getLogger(logger.name+".BaseDB")
    if name and host and user and pw:
      self.name= name
      self.host = host
      self.port = port
      self.user = user
      self.pw = pw
      self.connect()
    else:
      self.logger.error("__init__: host, user, pw and DB name required")
      raise MysqlException("Missing arguments host=%s, user=%s, name=%s or missing pw?",
                           host, user, name)

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
                              db=self.name,
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

  def send_query(self,query):
    """
    Send a MySQL query

    @param crsr : cursor object for a database

    @param query : MySQL query
    @type  query : str

    @return: str with query response
    """
    query += ";"
    self.c.execute(query)
    response = self.c.fetchall()
    return response

  def commit(self):
    """
    Commits the most recent database transaction
    """
    return self.db.commit()
        
  def insertRecord(self, table, rec):
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
    self.c.execute("SELECT * FROM " + table + " ORDER BY ID DESC LIMIT 1;")
    res = self.c.fetchone()
    # This returns the column names
    descr = [x[0] for x in self.c.description]
    # This returns the row as a dictionary
    return dict(list(zip(descr,res)))
  
  def getRecordById(self,table,rec_id):
    """
    Get the record with the given ID
    
    @param table : table name
    @type  table : str
    
    @param id : row ID
    @type  id : int

    @return: dict
    """
    self.checkDB()
    self.c.execute("SELECT * FROM " + table + " WHERE ID = %s;",(rec_id,))
    res = self.c.fetchone()
    descr = [x[0] for x in self.c.description]
    return dict(list(zip(descr,res)))
    
  def get(self,*args):
    """
    Executes a query of the database
    
    @param args : query to be executed
    
    @return: record (dict)
    """
    self.checkDB()
    self.c.execute(*args)
    result = np.array(self.c.fetchall())
    return result
    
  def get_as_dict(self, *args, **kwargs):
    """
    Executes a query of the database and returns the result as a dict
  
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
    res = self.get(*args)
    # an empty result will have a shape of (1,)
    if len(res.shape) > 1:
      descr = [x[0] for x in self.c.description]
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
        
  def updateValues(self, vald, table):
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
    List the table names in the database.

    @return: tuple of tuples of str
    """
    try:
      self.c.execute("""SHOW TABLES;""")
      result = self.c.fetchall()
    except MySQLdb.Error as e:
      self.logger.error(
                  "get_public_tables: MySQLdb error: Cannot connect to server")
      self.logger.error("get_public_tables: error code:",e.args[0])
      self.logger.error("get_public_tables: error message:",e.args[1])
      result = None
    return result

  def get_columns(self, table):
    """
    Returns information about the columns of a table
    """
    self.c.execute("show columns from "+table+";")
    return self.c.fetchall()
  
  def get_data_index(self):
    """
    """
    tables = self.get_public_tables()
    tableindex = {}
    for table in tables:
      tb = table[0]
      tableindex[tb] = self.get_columns(tb)
    return tableindex
    
  def report_table(self,table,columns):
    """
    Reports on the columns in a table
    
    Response has keys: Extra, Default, Field, Key, Null, Type

    @param table : table name
    @type  table : str

    @param columns : list of column names
    @type  columns : list of str

    @return: result of query
    """
    self.logger.info("report_table: showing name and type of columns for %s",
                     table)
    response = self.get("show columns from "+table+";")
    self.logger.debug("report_table: response: %s", response)
    report = []
    for row in response:
      name = row[0]
      for column in columns:
        if name == column:
          report.append(row)
    return report

  def get_rows_by_date(self, table, columns, year, doy):
    """
    Gets data from a table

    @param table : table name
    @type  table : str

    @param columns : list of columns to be selected
    @type  columns : list of str

    @type year : int

    @param doy : day of year
    @type  doy : int

    @return: dict of numpy arrays keyed on column name
    """
    columnstr = str(columns).lstrip('[').rstrip(']').replace("'","")
    try:
      response = self.get_as_dict("select " + columnstr
                        + " from "+table+" where year=%s and doy=%s",
                        (year,doy))
    except Mysql.MySQLdb.OperationalError as details:
      print("MySQLdb OperationalError:",details)
    else:
      return response

  def get_rows_by_time(self,table,columns,year,doy,utcs):
    """
    Queries a table for quantities in columns at designated times

    This loops over a list of UTs.  For each, it takes the first row it finds 
    matching the date and time.  So it has an effective resolution of one second.

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
      result = self.get_as_dict(fmt,(year,doy,utc))
      for col in columns:
        data[col].append(result[col][0])
    for col in columns:
      data[col] = np.array(data[col])
    return data

############################ Global Functions ##########################

def get_databases(host, user, pw):
  """
  Command line function to recall the database names on a host

  @param host : host.domain or IP address
  @type  host : str

  @param user : user login name
  @type  user : str

  @param pw : password
  @type  pw : str

  @return: nested tuple of tuples of str
  """
  conn = open_db("",host, user, pw)
  response = ask_db(conn,"""SHOW DATABASES;""")
  conn.close()
  return response

def show_databases(host, user, passwd):
  """
  Prints report on all the databases

  @param host : host name
  @type  host : str

  @param user : user name
  @type  user : str

  @param passwd : user's password
  @type  passwd : str

  @return: printed report all the tables in each database.
  """
  dbs = get_databases(host, user, passwd)
  print(("Databases on %s" % host))
  for db in dbs:
    print("  ",db[0])
  for line in dbs:
    if line[0] != 'information_schema' and line[0] != 'mysql' and \
      not re.search('wiki',line[0]):
      print("Tables in",line[0])
      db = open_db(line[0],host,user,passwd)
      tbs = get_public_tables(db)
      logging.debug(str(tbs))
      for tb in tbs:
        try:
          reserved_index = reserved.index(tb[0].upper())
          query = "SELECT COUNT(*) FROM `"+tb[0]+"`;"
        except ValueError as details:
          query = "SELECT COUNT(*) FROM "+tb[0]+";"
        result = ask_db(db,query)
        print("  ",tb[0],"has",result[0][0],"rows")
      db.close()
  logging.info("Disconnected from server on %s",host)

def check_database(host, user, passwd, database):
  """
  Does the database exist on the host?

  @param db : database to be opened
  @type  db : str

  @param host : name or IP address of database host
  @type  host : str

  @param user : an authorized database user
  @type  user : str

  @param passwd : user's password
  @type  passwd : str

  @return: True or False
  """
  logging.info("check_database: getting databases")
  dbs = get_databases(host, user, passwd)
  exists = False
  for line in dbs:
    if line[0] == database:
      exists = True
      break
  logging.warning("check_database: Does %s ? %s",database,str(exists))
  return exists

def report_insert_error(error, msg):
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
  host,user,pw = _validate_login(host, user, pw)
  if host and user and pw:
    GAVRTdb_login = (host,user,pw)
    fd = open(os.environ['HOME']+"/.GAVRTlogin.p", "wb")
    pickle.dump(GAVRTdb_login, fd)
    os.fchmod(fd, int(stat.S_IREAD | stat.S_IWRITE))
    fd.close()
    return True
  else:
    logger("create_GAVRT_mysql_auth: failed")
    return False

def _validate_login(host=None, user=None, pw=None):
  """
  """
  if host == None:
    host = _host # "kraken.gavrt.org"
  if user == None:
    user = _user # "ops"
  if pw == None:
    pw = _pw     # 'v3H!cle'
  return host, user, pw

############### not yet converted to BaseDB methods

def check_table(db_conn, table): ## convert to BaseDB method
  """
  Checks whether a table  exists.

  @param db_conn : database
  @type  db_conn : connection object

  @param table : table name
  @type  table : str

  @return: True or False
  """
  q = ask_db(db_conn,"SHOW TABLES LIKE '"+table+"';")
  if len(q) == 1:
    return True
  elif len(q) == 0:
    return False
  else:
    print("check_table: did not understand response")
    print(q)
    return False
    
def update_record(db, table, fields, condition): # not allowed in DSS28db
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
  for k,v in list(fields.items()):
    query += k + " = " + "%s" +", "
    values.append(v)
  query = query[:-2]    #strip final ,
  query += " WHERE " + condition[0] + " = " + str(condition[1]) + ";"
  try:
    c = db.cursor()
  except Exception as details:
    logging.error("Could not get cursor")
    return (-1,details)
  try:
    c.executemany(query,[tuple(values)])
    # This should return ()
    result = c.fetchall()
  except Exception as details:
    logging.error("insert_record: could not execute: "+query,tuple(values))
    result = (-2,details)
  logging.debug("update_record: "+result)
  c.close()
  db.commit()
  return result
  
def create_table(database,name,keys):
  """
  Creates a table with the specified column names and types.
  
  Here's an example::
  
    name = "customers"
    keys = {'name':    'CHAR(20 NOT NULL'),
            'job':     'VARCHAR(20)',
            'sex':     "ENUM('M','F')",
            'hobbies': "SET('chess','sailing','reading','knitting')",
            'birth':   'DATE',
            'balance': 'FLOAT'}
    create_table(someDatabase,name,keys)
    
  When a column name conflicts with a MySQL reserved word, the name is
  put in backward quotes.  This keeps the column name possibilities from
  being too restricted, e.g., `dec`.

  This exists from a time prior to when I started using 'tedia2sql'
  to create tables and may still be useful in simple situations.
  
  @param database : the database for the new table
  @param name :     the table name
  @param keys :     a dictionary with the column names and formats.
  """
  tbname = name
  key_types = keys
  # get a cursor
  cursor = database.cursor()
  # assemble the CREATE string
  create_string = 'CREATE TABLE '+tbname+' ('
  for key in list(key_types.keys()):
    try:
      colname = sql.reserved(key)
    except:
      colname = '`'+key+'`'
    create_string = create_string + ' '+colname+' '+key_types[key]+','
  create_string = create_string.rstrip(',') + ');'
  if diag:
    print("create_table:",create_string)
  try:
    cursor.execute(create_string)
    return tbname
  except Exception as detail:
    print("create_table: execution failed")
    print(detail)
    return None
  
def make_table(database,table_data):
  """
  Adds or populates a table in the specified database from a dictionary.

  Makes a table from the dictionary 'table_data' with the keys:: 
  
    name - the name of the table, a simple string with no spaces
    keys - a dictionary whose keys are the (lower-case) keys for the table
           and whose values are the corresponding data types.
    data - a list of dictionaries.  Each dictionary corresponds to
           a row in the table. The keys correspond to the keys defined in
           'keys'.
           
  This is useful when the table is created once and for all, or
  completely replaced. Here's a small example::
  
    {name: "customers",
     keys: {'name':    'CHAR(20 NOT NULL'),
            'job':     'VARCHAR(20)',
            'sex':     "ENUM('M','F')",
            'hobbies': "SET('chess','sailing','reading','knitting')",
            'birth':   'DATE',
            'balance': 'FLOAT'},
    data: [{name: "Sam", job: "carpenter", birth: 1950-07-21, balance: 5.25},
           {job: "nurse", name: "Sally", balance: 8.50, birth: 1960-3-15}]
  """
  tbname = table_data['name']
  key_types = table_data['keys']
  data = table_data['data']
  if not sql.check_table(database,tbname):
    if create_table(database,tbname,keys) == None:
      return None
  else:
    if diag:
      print("make_table:",tbname,"exists")
    # get a cursor
    cursor = database.cursor()
  # Insert the data
  for datum in data:
    try:
      sql.insert_record(database,tbname,datum)
    except Exception as detail:
      # row exists
      print("make_table: Error for",datum)
      print(detail)
      errtype, value, traceback = sys.exc_info()
      sys.excepthook(errtype,value,traceback)
  return tbname

def get_table_columns(db_conn,table):
  """
  Returns information about the columns in a table.
  Inputs::
  
    db_conn (connection object)
    table   (string)
    
  Outputs::
  
    a sequence of sequences with column information.  Each column sequence
    consisting of:
      column name       (string)
      column type       (string) such as 'tinyint(3) unsigned' or 'float'
      contains nulls    (string) 'YES' or 'NO'
      key               (string) is the column indexed?
      default value              e.g. 'None'
      extra information (string)
  """
  q = sql.ask_db(db_conn,"SHOW COLUMNS FROM "+table+";")
  return q
