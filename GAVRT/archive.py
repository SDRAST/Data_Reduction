"""
Mysql - stuff to make using MySQLdb a little easier to use

The module's functions, mostly written by Glenn Jones with some mods by 
Tom Kuiper, can be used to access a database without OO programming, 
which some people may like.  You'll need 'open_db' to get a connection 
object. The functions defined here are::
  insert_record(db,table,fields):
  update_record(db,table,fields,condition):
  get_as_dict(db,*args,**kwargs):
  get_last_id(db,table):
  get_last_record(db,table):
  get_record_by_id(db,table,id)
  rt_error(error,msg)
"""
import logging
import MySQLdb
import numpy as np

from Data_Reduction.GAVRT.mysql import _validate_login

logger = logging.getLogger(__name__)

def open_db(database='dss28_eac'):
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

def open_db(db, host, user, passwd): ### merge with the above
  """
  Open a MySQL database.

  Useful for simple command line SQL commands

  @param db : database to be opened
  @type  db : str

  @param host : name or IP address of database host
  @type  host : str

  @param user : an authorized database user
  @type  user : str

  @param passwd : user's password
  @type  passwd : str

  @return: a connection object with creating a database class
  """
  try:
    conn = MySQLdb.connect(db=db, host=host, user=user, passwd=passwd)
    if db == "":
      database = "mysql server"
    else:
      database = db
    logging.info("open_db: Connected to %s on %s",database,host)
  except MySQLdb.Error, e:
    print "open_db: Cannot connect to server"
    print "Error code:",e.args[0]
    print "Error message:",e.args[1]
    conn = None
  if diag:
    print "open_db: connection:",conn
  return conn

def ask_db(connection, query_string):
  """
  Command line function for a simple SQL query.

  @param connection : see open_db()
  @type  connection : database connection object
  @param query_string : MySQL query

  @return: the query result.
  """
  try:
    cursor = connection.cursor()
  except MySQLdb.Error, e:
    print "ask_db: Cannot get a cursor for this connection"
    print "Error code:",e.args[0]
    print "Error message:",e.args[1]
    result = None
  try:
    cursor.execute(query_string)
    result = cursor.fetchall()
  except MySQLdb.Error, e:
    print "ask_db: query failed"
    print query_string
    print "Error code:",e.args[0]
    print "Error message:",e.args[1]
    result = None
  cursor.close()
  connection.commit()
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

