# -*- coding: utf-8 -*-
"""
Global Python convenience functions for MySQL

The first section has functions that are handy on the command line
for getting a peek at the structure of databases. (It's hard to imagine
programming uses for 'get_databases' and 'get_public_tables')::
  get_databases(host,user,pw)
  show_databases(host,user,passwd)
  check_database(host,user,passwd,database)

In the second group a connection is created and used until closed.
'open_db' and 'ask_db' handle errors. These are useful for trying queries by
hand before serious coding. The functions defined in this section are::
  open_db(db,host,user,passwd)
  ask_db(connection,query_string)
  check_table(db_conn,table)
"""
import logging
from Data_Reduction.DSN.GAVRT.Mysql import open_db

loglevel = logging.WARNING

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
  print("Databases on %s" % host)
  for db in dbs:
    print "  ",db[0]
  for line in dbs:
    if line[0] != 'information_schema' and line[0] != 'mysql' and \
      not re.search('wiki',line[0]):
      print "Tables in",line[0]
      db = open_db(line[0],host,user,passwd)
      tbs = get_public_tables(db)
      logging.debug(str(tbs))
      for tb in tbs:
        try:
          reserved_index = reserved.index(tb[0].upper())
          query = "SELECT COUNT(*) FROM `"+tb[0]+"`;"
        except ValueError,details:
          query = "SELECT COUNT(*) FROM "+tb[0]+";"
        result = ask_db(db,query)
        print "  ",tb[0],"has",result[0][0],"rows"
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

def check_database(host, user, passwd, database):
  """
  Does the database exist on the host?

  @param host : name or IP address of database host
  @type  host : str

  @param user : an authorized database user
  @type  user : str

  @param passwd : user's password
  @type  passwd : str

  @param database : database to be opened
  @type  database : str

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
#######################################################################

def open_db(db, host, user, passwd):
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

def check_table(db_conn, table):
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
    print "check_table: did not understand response"
    print q
    return False
