"""
mysql_tables.py - tools for making MySQL tables

In general, it should not be necessary to create tables with Python.  'dia'
should have been used to create an SQL file which creates all the tables as
part of an overall scheme.

Contains::
  create_table(database,name,keys)
  make_table(database,table_data)
  get_table_columns(db_conn,table)
"""
diag = False

import MySQLdb
import Data_Reduction.DSN.GAVRT.Mysql as sql
import sys
  
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
  for key in key_types.keys():
    try:
      colname = sql.reserved(key)
    except:
      colname = '`'+key+'`'
    create_string = create_string + ' '+colname+' '+key_types[key]+','
  create_string = create_string.rstrip(',') + ');'
  if diag:
    print "create_table:",create_string
  try:
    cursor.execute(create_string)
    return tbname
  except Exception, detail:
    print "create_table: execution failed"
    print detail
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
      print "make_table:",tbname,"exists"
    # get a cursor
    cursor = database.cursor()
  # Insert the data
  for datum in data:
    try:
      sql.insert_record(database,tbname,datum)
    except Exception, detail:
      # row exists
      print "make_table: Error for",datum
      print detail
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
