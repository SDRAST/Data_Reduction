# -*- coding: utf-8 -*-
import pickle
import os

import Data_Reduction.DSN.GAVRT.Mysql as Mysql

host,user,pw = pickle.load(open(os.environ['HOME']+"/.GAVRTlogin.p", "rb" ))

def connect_gavrt():
  db = Mysql.BaseDB(host, user, pw, 'dss28_eac')
  print ("Connected to",db.host,"as",db.user,"with 'db'")
  return db
