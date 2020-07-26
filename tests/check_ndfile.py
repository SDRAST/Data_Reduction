"""
check structure of a text file

This test uses a t-file but any columnar text data file will work.
This tests the base Observation class.
"""

from Data_Reduction import DataGetterMixin, Map, Observation
from support.mixin import mixIn

import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# here we initiate a base class after mixing in the data getter
mixIn(Observation, DataGetterMixin)
obs = Observation(dss=28, date="2012/127", project="SolarPatrol")
obs.open_datafile('t12127.10', 
                          delimiter=[17,16,3,11,7,9,8,2,6], skip_header=1,
                      names="UTC Epoch Chan Tsys Int Az El Diode Level".split())
logger.info("Observation 'obs' opened with column widths specified")

# the data getter is already mixed in to Observation
obs2 = Observation(dss=28, date="2012/127", project="SolarPatrol")
obs2.open_datafile('t12127.10', skip_header=1,
             names="Year DOY UTC Epoch Chan Tsys Int Az El Diode Level".split())
logger.info("Observation 'obs2' opened with assigned column names")

# the class Map inherits from DataGetterMixin, so no explicit mixin required.
obsmap = Map(dss=84, date="2020/163", project="SolarPatrol")
obsmap.initialize('sim-venus.dat', source="Venus")
logger.info("Map 'obsmap' initialized")

