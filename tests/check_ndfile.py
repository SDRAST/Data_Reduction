"""
check structure of a text file

This test uses a t-file but any columnar text data file will work.
This tests the base Observation class.
"""

from Data_Reduction import Map, Observation

import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)

obs = Observation(dss=28, date="2012/127", project="SolarPatrol")
obs.open_datafile('t12127.10', 
                          delimiter=[17,16,3,11,7,9,8,2,6], skip_header=1,
                      names="UTC Epoch Chan Tsys Int Az El Diode Level".split())
logger.info("Observation 'obs' opened with column widths specified")

obs2 = Observation(dss=28, date="2012/127", project="SolarPatrol")
obs2.open_datafile('t12127.10', skip_header=1,
             names="Year DOY UTC Epoch Chan Tsys Int Az El Diode Level".split())
logger.info("Observation 'obs2' opened with assigned column names")

obsmap = Map(dss=84, date="2020/163", project="SolarPatrol")
obsmap.initialize('sim-venus.dat', source="Venus")
#data = obsmap.open_datafile('sim-venus.dat')
#metadata, signals = obsmap.get_data_channels(data)
#obsmap.make_data_struct(data, metadata, signals)
#obsmap.get_offsets(source="Venus")
#obsmap.make_channels(signals)
logger.info("Map 'obsmap' initialized")

