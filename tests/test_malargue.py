"""
Test Data_Reduction.Malargue
"""
import logging
from Astronomy.DSN_coordinates import DSS
from Data_Reduction.Malargue import Map, Observation

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

obs = Observation(dss=84, date="2020/163", project="SolarPatrol",
                  datafile="sim-venus.dat")
obs.get_offsets(source="Venus")

mapobs = Map(dss=84, date="2020/163", project="SolarPatrol",
             datafile="sim-venus.dat", source='Venus')
