"""
Test Data_Reduction.Malargue
"""
import logging
from Astronomy.DSN_coordinates import DSS
from Data_Reduction.Malargue import Observation

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

obs = Observation(dss=84, date="2020/163", project="SolarPatrol",
                  datafile="sim-venus.dat")


