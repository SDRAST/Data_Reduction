"""
check structure of a text file

This test uses a t-file but any columnar text data file will work.
This tests the base Observation class.
"""

from Data_Reduction import Observation
from Astronomy.DSN_coordinates import DSS
import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

tel = DSS(28)
obs = Observation(dss=28, date="2012/127", project="SolarPatrol")
obs.extended_init(datafile='t12127.10', 
                          delimiter=[17,16,3,11,7,9,8,2,6], skip_header=1,
                      names="UTC Epoch Chan Tsys Int Az El Diode Level".split())
                      
obs2 = Observation(dss=28, date="2012/127", project="SolarPatrol")
obs.extended_init(datafile='t12127.10', skip_header=1,
             names="Year DOY UTC Epoch Chan Tsys Int Az El Diode Level".split())


