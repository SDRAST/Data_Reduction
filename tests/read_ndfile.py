import logging

from Data_Reduction import Observation

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

datadir = "/usr/local/projects/SolarPatrol/Observations/dss28/2012/127/"
obs = Observation(dss=28, datafile=datadir+'t12127.10', 
                  delimiter=[17,16,3,11,7,9,8,2,6], 
                  names="UTC Epoch Chan Tsys Int Az El Diode Level".split(), 
                  skip_header=1)
