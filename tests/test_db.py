import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

from Data_Reduction.GAVRT import Observation, Session

ses = Session(None, year=2019, doy=242)
ses.list_maps()

