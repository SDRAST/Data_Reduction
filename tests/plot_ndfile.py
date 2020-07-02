"""
Tests the base MapPlotter

Gets data from a text file
"""
import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)

from Data_Reduction.plotting import MapPlotter

mapven = MapPlotter(dss=84, date="2020/163", project="SolarPatrol")
mapven.initialize('sim-venus.dat', source="Venus")
logger.info("MapPlotter 'mapven' is now initialized")

mapsun = MapPlotter(dss=84, date="2020/163", project="SolarPatrol")
mapsun.initialize('sim-sol.dat', source="Sun")
logger.info("MapPlotter 'mapsun' is now initialized")
