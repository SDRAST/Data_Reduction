"""
Data_Reduction.DSN
==================

Subclasses for reducing data taken with DSN-like open loop recorders.

Open-loop recorders are raw IF voltage recorders that are not synchronized with
the communications between the spacecraft and the ground station.  As such, they
are the most basic kind of recorder possible in radio astronomy, equivalent to
VLBI recorders.  Indeed, an early implementation was known as the "VLBI Science
Recorder" (VSR), followed by later varieties of VSR and eventually, the OSR.

OLR recordings at different stations are indeed combined for VLBI measurements
of spacecraft with respect to distant radio sources, a powerful navigation tool.

Raw IF recordings can be computational converted into any of the standard
signal types used in radio astronomy -- square-law detected power, spectra, 
Stokes parameters, VLBI U-V maps, high time ans spectral resolution pulsar data,
*etc.*
"""
import logging

import Data_Reduction as DR

logger =  logging.getLogger(__name__)

class Observation(DR.Observation):
  """
  class for observations based on open-loop recordings
  """
  def __init__(self):
    """
    """
    pass

class Map(Observation):
  """
  """
  def __init__(self):
    """
    """
    pass
