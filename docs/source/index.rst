.. Data_Reduction documentation master file, created by
   sphinx-quickstart on Sat Jun 5 14:57:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Single Dish Radio Astronomy Software Tools
==========================================

For an overview of SDRAST and the current status please visit https://sdrast.github.io/.

Note
====
These packages are in a very preliminary state and need lots of work.  Some code
is obsolete or buggy from not having been used and therefore not maintained a 
long time.

Base Class Module
=================
.. automodapi:: Data_Reduction
.. automodapi:: Data_Reduction.plotting

Context Modules
===============
.. automodapi:: Data_Reduction.DSN
.. automodapi:: Data_Reduction.DSN.database
.. automodapi:: Data_Reduction.DSN.logs
.. automodapi:: Data_Reduction.DSN.old_VSR
.. automodapi:: Data_Reduction.DSN.OLSR
.. automodapi:: Data_Reduction.DSN.STATS
.. automodapi:: Data_Reduction.DSN.SAO
.. automodapi:: Data_Reduction.DSN.Tid_ASCII_data
.. automodapi:: Data_Reduction.DSN.Tid_data
.. automodapi:: Data_Reduction.DSN.tipping
.. automodapi:: Data_Reduction.FITS
.. automodapi:: Data_Reduction.GAVRT
.. automodapi:: Data_Reduction.GAVRT.mysql
.. automodapi:: Data_Reduction.GAVRT.plotter
.. automodapi:: Data_Reduction.Malargue
.. automodapi:: Data_Reduction.OldGAVRT
.. automodapi:: Data_Reduction.OldGAVRT.plotting
.. automodapi:: Data_Reduction.OldGAVRT.solar

Supporting Modules
==================
.. automodapi:: Data_Reduction.boresight_fitter
.. automodapi:: Data_Reduction.hyperfine_fit
.. automodapi:: Data_Reduction.maps
.. automodapi:: Data_Reduction.SLATool
.. automodapi:: Data_Reduction.SLAPlotter
.. automodapi:: Data_Reduction.TAMS_tipping
.. automodapi:: Data_Reduction.tipping


.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
