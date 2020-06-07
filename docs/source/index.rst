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
is obsolete or buggy from not having been maintained a long time.  Some modules
cannot be processed by `automodapi` right now.  Functions should not be imported
from other modules; rather the modules should be imported and the extended function
name used, *e.g.* `argmax()`.

.. automodapi:: Data_Reduction
.. automodapi:: Data_Reduction.DSN
Data_Reduction.DSN.GAVRT
------------------------
cannot be handled by ``automodapi`` right now.  Needs fixing.
.. automodapi:: Data_Reduction.DSN.SAO
.. automodapi:: Data_Reduction.FITS
.. automodapi:: Data_Reduction.boresight_fitter
.. automodapi:: Data_Reduction.hyperfine_fit
Data_Reduction.maps
-------------------
cannot be handled by ``automodapi`` right now.  Needs fixing.
.. automodapi:: Data_Reduction.SLATool
Data_Reduction.SLAPlotter
-------------------------
cannot be handled by ``automodapi`` right now.  Needs fixing.

Data_Reduction.TAMS_tipping
---------------------------
cannot be handled by ``automodapi`` right now.  Needs fixing.
.. automodapi:: Data_Reduction.tipping


.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
