Specification for storage of boresight data
===========================================

Dean Shaff
~~~~~~~~~~

Boresight data are complicated. A single boresight run might have at minimum
data from two difference position offset axes, running in a single direction,
from 4 different power meter channels. In order to cope with the complexity of
these data, I propose a specification for storing and serializing boresight data.
The following is a list of definitions used to describe how boresight data are
collected, and how they are stored.

* data: Offset and system temperature data.
* meta_data: Information describing when the boresight was conducted, at what
azimuth and elevation angles, and using what calibration source (if any).
* axis: Position offset axis. A single boresight run involves scanning or
stepping along two different axes. Usually this is either elevation ``el``
and cross elevation ``xel``.
* direction: Describes moving in positive or negative direction when scanning
or stepping.
    * ``right``: Moving from negative offset to positive offset
    * ``left``: Moving from position offset to negative offset
* channel: power meter channel.

Unless otherwise noted, whenever these terms are used as means of accessing an
object, they should be lowercase, using snake case. For example, when accessing
data from a Python object representation of boresight data, I might do the
following:

.. code-block:: python

    boresight_obj["data"]["el"]["right"]
    boresight_obj["meta_data"]["timestamp"]
    # the following are not acceptable:
    boresight_obj["Data"]["EL"]["Right"]
    boresight_obj["metaData"]["timeStamp"]

The main distinction between ``meta_data`` and ``data`` is that the ``data`` field
contains data that can be easily converted to binary data, while the ``meta_data``
field contains data that is more descriptive, thus making it harder to convert
to binary. For some data formats, namely JSON, this distinction isn't particularly
meaningful, while for some it is very meaningful (HDF5 and FITS).

``meta_data`` structure
---------------------

Note that I use Python data type names when describing the data types of
``meta_data`` members.

A boresight file/object should have at minimum the following members in the
``meta_data`` field:

* timestamp (str): A UTC timestamp indicating when the boresight run began. It should
be of the following format: ``YYYY-jjj-HHhMMmSSs`` OR ``YYYY-mm-ddTHH:MM:SS``.
Using Python datetime formatting convention, these formats look as follows:
``"%Y-%j-%Hh%Mm%Ss"`` and ``"%Y-%m-%dT%H:%M:%S"``, respectively.
* initial_el (float): Elevation position offset before boresight run
* initial_xel (float): cross elevation position offset before boresight run
* boresight_type (str): "scanning" or "stepping" boresight.

The following members are optional, but highly recommended, as they provide
useful descriptive information about the boresight run:

* az (float): Antenna azimuth angle when boresight run started.
* el (float): Antenna elevation angle when boresight run started.
* name (str): The name of the calibrator source used for boresight.
* flux (float): calibrator flux.
* ra_J2000 (float): RA J2000 of calibrator source, in radians.
* dec_J2000 (float): DEC J2000 of calibrator source, in radians
* xt (float): XT subreflector position.
* yt (float): YT subreflector position.
* x (float): X subreflector position.
* y (float): Y subreflector position.
* z (float): Z subreflector position.
* rate (float): Offset rate (scanning boresight only)
* sample_rate (float): power meter sample rate (scanning boresight only)
* limit (float): How far on either side of initial axis
point to start and finish (scanning boresight only)
* integration_time (float): Amount of time for which to integrate
power meter data at each boresight step (stepping boresight only)
* settle_time (float): Maximum amount of time to settle between
boresight steps. (stepping boresight only)


``data`` structure
------------------

* axis
    * direction
        * offset_data
        * channels
            * channel
                * tsys_data
                * fit (optional)

For example, I might have done a scanning boresight run in elevation and
cross elevation, in both right and left directions. My data field would look as
follows:

* el
    * right
        * channels
            * 0
                * offset_data (array of float)
                * tsys_data (array of float)
                * fit
                    * chi_sqr (float): Total chi square
                    * chi_sqr_dof (float): chi square per degree of freedom
                    * popt (array of float): optimal fit parameters
                    * pcov (2D array of float): optimal fit parameter covariance
                        matrix
                    * amplitude (float): amplitude of gaussian fit
                    * offset (float): mean of gaussian fit
                    * sigma (float): sigma of gaussina fit
                    * amplitude_err (float): Error in ``amplitude`` parameter
                    * offset_err (float): Error in ``offset`` parameter
                    * sigma_err (float): Error in ``sigma`` parameter
            * 1
                * ...
            * 2
                * ...
            * 3
                * ...
    * left
        * ...
* xel
    * ...

``fit`` data is optional, and implementation dependent.
Normally boresight data is fit with a gaussian plus a linear funciton.
In the above example, the ``popt`` and ``pcov`` arrays are the return values of
the ``scipy.optimize.curve_fit`` routine.

Implementations
---------------
