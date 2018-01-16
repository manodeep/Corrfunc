.. _which_corrfunc:

***********************************
Which correlation function to use? 
***********************************
Corrfunc has a variety of correlation functions to cover a broad range of Science applications. The basic distinction occurs if the input particles are directly
from a simulation or from an observational survey (or equivalently, a simulation that has been processed to look like a survey). For simulation data, referred throughout
as `theory`, the assumption is that the particle positions are Cartesian, co-moving XYZ.  For survey data, referred throughout as `mocks`, the assumption is that
particle positions are `Right Ascension` (0 -- 360 deg), `Declination` (-90 -- 90 deg) and `CZ` (speed of light multiplied by the redshift). Depending on the exact
type of data, **and** the desired correlation function you want, the following table should help you figure out which code you should use.


+-------------------+---------------+-----------------+-----------------------------------------+-------------------------------+---------------------------------------+
| Input Data        | Periodic      | Particle domain |    Desired correlation function         |   Returns                     | Python code                           |
+===================+===============+=================+=========================================+===============================+=======================================+
| X, Y, Z           | True          | Cube (box)      | wp(:math:`r_p`)                         | 2-D Projected Correlation     |:py:mod:`Corrfunc.theory.wp`           |
|                   |               |                 +-----------------------------------------+-------------------------------+---------------------------------------+
|                   |               |                 | :math:`\xi(r)`                          | 3-D Real-space Correlation    |:py:mod:`Corrfunc.theory.xi`           |
+-------------------+---------------+-----------------+-----------------------------------------+-------------------------------+---------------------------------------+
| X, Y, Z           | True or False | Arbitrary       | :math:`\xi(r)`                          | Pair-counts in 3-D real-space |:py:mod:`Corrfunc.theory.DD`           |
|                   |               |                 +-----------------------------------------+-------------------------------+---------------------------------------+
|                   |               |                 | :math:`\xi(r_p, \pi)`                   | Pair-counts in 2-D            |:py:mod:`Corrfunc.theory.DDrppi`       |
|                   |               |                 +-----------------------------------------+-------------------------------+---------------------------------------+
|                   |               |                 | :math:`\xi(s, \mu)`                     | Pair-counts in 2-D            |:py:mod:`Corrfunc.theory.DDsmu`        |
+-------------------+---------------+-----------------+-----------------------------------------+-------------------------------+---------------------------------------+
| ra, dec, cz       | False         | Arbitrary       | :math:`\xi(r_p, \pi)`                   | Pair-counts in 2-D            |:py:mod:`Corrfunc.mocks.DDrppi_mocks`  |
|                   |               |                 +-----------------------------------------+-------------------------------+---------------------------------------+
|                   |               |                 | :math:`\xi(s, \mu)`                     | Pair-counts in 2-D            |:py:mod:`Corrfunc.mocks.DDsmu_mocks`   |
+-------------------+---------------+-----------------+-----------------------------------------+-------------------------------+---------------------------------------+
| ra, dec           | False         | Arbitrary       | :math:`\omega(\theta)`                  | Pair-counts in angular space  |:py:mod:`Corrfunc.mocks.DDtheta_mocks` |
+-------------------+---------------+-----------------+-----------------------------------------+-------------------------------+---------------------------------------+

In all cases where only pair-counts are returned (e.g., all of the `mocks` routines), you will need to compute at least
an additional `RR` term. Please see :py:mod:`Corrfunc.utils.convert_3d_counts_to_cf` to
convert 3-D pair-counts (or angular pair counts) into a correlation
function. For 2-D pair-counts, please use :py:mod:`Corrfunc.utils.convert_rp_pi_counts_to_wp`
to convert into a projected correlation function. If you want to compute
the :math:`\xi(r_p, \pi)` from the 2-D pair-counts, then simply call
:py:mod:`Corrfunc.utils.convert_3d_counts_to_cf` with the arrays.

Also, see :ref:`commandline-interface` for a detailed list of the clustering statistics and the various available API interfaces.
    
    
