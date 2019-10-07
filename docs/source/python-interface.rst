.. _python-interface:

****************************************
Using the python extensions in Corrfunc
****************************************

This guide assumes that you already followed the :ref:`step_by_step_install`
section of the documentation to get the package and its dependencies set
up on your machine. Rest of document also assumes that you have installed
the C extensions for python.


Importing Corrfunc
===================

After installing Corrfunc you can open up a python terminal and import the
base package by:

    >>> import Corrfunc

All of the functionality is divided into ``theory`` routines and ``mocks``
routines. These routines can be independently imported by using:

    >>> from Corrfunc.theory import *
    >>> from Corrfunc.mocks import *

You can access the full API documentation by simply typing:

.. code:: python
          
          help(DD)              # theory pair-counter in 3-D separation (r)
          help(DDrppi_mocks)    # mocks pair-counter in 2-D (rp, pi)

.. _first_steps:

First steps with Corrfunc
============================

Overview of Corrfunc inputs
------------------------------

Broadly speaking, Corrfunc requires these following inputs:

* (At least) 3 arrays specifying the positions for the particles
  
  - For ``Corrfunc.theory`` routines, these positions are Cartesian XYZ in
    co-moving ``Mpc/h`` units.

  - For ``Corrfunc.mocks`` routines, these positions are ``Right Ascension``,
    ``Declination``, and ``Speed of Light * Redshift`` or ``Co-moving
    distance``. The angles are expected in degrees, while the distance is
    expected in co-moving ``Mpc/h``.

    See :ref:`read_catalog` for details on how to read in arrays from a file.

* A boolean flag specifying in an auto-correlation or cross-correlation is
  being performed. In case of cross-correlations, another set of 3 arrays
  **must** be passed as input. This second set of arrays typically represents
  randoms for ``Corrfunc.mocks``.
 
* A file containing the bins for the clustering statistic (where
  relevant). Look at ``theory/tests/bins`` for an example of the contents of
  the file for spatial bins. See ``mocks/tests/angular_bins`` for an example
  containing angular bins for mocks routines. Passing a filename is the most
  general way of specifying bins in Corrfunc. However, you can also pass in a
  1-D array for the bins.   
  
  See :ref:`generate_bins` for details on how to specify the bins as a file as
  well as an array

See :ref:`function_usage_examples` for a broad overview of the typical tasks
associated with computing correlation functions. Read on for the various
pair-counters available within the python interfaces of Corrfunc. 

Calculating spatial clustering statistics in simulation boxes
==============================================================

Corrfunc can compute a range of spatial correlation functions and the
counts-in-cells. For all of these calculations a few inputs are required. The
following code section sets up the default inputs that are used later on in the
clustering functions:


.. code:: python
               
               >>> import numpy as np
               >>> from Corrfunc.io import read_catalog
          
               # Read the default galaxies supplied with
               # Corrfunc. ~ 1 million galaxies on a 420 Mpc/h cube
               >>> X, Y, Z = read_catalog()
               
               # Specify boxsize for the XYZ arrays
               >>> boxsize = 420.0
               
               # Number of threads to use
               >>> nthreads = 2

               # Create the bins array
               >>> rmin = 0.1
               >>> rmax = 20.0
               >>> nbins = 20
               >>> rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
          
               # Specify the distance to integrate along line of sight
               >>> pimax = 40.0

               # Specify the max. of the cosine of the angle to the LOS for
               # DD(s, mu) 
               >>> mu_max = 1.0

               # Specify the number of linear bins in `mu`
               >>> nmu_bins = 20
               
               # Specify that an autocorrelation is wanted
               >>> autocorr = 1

.. testsetup:: theory

          import numpy as np
          from Corrfunc.io import read_catalog
          X, Y, Z = read_catalog()
          boxsize = 420.0
          nthreads = 2
          rmin = 0.1
          rmax = 20.0
          nbins = 20
          rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
          pimax = 40.0
          mu_max = 1.0
          nmu_bins = 20
          autocorr = 1

          
Calculating 2-D projected auto-correlation (``Corrfunc.theory.wp``)
---------------------------------------------------------------------

Corrfunc can directly compute the projected auto-correlation function,
:math:`w_p(r_p)`. This calculation sets periodic boundary conditions. Randoms
are calculated analytically based on the supplied boxsize. The projected
separation, :math:`r_p` is calculated in the X-Y plane while the line-of-sight
separation, :math:`\pi` is calculated in the Z plane. Only pairs with
:math:`\pi` separation less than :math:`\pi_{max}` are counted.

.. testcode:: theory

          from Corrfunc.theory.wp import wp
          results_wp = wp(boxsize, pimax, nthreads, rbins, X, Y, Z)
          

Calculating 3-D autocorrelation (``Corrfunc.theory.xi``)
------------------------------------------------------------

Corrfunc can also compute the 3-D auto-correlation function,
:math:`\xi(r)`. Like :math:`w_p(r_p)`, this calculation also enforces periodic
boundary conditions and an auto-correlation. Randoms are calculated
analytically on the supplied boxsize. 

.. testcode:: theory

          from Corrfunc.theory.xi import xi
          results_xi = xi(boxsize, nthreads, rbins, X, Y, Z)

   
Calculating 3-D pair-counts (``Corrfunc.theory.DD``)
-----------------------------------------------------

Corrfunc can return the pair counts in 3-D real-space for a set of arrays. The
calculation can be either auto or cross-correlation, *and* with or without periodic
boundaries. The pairs are always double-counted. Additionally, if the smallest
bin is ``0.0`` for an autocorrelation, then the self-pairs *will* be counted.

.. testcode:: theory

          from Corrfunc.theory.DD import DD
          results_DD = DD(autocorr, nthreads, rbins, X, Y, Z, boxsize=boxsize)
          

Calculating 2-D pair-counts (``Corrfunc.theory.DDrppi``)
--------------------------------------------------------

Corrfunc can return the pair counts in 2-D real-space for a set of arrays. The
calculation can be either auto or cross-correlation, *and* with or without periodic
boundaries. The projected separation, :math:`r_p` is calculated in the X-Y plane while the
line-of-sight separation, :math:`\pi` is calculated in the Z plane.

The pairs are always double-counted. Additionally, if the smallest
bin is ``0.0`` for an autocorrelation, then the self-pairs *will* be counted.

.. testcode:: theory

          from Corrfunc.theory.DDrppi import DDrppi
          results_DDrppi = DDrppi(autocorr, nthreads, pimax, rbins, X, Y, Z, boxsize=boxsize)

Calculating 2-D pair-counts (``Corrfunc.theory.DDsmu``)
--------------------------------------------------------

Corrfunc can return the pair counts in 2-D real-space for a set of arrays. The
calculation can be either auto or cross-correlation, *and* with or without periodic
boundaries. The spatial separation, :math:`s` is calculated in 3-D while 
:math:`mu` is the cosine of angle to the line-of-sight and is calculated
assuming that the Z-axis is the line-of-sight.

.. math::
   
   \mathbf{s}  &= \mathbf{v_1} - \mathbf{v_2}, \\
   {\mu} &= \frac{\left(z_1 - z_2 \right)}{\Vert\mathbf{s}\Vert}

where, :math:`\mathbf{v_1}:=(x_1, y_1, z_1)` and :math:`\mathbf{v_2}:=(x_2, y_2, z_2)` are the vectors for the
two points under consideration, and, :math:`\Vert\mathbf{s}\Vert=\sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}`

The pairs are always double-counted. Additionally, if the smallest
bin is ``0.0`` for an autocorrelation, then the self-pairs *will* be counted.

.. testcode:: theory

          from Corrfunc.theory.DDsmu import DDsmu
          results_DDsmu = DDsmu(autocorr, nthreads, rbins, mu_max, nmu_bins, X, Y, Z, boxsize=boxsize)
              

Calculating the Counts-in-Cells (``Corrfunc.theory.vpf``)
---------------------------------------------------------
Corrfunc can calculate the counts-in-cells statistics. The simplest example for
counts-in-cells is the Void Probability Function -- the probability that a
sphere of a certain size contains zero galaxies.

.. testcode:: theory

          from Corrfunc.theory.vpf import vpf

          # Maximum radius of the sphere in Mpc/h
          rmax = 10.0

          # Number of bins to cover up to rmax
          nbins = 10

          # Number of random spheres to place
          nspheres = 10000

          # Max number of galaxies in sphere (must be >=1)
          numpN = 6

          # Random number seed (used for choosing sphere centres)
          seed = 42

          results_vpf = vpf(rmax, nbins, nspheres, numpN, seed, X, Y, Z, boxsize=boxsize)


Calculating clustering statistics in mock catalogs
===================================================
In order to calculate clustering statistics in mock catalogs, the galaxy
positions are assumed to be specified as on-sky (``Right Ascension``, 
``Declination``, and ``speed of light * redshift``). The following code section
sets up the default arrays and parameters for the actual clustering calculations:


.. code:: python
          
          import numpy as np
          import Corrfunc
          from os.path import dirname, abspath, join as pjoin
          from Corrfunc.io import read_catalog

          # Mock catalog (SDSS-North) supplied with Corrfunc
          mock_catalog = pjoin(dirname(abspath(Corrfunc.__file__)), "../mocks/tests/data/", "Mr19_mock_northonly.rdcz.ff")
          RA, DEC, CZ = read_catalog(mock_catalog)

          # Randoms catalog (SDSS-North) supplied with Corrfunc
          randoms_catalog = pjoin(dirname(abspath(Corrfunc.__file__)), "../mocks/tests/data/", "Mr19_randoms_northonly.rdcz.ff")
          RAND_RA, RAND_DEC, RAND_CZ = read_catalog(randoms_catalog)
                                  
          # Number of threads to use
          nthreads = 2

          # Specify cosmology (1->LasDamas, 2->Planck)
          cosmology = 1 
          
          # Create the bins array
          rmin = 0.1
          rmax = 20.0
          nbins = 20
          rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
          
          # Specify the distance to integrate along line of sight
          pimax = 40.0

          # Specify the max. of the cosine of the angle to the LOS
          # for DD(s, mu)
          mu_max = 1.0

          # Specify the number of linear bins in `mu`
          nmu_bins = 20

          # Specify that an autocorrelation is wanted
          autocorr = 1


.. testsetup:: mocks

          import numpy as np
          import Corrfunc
          from os.path import dirname, abspath, join as pjoin
          from Corrfunc.io import read_catalog
          mock_catalog = pjoin(dirname(abspath(Corrfunc.__file__)), "../mocks/tests/data/", "Mr19_mock_northonly.rdcz.ff")
          RA, DEC, CZ = read_catalog(mock_catalog)
          randoms_catalog = pjoin(dirname(abspath(Corrfunc.__file__)), "../mocks/tests/data/", "Mr19_randoms_northonly.rdcz.ff")
          RAND_RA, RAND_DEC, RAND_CZ = read_catalog(randoms_catalog)
          nthreads = 2
          cosmology = 1 
          rmin = 0.1
          rmax = 20.0
          nbins = 20
          rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
          pimax = 40.0
          mu_max = 1.0
          nmu_bins = 20
          autocorr = 1

          
Calculating 2-D pair counts (``Corrfunc.mocks.DDrppi_mocks``)
-------------------------------------------------------------
Corrfunc can calculate pair counts for mock catalogs. The input positions are
expected to be ``Right Ascension``, ``Declination`` and ``CZ`` (speed of light
times redshift, in ``Mpc/h``). Cosmology has to be specified since ``CZ`` needs
to be converted into co-moving distance. If you want to calculate in arbitrary
cosmology, then convert ``CZ`` into co-moving distance, and then pass the
converted array while setting the option ``is_comoving_dist=True``. The
projected and line of sight separations are calculated using the following
equations from `Zehavi et al. 2002 <http://adsabs.harvard.edu/abs/2002ApJ...571..172Z>`_

.. math::
   
   \mathbf{s} &= \mathbf{v_1} - \mathbf{v_2}, \\
   \mathbf{l} &= \frac{1}{2}\left(\mathbf{v_1} + \mathbf{v_2}\right), \\
   \pi &= \left(\mathbf{s} \cdot \mathbf{l}\right)/\Vert\mathbf{l}\Vert, \\
   r_p^2 &= \mathbf{s} \cdot \mathbf{s} - \pi^2

where, :math:`\mathbf{v_1}:=(x_1, y_1, z_1)` and :math:`\mathbf{v_2}:=(x_2, y_2, z_2)` are the vectors for the
two points under consideration, and, :math:`\Vert\mathbf{s}\Vert=\sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}`.

Here is the python code to call ``Corrfunc.mocks.DDrppi_mocks``:

.. testcode:: mocks

          from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
          results_DDrppi_mocks = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, rbins, RA, DEC, CZ)
  

Calculating 2-D pair counts (``Corrfunc.mocks.DDsmu_mocks``)
-------------------------------------------------------------
Corrfunc can calculate pair counts for mock catalogs. The input positions are
expected to be ``Right Ascension``, ``Declination`` and ``CZ`` (speed of light
times redshift, in ``Mpc/h``). Cosmology has to be specified since ``CZ`` needs
to be converted into co-moving distance. If you want to calculate in arbitrary
cosmology, then convert ``CZ`` into co-moving distance, and then pass the
converted array while setting the option ``is_comoving_dist=True``. The
projected and line of sight separations are calculated using the following
equations from `Zehavi et al. 2002 <http://adsabs.harvard.edu/abs/2002ApJ...571..172Z>`_

.. math::
   
   \mathbf{s} &= \mathbf{v_1} - \mathbf{v_2}, \\
   \mathbf{l} &= \frac{1}{2}\left(\mathbf{v_1} + \mathbf{v_2}\right), \\
   \mu &= \left(\mathbf{s} \cdot \mathbf{l}\right)/\left(\Vert\mathbf{l}\Vert \Vert\mathbf{s}\Vert \right)
   
where, :math:`\mathbf{v_1}:=(x_1, y_1, z_1)` and :math:`\mathbf{v_2}:=(x_2, y_2, z_2)` are the vectors for the
two points under consideration, and, :math:`\Vert\mathbf{s}\Vert=\sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}`
   
Here is the python code to call ``Corrfunc.mocks.DDsmu_mocks``:

.. testcode:: mocks

          from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
          results_DDsmu_mocks = DDsmu_mocks(autocorr, cosmology, nthreads, mu_max, nmu_bins, rbins, RA, DEC, CZ)
  
              
Calculating angular pair-counts (``Corrfunc.mocks.DDtheta_mocks``)
-------------------------------------------------------------------
Corrfunc can compute angular pair counts for mock catalogs. The input positions
are expected to be ``Right Ascension`` and ``Declination``. Since all
calculations are in angular space, cosmology is not required.

.. testcode:: mocks

          from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
          results_DDtheta_mocks = DDtheta_mocks(autocorr, nthreads, rbins, RA, DEC)

          
Calculating the Counts-in-Cells (``Corrfunc.mocks.vpf_mocks``)
---------------------------------------------------------------
Corrfunc can calculate the counts-in-cells statistics. The simplest example for
counts-in-cells is the Void Probability Function -- the probability that a
sphere of a certain size contains zero galaxies.

.. testcode:: mocks

          from Corrfunc.mocks.vpf_mocks import vpf_mocks

          # Maximum radius of the sphere in Mpc/h
          rmax = 10.0

          # Number of bins to cover up to rmax
          nbins = 10

          # Number of random spheres to place
          nspheres = 10000

          # Max number of galaxies in sphere (must be >=1)
          numpN = 6

          # Minimum number of random points needed in a ``rmax`` sphere
          # such that it is considered to be entirely within the mock
          # footprint. Does not matter in this case, since we already
          # have the centers for the fully enclosed spheres
          threshold_ngb = 1                                                                                                        

          # File with sphere centers (centers such that spheres with size
          # rmax=10 Mpc/h are completely inside the survey)
          centers_file = pjoin(dirname(abspath(Corrfunc.__file__)), "../mocks/tests/data/", "Mr19_centers_xyz_forVPF_rmax_10Mpc.txt")

          results_vpf_mocks = vpf_mocks(rmax, nbins, nspheres, numpN, threshold_ngb, centers_file, cosmology, RA, DEC, CZ, RAND_RA, RAND_DEC, RAND_CZ)

See the complete reference here :py:mod:`Corrfunc`.
