|logo|

|Release| |PyPI| |MIT licensed| |Travis Build| |GitHub CI| |RTD| |Issues|

|CoreInfra| |FAIRSoft|

|Paper I| |Paper II|


Description
===========

This repo contains a suite of codes to calculate correlation functions and
other clustering statistics for **simulated** galaxies in a cosmological box (co-moving XYZ)
and on **observed** galaxies with on-sky positions (RA, DEC, CZ). Read the
documentation on `corrfunc.rtfd.io <http://corrfunc.rtfd.io/>`_.

Why Should You Use it
======================

1. **Fast** Theory pair-counting is **7x** faster than ``SciPy cKDTree``, and at least **2x** faster than all existing public codes.
2. **OpenMP Parallel** All pair-counting codes can be done in parallel (with strong scaling efficiency >~ 95% up to 10 cores)
3. **Python Extensions** Python extensions allow you to do the compute-heavy bits using C while retaining all of the user-friendliness of Python.
4. **Weights** All correlation functions now support *arbitrary, user-specified* weights for individual points
5. **Modular** The code is written in a modular fashion and is easily extensible to compute arbitrary clustering statistics.
6. **Future-proof** As we get access to newer instruction-sets, the codes will get updated to use the latest and greatest CPU features.

*If you use the codes for your analysis, please star this repo -- that helps us keep track of the number of users.*

Benchmark against Existing Codes
================================

Please see this
`gist <https://gist.github.com/manodeep/cffd9a5d77510e43ccf0>`__ for
some benchmarks with current codes. If you have a pair-counter that you would like to compare, please add in a corresponding function and update the timings.

Installation
============

Pre-requisites
--------------

1. ``make >= 3.80``
2. OpenMP capable compiler like ``icc``, ``gcc>=4.6`` or ``clang >= 3.7``. If
   not available, please disable ``USE_OMP`` option option in
   ``theory.options`` and ``mocks.options``. On a HPC cluster, consult the cluster
   documentation for how to load a compiler (often ``module load gcc`` or similar).
   If you are using Corrfunc with Anaconda Python, then ``conda install gcc`` (MAC/linux)
   should work.  On MAC, ``(sudo) port install gcc5`` is also an option.
3. ``gsl >= 2.4``.  On an HPC cluster, consult the cluster documentation
   (often ``module load gsl`` will work).  With Anaconda Python, use
   ``conda install -c conda-forge gsl`` (MAC/linux).  On MAC, you can use
   ``(sudo) port install gsl`` (MAC) if necessary.
4. ``python >= 2.7`` or ``python>=3.4`` for compiling the CPython extensions.
5. ``numpy>=1.7`` for compiling the CPython extensions.

Method 1: Source Installation (Recommended)
-------------------------------------------

::

    $ git clone https://github.com/manodeep/Corrfunc.git
    $ cd Corrfunc
    $ make
    $ make install
    $ python -m pip install . [--user]
    
    $ make tests  # run the C tests
    $ python -m pip install pytest
    $ python -m pytest  # run the Python tests

Assuming you have ``gcc`` in your ``PATH``, ``make`` and
``make install`` should compile and install the C libraries + Python
extensions within the source directory. If you would like to install the
CPython extensions in your environment, then
``python -m pip install . [--user]`` should be sufficient. If you are primarily
interested in the Python interface, you can condense all of the steps
by using ``python -m pip install . [--user] --install-option="CC=yourcompiler"``
after ``git clone [...]`` and ``cd Corrfunc``.

Compilation Notes
~~~~~~~~~~~~~~~~~

- If Python and/or numpy are not available, then the CPython extensions will not be compiled.

- ``make install`` simply copies files into the ``lib/bin/include`` sub-directories. You do not need ``root`` permissions

- Default compiler on MAC is set to ``clang``, if you want to specify a different compiler, you will have to call ``make CC=yourcompiler``,  ``make install CC=yourcompiler``, ``make tests CC=yourcompiler`` etc. If you want to permanently change the default compiler, then please edit the `common.mk <common.mk>`__ file in the base directory.

- If you are directly using ``python -m pip install . [--user] --install-option="CC=yourcompiler"``, please run a ``make distclean`` beforehand (especially if switching compilers)

- Please note that Corrfunc is compiling with optimizations for the architecture
  it is compiled on.  That is, it uses ``gcc -march=native`` or similar.
  For this reason, please try to compile Corrfunc on the architecture it will
  be run on (usually this is only a concern in heterogeneous compute environments,
  like an HPC cluster with multiple node types).  In many cases, you can
  compile on a more capable architecture (e.g. with AVX-512 support) then
  run on a less capable architecture (e.g. with only AVX2), because the
  runtime dispatch will select the appropriate kernel.  But the non-kernel
  elements of Corrfunc may emit AVX-512 instructions due to ``-march=native``.
  If an ``Illegal instruction`` error occurs, then you'll need to recompile
  on the target architecture.

Installation notes
~~~~~~~~~~~~~~~~~~

If compilation went smoothly, please run ``make tests`` to ensure the
code is working correctly. Depending on the hardware and compilation
options, the tests might take more than a few minutes. *Note that the
tests are exhaustive and not traditional unit tests*.

For Python tests, please run ``python -m pip install pytest`` and ``python -m pytest``
from the Corrfunc root dir.

While we have tried to ensure that the package compiles and runs out of
the box, cross-platform compatibility turns out to be incredibly hard.
If you run into any issues during compilation and you have all of the
pre-requisites, please see the `FAQ <FAQ>`__ or `email
the Corrfunc mailing list <mailto:corrfunc@googlegroups.com>`__. Also, feel free to create a new issue
with the ``Installation`` label.


Method 2: pip installation
--------------------------

The Python package is directly installable via ``python -m pip install Corrfunc``. However, in that case you will lose the ability to recompile the code.  This usually fine if you are only using the Python interface and are on a single machine, like a laptop.  For usage on a cluster or other environment with multiple CPU architectures, you may find it more useful to use the Source Installation method above in case you need to compile for a different architecture later.

Testing a pip-installed Corrfunc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can check that a pip-installed Corrfunc is working with:

::

   $ python -m pytest --pyargs Corrfunc


The pip installation does not include all of the test data contained in the main repo,
since it would total over 100 MB and the tests that generate on-the-fly data are similarly
exhaustive.  pytest will mark tests where the data files are not availabe as "skipped".
If you would like to run the data-based tests, please use the Source Installation method.


OpenMP on OSX
--------------

Automatically detecting OpenMP support from the compiler and the runtime is a
bit tricky. If you run into any issues compiling (or running) with OpenMP,
please refer to the `FAQ <FAQ>`__ for potential solutions.


Clustering Measures on simulated galaxies
=========================================

Input data
----------

The input galaxies (or any discrete distribution of points) are derived from a
simulation. For instance, the galaxies could be a result of an Halo Occupation
Distribution (HOD) model, a Subhalo Abundance matching (SHAM) model, a
Semi-Empirical model (SEM), or a Semi-Analytic model (SAM) etc. The input set of
points can also be the dark matter halos, or the dark matter particles from
a cosmological simulation. The input set of points are expected to have
positions specified in Cartesian XYZ.

Types of available clustering statistics
----------------------------------------

All codes that work on cosmological boxes with co-moving positions are
located in the ``theory`` directory. The various clustering measures
are:

1. ``DD`` -- Measures auto/cross-correlations between two boxes.
   The boxes do not need to be cubes.

2. ``xi`` -- Measures 3-d auto-correlation in a cubic cosmological box.
   Assumes PERIODIC boundary conditions.

3. ``wp`` -- Measures auto 2-d point projected correlation function in a
   cubic cosmological box. Assumes PERIODIC boundary conditions.

4. ``DDrppi`` -- Measures the auto/cross correlation function between
   two boxes. The boxes do not need to be cubes.

5. ``DDsmu`` -- Measures the auto/cross correlation function between
   two boxes. The boxes do not need to be cubes.

6. ``vpf`` -- Measures the void probability function + counts-in-cells.

Clustering measures on observed galaxies
========================================

Input data
----------

The input galaxies are typically observed galaxies coming from a large-scale
galaxy survey. In addition, simulated galaxies that have been projected onto the sky
(i.e., where observational systematics have been incorporated and on-sky
positions have been generated) can also be used. We generically refer to both
these kinds of galaxies as "mocks".


The input galaxies are expected to have positions specified in spherical
co-ordinates with at least right ascension (RA) and declination (DEC).
For spatial correlation functions, an approximate "co-moving" distance
(speed of light multiplied by redshift, CZ) is also required.


Types of available clustering statistics
----------------------------------------

All codes that work on mock catalogs (RA, DEC, CZ) are located in the
``mocks`` directory. The various clustering measures are:

1. ``DDrppi_mocks`` -- The standard auto/cross correlation between two data
   sets. The outputs, DD, DR and RR can be combined using ``wprp`` to
   produce the Landy-Szalay estimator for `wp(rp)`.

2. ``DDsmu_mocks`` -- The standard auto/cross correlation between two data
   sets. The outputs, DD, DR and RR can be combined using the Python utility
   ``convert_3d_counts_to_cf`` to produce the Landy-Szalay estimator for `xi(s, mu)`.

3. ``DDtheta_mocks`` -- Computes angular correlation function between two data
   sets. The outputs from ``DDtheta_mocks`` need to be combined with
   ``wtheta`` to get the full `\omega(\theta)`

4. ``vpf_mocks`` -- Computes the void probability function on mocks.

Science options
===============

If you plan to use the command-line, then you will have to specify the
code runtime options at compile-time. For theory routines, these options
are in the file `theory.options <theory.options>`__ while for the mocks, these options are
in file `mocks.options <mocks.options>`__.

**Note** All options can be specified at
runtime if you use the Python interface or the static libraries. Each one of
the following ``Makefile`` option has a corresponding entry for the runtime
libraries.

Theory (in `theory.options <theory.options>`__)
-------------------------------------------------

1. ``PERIODIC`` (ignored in case of wp/xi) -- switches periodic boundary
   conditions on/off. Enabled by default.

2. ``OUTPUT_RPAVG`` -- switches on output of ``<rp>`` in each ``rp``
   bin. Can be a massive performance hit (~ 2.2x in case of wp).
   Disabled by default.

Mocks (in `mocks.options <mocks.options>`__)
----------------------------------------------

1. ``OUTPUT_RPAVG`` -- switches on output of ``<rp>`` in each ``rp``
   bin for ``DDrppi_mocks``. Enabled by default.

2. ``OUTPUT_THETAAVG`` -- switches on output of in each theta bin. Can
   be extremely slow (~5x) depending on compiler, and CPU capabilities.
   Disabled by default.

3. ``LINK_IN_DEC`` -- creates binning in declination for ``DDtheta_mocks``. Please
   check that for your desired limits ``\theta``, this binning does not
   produce incorrect results (due to numerical precision). Generally speaking,
   if your ``\thetamax`` (the max. ``\theta`` to consider pairs within) is too
   small (probaly less than 1 degree), then you should check with and without
   this option. Errors are typically sub-percent level.

4. ``LINK_IN_RA`` -- creates binning in RA once binning in DEC has been
   enabled for ``DDtheta_mocks``. Same numerical issues as ``LINK_IN_DEC``

5. ``FAST_ACOS`` -- Relevant only when ``OUTPUT_THETAAVG`` is enabled for
   ``DDtheta_mocks``. Disabled by default. An ``arccos`` is required to
   calculate ``<\theta>``. In absence of vectorized ``arccos`` (intel compiler,
   ``icc`` provides one via intel Short Vector Math Library), this calculation is extremely slow. However, we can approximate
   ``arccos`` using polynomials (with `Remez Algorithm <https://en.wikipedia.org/wiki/Remez_algorithm>`_).
   The approximations are taken from implementations released by `Geometric Tools <http://geometrictools.com/>`_.
   Depending on the level of accuracy desired, this implementation of ``fast acos``
   can be tweaked in the file `utils/fast_acos.h <utils/fast_acos.h>`__. An alternate, less
   accurate implementation is already present in that file. Please check that the loss of
   precision is not important for your use-case.

6. ``COMOVING_DIST`` -- Currently there is no support in ``Corrfunc`` for different cosmologies. However, for the
   mocks routines like, ``DDrppi_mocks`` and ``vpf_mocks``, cosmology parameters are required to convert between
   redshift and co-moving distance. Both ``DDrppi_mocks`` and ``vpf_mocks`` expects to receive a ``redshift`` array
   as input; however, with this option enabled, the ``redshift`` array will be assumed to contain already converted
   co-moving distances. So, if you have redshifts and want to use an arbitrary cosmology, then convert the redshifts
   into co-moving distances, enable this option, and pass the co-moving distance array into the routines.

Common Code options for both Mocks and Theory
==============================================

1. ``DOUBLE_PREC`` -- switches on calculations in double
   precision. Calculations are performed in double precision when enabled. This
   option is disabled by default in theory and enabled by default in the mocks
   routines.

2. ``USE_OMP`` -- uses OpenMP parallelization. Scaling is great for DD
   (close to perfect scaling up to 12 threads in our tests) and okay (runtime
   becomes constant ~6-8 threads in our tests) for ``DDrppi`` and ``wp``.
   Enabled by default. The ``Makefile`` will compare the `CC` variable with
   known OpenMP enabled compilers and set compile options accordingly.
   Set in `common.mk <common.mk>`__ by default.

3. ``ENABLE_MIN_SEP_OPT`` -- uses some further optimisations based on the
   minimum separation between pairs of cells. Enabled by default.

4. ``COPY_PARTICLES`` -- whether or not to create a copy of the particle
   positions (and weights, if supplied). Enabled by default (copies of the
   particle arrays **are** created)

5. ``FAST_DIVIDE`` -- Disabled by default. Divisions are slow but required
   ``DDrppi_mocks(r_p,\pi)``, ``DDsmu_mocks(s, \mu)`` and ``DD(s, \mu)``.
   Enabling this option, replaces the divisions with a reciprocal
   followed by a Newton-Raphson. The code will run ~20% faster at the expense
   of some numerical precision. Please check that the loss of precision is not
   important for your use-case.

*Optimization for your architecture*

1. The values of ``bin_refine_factor`` and/or ``zbin_refine_factor`` in
   the ``countpairs\_\*.c`` files control the cache-misses, and
   consequently, the runtime. In trial-and-error methods, Manodeep has seen
   any values larger than 3 are generally slower for theory routines but
   can be faster for mocks. But some different
   combination of 1/2 for ``(z)bin_refine_factor`` might be faster on
   your platform.

2. If you are using the angular correlation function and need ``thetaavg``,
   you might benefit from using the INTEL MKL library. The vectorized
   trigonometric functions provided by MKL can provide significant speedup.


Running the codes
=================

Read the documentation on `corrfunc.rtfd.io <http://corrfunc.rtfd.io/>`_.


Using the command-line interface
--------------------------------

Navigate to the correct directory. Make sure that the options, set in
either `theory.options <theory.options>`__ or `mocks.options <mocks.options>`__ in the root directory are
what you want. If not, edit those two files (and possibly
`common.mk <common.mk>`__), and recompile. Then, you can use the command-line
executables in each individual subdirectory corresponding to the
clustering measure you are interested in. For example, if you want to
compute the full 3-D correlation function, ``\xi(r)``, then run the
executable ``theory/xi/xi``. If you run executables without any arguments,
the program will output a message with all the required arguments.

Calling from C
--------------

Look under the `run_correlations.c <theory/examples/run_correlations.c>`__ and
`run_correlations_mocks.c <mocks/examples/run_correlations_mocks.c>`__ to see examples of
calling the C API directly. If you run the executables,
``run_correlations`` and ``run_correlations_mocks``, the output will
also show how to call the command-line interface for the various
clustering measures.

Calling from Python
-------------------

If all went well, the codes can be directly called from ``python``.
Please see `call_correlation_functions.py <Corrfunc/call_correlation_functions.py>`__ and
`call_correlation_functions_mocks.py <Corrfunc/call_correlation_functions_mocks.py>`__ for examples on how to
use the CPython extensions directly. Here are a few examples:

.. code:: python

    from __future__ import print_function
    import os.path as path
    import numpy as np
    import Corrfunc
    from Corrfunc.theory import wp

    # Setup the problem for wp
    boxsize = 500.0
    pimax = 40.0
    nthreads = 4

    # Create a fake data-set.
    Npts = 100000
    x = np.float32(np.random.random(Npts))
    y = np.float32(np.random.random(Npts))
    z = np.float32(np.random.random(Npts))
    x *= boxsize
    y *= boxsize
    z *= boxsize

    # Setup the bins
    rmin = 0.1
    rmax = 20.0
    nbins = 20

    # Create the bins
    rbins = np.logspace(np.log10(0.1), np.log10(rmax), nbins + 1)

    # Call wp
    wp_results = wp(boxsize, pimax, nthreads, rbins, x, y, z, verbose=True, output_rpavg=True)

    # Print the results
    print("#############################################################################")
    print("##       rmin           rmax            rpavg             wp            npairs")
    print("#############################################################################")
    print(wp_results)


Author & Maintainers
=====================

Corrfunc was designed and implemented by `Manodeep Sinha <https://github.com/manodeep>`_,
with contributions from `Lehman Garrison <https://github.com/lgarrison>`_,
`Nick Hand <https://github.com/nickhand>`_, and `Arnaud de Mattia <https://github.com/adematti>`_.
Corrfunc is currently maintained by Manodeep Sinha and Lehman Garrison.

Citing
======

If you use ``Corrfunc`` for research, please cite using the MNRAS code paper with the following
bibtex entry:

::

   @ARTICLE{2020MNRAS.491.3022S,
       author = {{Sinha}, Manodeep and {Garrison}, Lehman H.},
       title = "{CORRFUNC - a suite of blazing fast correlation functions on
       the CPU}",
       journal = {\mnras},
       keywords = {methods: numerical, galaxies: general, galaxies:
       haloes, dark matter, large-scale structure of Universe, cosmology:
       theory},
       year = "2020",
       month = "Jan",
       volume = {491},
       number = {2},
       pages = {3022-3041},
       doi = {10.1093/mnras/stz3157},
       adsurl =
       {https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.3022S},
       adsnote = {Provided by the SAO/NASA
       Astrophysics Data System}
   }


If you are using ``Corrfunc v2.3.0`` or later, **and** you benefit from the
enhanced vectorised kernels, then please additionally cite this paper:

::

      @InProceedings{10.1007/978-981-13-7729-7_1,
          author="Sinha, Manodeep and Garrison, Lehman",
          editor="Majumdar, Amit and Arora, Ritu",
          title="CORRFUNC: Blazing Fast Correlation Functions with AVX512F SIMD Intrinsics",
          booktitle="Software Challenges to Exascale Computing",
          year="2019",
          publisher="Springer Singapore",
          address="Singapore",
          pages="3--20",
          isbn="978-981-13-7729-7",
          url={https://doi.org/10.1007/978-981-13-7729-7_1}
      }



Mailing list
============

If you have questions or comments about the package, please do so on the
mailing list: https://groups.google.com/forum/#!forum/corrfunc

LICENSE
=======

Corrfunc is released under the MIT license. Basically, do what you want
with the code, including using it in commercial application.

Project URLs
============

-  Documentation (http://corrfunc.rtfd.io/)
-  Source Repository (https://github.com/manodeep/Corrfunc)
-  Entry in the Astrophysical Source Code Library (ASCL) |ASCL|
-  Zenodo Releases |Zenodo|

.. |logo| image:: https://github.com/manodeep/Corrfunc/blob/master/corrfunc_logo.png
    :target: https://github.com/manodeep/Corrfunc
    :alt: Corrfunc logo
.. |Release| image:: https://img.shields.io/github/release/manodeep/Corrfunc.svg
   :target: https://github.com/manodeep/Corrfunc/releases/latest
   :alt: Latest Release
.. |PyPI| image:: https://img.shields.io/pypi/v/Corrfunc.svg
   :target: https://pypi.python.org/pypi/Corrfunc
   :alt: PyPI Release
.. |MIT licensed| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://raw.githubusercontent.com/manodeep/Corrfunc/master/LICENSE
   :alt: MIT License
.. |Travis Build| image:: https://travis-ci.com/manodeep/Corrfunc.svg?branch=master
   :target: https://travis-ci.com/manodeep/Corrfunc
   :alt: Build Status
.. |GitHub CI| image:: https://github.com/manodeep/Corrfunc/workflows/GitHub%20CI/badge.svg
   :target: https://github.com/manodeep/Corrfunc/actions
   :alt: GitHub Actions Status
.. |Issues| image:: https://img.shields.io/github/issues/manodeep/Corrfunc.svg
   :target: https://github.com/manodeep/Corrfunc/issues
   :alt: Open Issues
.. |RTD| image:: https://readthedocs.org/projects/corrfunc/badge/?version=master
   :target: http://corrfunc.readthedocs.io/en/master/?badge=master
   :alt: Documentation Status

.. |CoreInfra| image:: https://bestpractices.coreinfrastructure.org/projects/5037/badge
   :target: https://bestpractices.coreinfrastructure.org/en/projects/5037
   :alt: Core Infrastructure Best Practices Status

.. |FAIRSoft| image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
   :target: https://fair-software.eu
   :alt: Fair Software (EU) Compliance

.. |Paper I| image:: https://img.shields.io/badge/arXiv-1911.03545-%23B31B1B
   :target: https://arxiv.org/abs/1911.03545
   :alt: Corrfunc Paper I
.. |Paper II| image:: https://img.shields.io/badge/arXiv-1911.08275-%23B31B1B
   :target: https://arxiv.org/abs/1911.08275
   :alt: Corrfunc Paper II

.. |ASCL| image:: https://img.shields.io/badge/ascl-1703.003-blue.svg?colorB=262255
   :target: http://ascl.net/1703.003
   :alt: ascl:1703.003
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3634195.svg
   :target: https://doi.org/10.5281/zenodo.3634195
