|logo| 

|Release| |PyPI| |MIT licensed| |ASCL| |Travis Build| |Issues| |RTD|

Description
===========

This repo contains a set of codes to calculate correlation functions and 
other clustering statistics in a cosmological box (co-moving XYZ)
or on a mock (RA, DEC, CZ). Read the documentation on `corrfunc.rtfd.io <http://corrfunc.rtfd.io/>`_. 

Why Should You Use it
======================

1. **Fast** Theory pair-counting is **7x** faster than ``SciPy cKDTree``, and at least **2x** faster than all existing public codes.
2. **OpenMP Parallel** All pair-counting codes can be done in parallel (with strong scaling efficiency >~ 95% up to 10 cores)
3. **Python Extensions** Python extensions allow you to do the compute-heavy bits using C while retaining all of the user-friendliness of python. 
4. **Weights** All correlation functions now support *arbitrary, user-specified* weights for individual points
5. **Modular** The code is written in a modular fashion and is easily extensible to compute arbitrary clustering statistics. 
6. **Future-proof** As I get access to newer instruction-sets, the codes will get updated to use the latest and greatest CPU features. 

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
   ``theory.options`` and ``mocks.options``. You might need to ask your
   sys-admin for system-wide installs of the compiler; if you prefer to
   install your own then ``conda install gcc`` (MAC/linux) or
   ``(sudo) port install gcc5`` (on MAC) should work. 
3. ``gsl >= 2.4``. Use either
   ``conda install -c conda-forge gsl``
   (MAC/linux) or ``(sudo) port install gsl`` (MAC) to install ``gsl``
   if necessary.
4. ``python >= 2.7`` or ``python>=3.4`` for compiling the C extensions.
5. ``numpy>=1.7`` for compiling the C extensions.

Preferred Install Method
-------------------------

::

    $ git clone https://github.com/manodeep/Corrfunc/
    $ make 
    $ make install
    $ python setup.py install (--user)
    $ make tests 

Assuming you have ``gcc`` in your ``PATH``, ``make`` and
``make install`` should compile and install the C libraries + python
extensions within the source directory. If you would like to install the
python C extensions in your environment, then
``python setup.py install (--user)`` should be sufficient. If you are primarily
interested in the ``python`` interface, you can condense all of the steps
by using ``python setup.py install CC=yourcompiler (--user)`` after ``git clone``.

Compilation Notes
------------------

- If python and/or numpy are not available, then the C extensions will not be compiled.

- ``make install`` simply copies files into the ``lib/bin/include`` sub-directories. You do not need ``root`` permissions

- Default compiler on MAC is set to ``clang``, if you want to specify a different compiler, you will have to call ``make CC=yourcompiler``,  ``make install CC=yourcompiler``, ``make tests CC=yourcompiler`` etc. If you want to permanently change the default compiler, then please edit the `common.mk <common.mk>`__ file in the base directory.

- If you are directly using ``python setup.py install CC=yourcompiler (--user)``, please run a ``make distclean`` beforehand (especially if switching compilers)


Alternate Install Method
-------------------------

The python package is directly installable via ``pip install Corrfunc``. However, in that case you will lose the ability to recompile the code according to your needs. Installing via ``pip`` is **not** recommended, please open an install issue on this repo first; doing so helps improve the code-base and saves future users from running into similar install issues. 

Installation notes
------------------

If compilation went smoothly, please run ``make tests`` to ensure the
code is working correctly. Depending on the hardware and compilation
options, the tests might take more than a few minutes. *Note that the
tests are exhaustive and not traditional unit tests*.

While I have tried to ensure that the package compiles and runs out of
the box, cross-platform compatibility turns out to be incredibly hard.
If you run into any issues during compilation and you have all of the
pre-requisites, please see the `FAQ <FAQ>`__ or `email
the Corrfunc mailing list <mailto:corrfunc@googlegroups.com>`__. Also, feel free to create a new issue
with the ``Installation`` label.

Clustering Measures on a Cosmological box
-----------------------------------------

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

Clustering measures on a Mock
-----------------------------

All codes that work on mock catalogs (RA, DEC, CZ) are located in the
``mocks`` directory. The various clustering measures are:

1. ``DDrppi_mocks`` -- The standard auto/cross correlation between two data
   sets. The outputs, DD, DR and RR can be combined using ``wprp`` to
   produce the Landy-Szalay estimator for `wp(rp)`.

2. ``DDsmu_mocks`` -- The standard auto/cross correlation between two data
   sets. The outputs, DD, DR and RR can be combined using the python utility 
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
runtime if you use the python interface or the static libraries. Each one of
the following ``Makefile`` option has a corresponding entry for the runtime
libraries. 

Theory (in `theory.options <theory.options>`__)
-------------------------------------------------

1. ``PERIODIC`` (ignored in case of wp/xi) -- switches periodic boundary
   conditions on/off. Enabled by default.

2. ``OUTPUT_RPAVG`` -- switches on output of ``<rp>`` in each ``rp``
   bin. Can be a massive performance hit (~ 2.2x in case of wp).
   Disabled by default. 

3. ``DOUBLE_PREC`` -- switches on calculations in double precision. Disabled
   by default (i.e., calculations are performed in single precision by default).
   
Mocks (in `mocks.options <mocks.options>`__)
----------------------------------------------

1. ``OUTPUT_RPAVG`` -- switches on output of ``<rp>`` in each ``rp``
   bin for ``DDrppi_mocks``. Enabled by default.

2. ``OUTPUT_THETAAVG`` -- switches on output of in each theta bin. Can
   be extremely slow (~5x) depending on compiler, and CPU capabilities.
   Disabled by default.

3. ``DOUBLE_PREC`` -- switches on calculations in double precision. Disabled
   by default (i.e., calculations are performed in single precision by default).
   
4. ``LINK_IN_DEC`` -- creates binning in declination for ``DDtheta``. Please
   check that for your desired limits ``\theta``, this binning does not 
   produce incorrect results (due to numerical precision). Generally speaking,
   if your ``\thetamax`` (the max. ``\theta`` to consider pairs within) is too
   small (probaly less than 1 degree), then you should check with and without
   this option. Errors are typically sub-percent level. 

5. ``LINK_IN_RA`` -- creates binning in RA once binning in DEC has been
   enabled. Same numerical issues as ``LINK_IN_DEC``

6. ``FAST_DIVIDE`` -- Disabled by default. Divisions are slow but required
   ``DD(r_p,\pi)``. Enabling this option, replaces
   the divisions with a reciprocal followed by a Newton-Raphson. The code
   will run ~20% faster at the expense of some numerical precision.
   Please check that the loss of precision is not important for your
   use-case. 

7. ``FAST_ACOS`` -- Relevant only when ``OUTPUT_THETAAVG`` is enabled. Disabled 
   by default. An ``arccos`` is required to calculate ``<\theta>``. In absence of vectorized
   ``arccos`` (intel compiler, ``icc`` provides one via intel Short Vector Math 
   Library), this calculation is extremely slow. However, we can approximate
   ``arccos`` using polynomials (with `Remez Algorithm <https://en.wikipedia.org/wiki/Remez_algorithm>`_).
   The approximations are taken from implementations released by `Geometric Tools <http://geometrictools.com/>`_.
   Depending on the level of accuracy desired, this implementation of ``fast acos`` 
   can be tweaked in the file `utils/fast_acos.h <utils/fast_acos.h>`__. An alternate, less
   accurate implementation is already present in that file. Please check that the loss of 
   precision is not important for your use-case. 

8. ``COMOVING_DIST`` -- Currently there is no support in ``Corrfunc`` for different cosmologies. However, for the
   mocks routines like, ``DDrppi_mocks`` and ``vpf_mocks``, cosmology parameters are required to convert between
   redshift and co-moving distance. Both ``DDrppi_mocks`` and ``vpf_mocks`` expects to receive a ``redshift`` array 
   as input; however, with this option enabled, the ``redshift`` array will be assumed to contain already converted
   co-moving distances. So, if you have redshifts and want to use an arbitrary cosmology, then convert the redshifts
   into co-moving distances, enable this option, and pass the co-moving distance array into the routines. 

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
use the C extensions directly. Here are a few examples:

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
                                                        

Common Code options for both Mocks and Cosmological Boxes
=========================================================

1. ``USE_OMP`` -- uses OpenMP parallelization. Scaling is great for DD
   (close to perfect scaling up to 12 threads in our tests) and okay (runtime
   becomes constant ~6-8 threads in our tests) for ``DDrppi`` and ``wp``.
   Enabled by default. The ``Makefile`` will compare the `CC` variable with
   known OpenMP enabled compilers and set compile options accordingly. 
   Set in `common.mk <common.mk>`__ by default. 

*Optimization for your architecture*

1. The values of ``bin_refine_factor`` and/or ``zbin_refine_factor`` in
   the ``countpairs\_\*.c`` files control the cache-misses, and
   consequently, the runtime. In my trial-and-error methods, I have seen
   any values larger than 3 are always slower. But some different
   combination of 1/2 for ``(z)bin_refine_factor`` might be faster on
   your platform.

2. If you have AVX2/AVX-512/KNC, you will need to add a new kernel within
   the ``*_kernels.c`` and edit the runtime dispatch code to call this new
   kernel. 

Author & Maintainers 
=====================

Corrfunc was designed by Manodeep Sinha and is currently maintained by
`Lehman Garrison <https://github.com/lgarrison>`_ and `Manodeep Sinha <https://github.com/manodeep>`_

Citing
======

If you use the code, please cite using the `ascl entry <http://ascl.net/1703.003>`_ as indexed by `ADS <http://adsabs.harvard.edu/abs/2017ascl.soft03003S>`_. The BibTex entry for the code is

::

      @misc{2017ascl.soft03003S,
         author = {{Sinha}, M. and {Garrison}, L.},
         title = "{Corrfunc: Blazing fast correlation functions on the CPU}",
         keywords = {Software},
         howpublished = {Astrophysics Source Code Library},
         year = 2017,
         month = mar,
         archivePrefix = "ascl",
         eprint = {1703.003},
         adsurl = {http://adsabs.harvard.edu/abs/2017ascl.soft03003S},
         adsnote = {Provided by the SAO/NASA Astrophysics Data System}
      }


Mailing list
============

If you have questions or comments about the package, please do so on the
mailing list: https://groups.google.com/forum/#!forum/corrfunc

LICENSE
=======

Corrfunc is released under the MIT license. Basically, do what you want
with the code, including using it in commercial application.

Project URL
===========

-  documentation (http://corrfunc.rtfd.io/)   
-  version control (https://github.com/manodeep/Corrfunc)

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
.. |ASCL| image:: https://img.shields.io/badge/ascl-1703.003-blue.svg?colorB=262255
   :target: http://ascl.net/1703.003
   :alt: ascl:1703.003
.. |Travis Build| image:: https://travis-ci.org/manodeep/Corrfunc.svg?branch=master
   :target: https://travis-ci.org/manodeep/Corrfunc
   :alt: Build Status
.. |Issues| image:: https://img.shields.io/github/issues/manodeep/Corrfunc.svg
   :target: https://github.com/manodeep/Corrfunc/issues
   :alt: Open Issues
.. |RTD| image:: https://readthedocs.org/projects/corrfunc/badge/?version=master
   :target: http://corrfunc.readthedocs.io/en/master/?badge=master
   :alt: Documentation Status
