|Release| |MIT licensed| |DOI| |Travis Build| |Issues| |Coverity|

Description
===========

This repo contains a set of codes to measure the following OpenMP
parallelized clustering measures in a cosmological box (co-moving XYZ)
or on a mock (RA, DEC, CZ). Also, contains the associated paper to be
published in Astronomy & Computing Journal (at some point).

Why Should You Use it
======================

1. **Fast** All theory pair-counting is at least an order of magnitude faster than all existing public codes. Particularly suited for MCMC. 
2. **Python Extensions** Python extensions allow you to do the compute-heavy bits using C while retaining all of the user-friendliness of python. 
3. **Modular** The code is written in a modular fashion and is easily extensible to compute arbitrary clustering statistics. 
4. **Future-proof** As I get access to newer instruction-sets, the codes will get updated to use the latest and greatest CPU features. 

Installation
============

Pre-requisites
--------------

1. ``make >= 3.80``
2. OpenMP capable compiler like ``icc``, ``gcc`` or ``clang >= 3.7``. If
   not available, please disable ``USE_OMP`` option option in
   ``theory.options`` and ``mocks.options``. You might need to ask your
   sys-admin for system-wide installs of the compiler; if you prefer to
   install your own then ``conda install gcc`` (MAC/linux) or
   ``(sudo) port install gcc5`` (on MAC) should work. *Note ``gcc`` on
   macports defaults to ``gcc48`` and the portfile is currently broken
   on ``El Capitan``*.
3. ``gsl``. Use either
   ``conda install -c https://conda.anaconda.org/asmeurer gsl``
   (MAC/linux) or ``(sudo) port install gsl`` (MAC) to install ``gsl``
   if necessary.
4. ``python >= 2.6`` or ``python>=3.4`` for compiling the C extensions.
5. ``numpy>=1.7`` for compiling the C extensions.

*If python and/or numpy are not available, then the C extensions will
not be compiled*.

*Default compiler on MAC is set to ``clang``, if you want to specify a
different compiler, you will have to call ``make CC=yourcompiler``*

Preferred Method
----------------

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
``python setup.py install (--user)`` should be sufficient.

Alternative
-----------

The python package is directly installable via ``pip install Corrfunc``. However, in that case you will lose the ability to recompile the code according to your needs. Not recommended unless you are desperate (i.e., `email me <mailto:manodeep@gmail.com>`__ if you are having install issues). 

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
me <mailto:manodeep@gmail.com>`__. Also, feel free to create a new issue
with the ``Installation`` label.

Clustering Measures on a Cosmological box
-----------------------------------------

All codes that work on cosmological boxes with co-moving positions are
located in the ``xi_theory`` directory. The various clustering measures
are:

1. ``xi_of_r`` -- Measures auto/cross-correlations between two boxes.
   The boxes do not need to be cubes.

2. ``xi`` -- Measures 3-d auto-correlation in a cubic cosmological box.
   Assumes PERIODIC boundary conditions.

3. ``wp`` -- Measures auto 2-d point projected correlation function in a
   cubic cosmological box. Assumes PERIODIC boundary conditions.

4. ``xi_rp_pi`` -- Measures the auto/cross correlation function between
   two boxes. The boxes do not need to be cubes.

5. ``vpf`` -- Measures the void probability function + counts-in-cells.

Clustering measures on a Mock
-----------------------------

All codes that work on mock catalogs (RA, DEC, CZ) are located in the
``xi_mocks`` directory. The various clustering measures are:

1. ``DDrppi`` -- The standard auto/cross correlation between two data
   sets. The outputs, DD, DR and RR can be combined using ``wprp`` to
   produce the Landy-Szalay estimator for :math:`w_p(r_p)`.

2. ``wtheta`` -- Computes angular correlation function between two data
   sets. The outputs from ``DDtheta_mocks`` need to be combined with
   ``wtheta`` to get the full :math:`\omega(\theta)`

3. ``vpf`` -- Computes the void probability function on mocks.

Science options
===============

1. ``PERIODIC`` (ignored in case of wp/xi) -- switches periodic boundary
   conditions on/off. Enabled by default.

2. ``OUTPUT_RPAVG`` -- switches on output of ``<rp>`` in each ``rp``
   bin. Can be a massive performance hit (~ 2.2x in case of wp).
   Disabled by default. Needs code option ``DOUBLE_PREC`` to be enabled
   as well. For the mocks, ``OUTPUT_RPAVG`` causes only a mild increase
   in runtime and is enabled by default.

3. ``OUTPUT_THETAAVG`` -- switches on output of in each theta bin. Can
   be extremely slow (~5x) depending on compiler, and CPU capabilities.
   Disabled by default.

Mocks
-----

1. ``LINK_IN_DEC`` -- creates binning in declination for mocks. Please
   check that for your desired binning in :math:`r_p`/:math:`\theta`,
   this binning does not produce incorrect results (due to numerical
   precision).

2. ``LINK_IN_RA`` -- creates binning in RA once binning in DEC has been
   enabled. Same numerical issues as ``LINK_IN_DEC``

3. ``FAST_DIVIDE`` -- Divisions are slow but required
   :math:`DD(r_p,\pi)`. This ``Makefile`` option (in ``mocks.options``) replaces
   the divisions to a reciprocal followed by a Newton-Raphson. The code
   will run ~20% faster at the expense of some numerical precision.
   Please check that the loss of precision is not important for your
   use-case. Also, note that the mocks tests for :math:`DD(r_p, \pi)`
   *will fail* if you enable ``FAST_DIVIDE``.

Running the codes
=================

The documentation is lacking currently but I am actively working on it.

Using the command-line interface
--------------------------------

Navigate to the correct directory. Make sure that the options, set in
either ``theory.options`` or ``mocks.options`` in the root directory are
what you want. If not, edit those two files (and possibly
``common.mk``), and recompile. Then, you can use the command-line
executables in each individual subdirectory corresponding to the
clustering measure you are interested in. For example, if you want to
compute the full 3-D correlation function, ``\xi(r)``, then navigate to
``xi_theory/xi`` and run the executable ``xi``. If you run executables
without any arguments, the message will you tell you all the required
arguments.

Calling from C
--------------

Look under the ``xi_theory/examples/run_correlations.c`` and
``xi_mocks/examples/run_correlations_mocks.c`` to see examples of
calling the C API directly. If you run the executables,
``run_correlations`` and ``run_correlations_mocks``, the output will
also show how to call the command-line interface for the various
clustering measures.

Calling from Python
-------------------

If all went well, the codes can be directly called from ``python``.
Please see ``Corrfunc/call_correlation_functions.py`` and
``Corrfunc/call_correlation_functions_mocks.py`` for examples on how to
use the Python interface. Here are a few examples:

.. code:: python

    from __future__ import print_function
    import os.path as path
    import numpy as np
    import Corrfunc
    from Corrfunc._countpairs import countpairs_wp as wp

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

    # Use a file with histogram bins, containing Nbins pairs of (rmin rmax)
    binfile = path.join(path.dirname(path.abspath(Corrfunc.__file__)), "../xi_theory/tests/", "bins")

    # Call wp
    wp_results = wp(boxsize, pimax, nthreads, binfile, x, y, z)

    # Print the results
    print("###########################################")
    print("##   rmin       rmax        wp       npairs")
    print("###########################################")
    for wp in wp_results:
        print("{0:10.4f} {1:10.4f} {2:12.6f} {3:8d}"
              .format(wp[0], wp[1], wp[3], wp[4]))
                                                        

Benchmark against Existing Codes
================================

Please see this
`gist <https://gist.github.com/manodeep/cffd9a5d77510e43ccf0>`__ for
some benchmarks with current codes.

Common Code options for both Mocks and Cosmological Boxes
=========================================================

1. ``DOUBLE_PREC`` -- does the calculations in double precision.
   Disabled by default.

2. ``USE_AVX`` -- uses the AVX instruction set found in Intel/AMD CPUs
   >= 2011 (Intel: Sandy Bridge or later; AMD: Bulldozer or later).
   Enabled by default - code will run much slower if the CPU does not
   support AVX instructions. The ``Makefile`` will automatically check
   for "AVX" support and disable this option for unsupported CPUs. 

3. ``USE_OMP`` -- uses OpenMP parallelization. Scaling is great for DD
   (perfect scaling up to 12 threads in my tests) and okay (runtime
   becomes constant ~6-8 threads in my tests) for ``DDrppi`` and ``wp``.
   Enabled by default. The ``Makefile`` will compare the `CC` variable with
   known OpenMP enabled compilers and set compile options accordingly. 

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

Author
======

Corrfunc is written/maintained by Manodeep Sinha. Please contact the
`author <mailto:manodeep@gmail.com>`__ in case of any issues.

Citing
======

If you use the code, please cite using the Zenodo DOI. The BibTex entry
for the code is

::

   @misc{manodeep_sinha_2016_55161,
       author       = {Manodeep Sinha},
       title        = {Corrfunc: Corrfunc-1.1.0},
       month        = jun,
       year         = 2016,
       doi          = {10.5281/zenodo.55161},
       url          = {http://dx.doi.org/10.5281/zenodo.55161}
   }
       
Mailing list
============

If you have questions or comments about the package, please do so on the
mailing list: https://groups.google.com/forum/#!forum/corrfunc

LICENSE
=======

Corrfunc is released under the MIT license. Basically, do what you want
with the code including using it in commercial application.

Project URL
===========

-  website (https://manodeep.github.io/Corrfunc/)
-  version control (https://github.com/manodeep/Corrfunc)

.. |Release| image:: https://img.shields.io/github/release/manodeep/Corrfunc.svg
   :target: https://github.com/manodeep/Corrfunc/releases/latest
.. |MIT licensed| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://raw.githubusercontent.com/manodeep/Corrfunc/master/LICENSE
.. |DOI| image:: https://zenodo.org/badge/19184/manodeep/Corrfunc.svg
   :target: https://zenodo.org/badge/latestdoi/19184/manodeep/Corrfunc
.. |Travis Build| image:: https://travis-ci.org/manodeep/Corrfunc.svg?branch=master
   :target: https://travis-ci.org/manodeep/Corrfunc
.. |Issues| image:: https://img.shields.io/github/issues/manodeep/Corrfunc.svg
   :target: https://github.com/manodeep/Corrfunc/issues
.. |Coverity| image:: https://img.shields.io/coverity/scan/6982.svg
   :target: https://scan.coverity.com/projects/manodeep-corrfunc
