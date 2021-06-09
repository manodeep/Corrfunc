.. _read_catalog:

******************************
Reading Catalogs for Corrfunc
******************************

All of the ``Corrfunc`` routines require some sort of
position arrays, X/Y/Z, as input. These arrays are
expected to be 1-D arrays of type ``np.array``. If
you already have have the required ``numpy`` arrays,
then you can just pass them straight to ``Corrfunc``.
If you need to read the arrays in from disk, then read
on. For the command-line interface, the input files can only
be in ASCII or fast-food format (for description of fast-food
binaries, see :ref:`fast_food_binary`).

.. toctree::
   :maxdepth: 1

   fast_food_binary


Reading from ASCII files
========================

This is the most straight forward way -- you need an ASCII
file with columns X/Y/Z (white-space separated).

Using ``numpy.genfromtxt``
---------------------------

.. code:: python

          import numpy as np
          fname = "myfile_containing_xyz_columns.dat"

          # For double precision calculations
          dtype = np.float64  ## change to np.float32 for single precision

          X, Y, Z = np.genfromtxt(fname, dtype=dtype, unpack=True)


.. note:: :py:mod:`Corrfunc.read_catalog` uses this exact code-snippet to read in ASCII files in python.


Reading from fast-food files
=============================

If you are using the command-line interface, then the code will **have** to
read the arrays from files. While ``Corrfunc`` natively supports both
ASCII and fast-food formats (for description of fast-food binaries, see
:ref:`fast_food_binary`), the following python utility is intended to
read both these types of files.


Using utility: :py:mod:`Corrfunc.io.read_catalog`
-------------------------------------------------

:py:mod:`Corrfunc.io.read_catalog` can directly read ASCII files or fast-food binary
files.

.. code:: python

          from Corrfunc.io import read_catalog

          # Read the standard theory catalog (on a box)
          # supplied with Corrfunc
          X, Y, Z = read_catalog()

          # Read some other format -> have to specify
          # filename
          fname = "myfile_containing_xyz_columns.dat"
          X, Y, Z = read_catalog(fname)
