.. _step_by_step_install:

************************
Package Installation
************************

To install Corrfunc, you can either use pip or clone the repo from GitHub and build the source code.
Either way, be sure to read the :ref:`Corrfunc_dependencies` section prior to installation.

Using pip
====================

The simplest way to install the latest release of the code is with pip. Before installation, be sure you have installed the package dependencies described in the :ref:`Corrfunc_dependencies` section

.. code:: python

    python -m pip install Corrfunc

This will install the latest official release of the code.
If you want the latest master branch,
you will need to build the code from source following the instructions in the next section.

Building from source
====================

If you don't install the latest release using pip,
you can instead clone the cource code and call the setup file.
Before installation, be sure you have installed the package dependencies
described in the :ref:`corrfunc_dependencies` section.
The first step is to clone the Corrfunc repository

.. code::

      git clone https://github.com/manodeep/Corrfunc.git
      cd Corrfunc
      make install
      python -m pip install . (--user)


.. _corrfunc_dependencies:

Dependencies
============

The command-line version of Corrfunc needs the following packages to be installed:

- `make <https://www.gnu.org/software/make/>`_: 3.80 or later
- `C compiler <https://gcc.gnu.org/>`_: gcc >=4.6, clang, icc. Multi-threading
  will be disabled if the compiler does not support OpenMP.
- `gsl <https://www.gnu.org/software/gsl/>`_: any recent version


If you plan to use the C extensions, then the following are required:

- `Python <http://www.python.org/>`_: 2.7 or later
- `Numpy <http://www.numpy.org/>`_: 1.7 or later

Any of the above can be installed with either pip or conda.

.. _verifying_your_installation:

Verifying your installation
==============================

After installing Corrfunc, you should run the integrated test suite to make
sure that the package was installed correctly. If you installed from source,
then run the following sequence from the repository root directory:

.. code:: console

    make tests  # run the C tests
    python -m pip install pytest
    python -m pytest  # run the Python tests

If you installed using pip/conda, then use the following to run the tests
from your shell:

.. code:: console

    python -m pytest --pyargs Corrfunc
          

Once you have installed the package, see :ref:`quickstart` for instructions on how to get up and running.
