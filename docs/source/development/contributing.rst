.. _contributing:

=========================
Contributing to Corrfunc
=========================
Corrfunc is written in a very modular fashion with minimal interaction between
the various calculations. The algorithm presented in Corrfunc is applicable to
a broad-range of astrophysical problems, viz., any situation that requires
looking at *all* objects around a target and performing some analysis with
this group of objects.

Here are the basic steps to get your statistic into the Corrfunc package:

* Fork the repo and add your statistic
* Add exhaustive tests. The output of your statistic should **exactly** agree with a
  brute-force implementation (under double-precision). Look at ``test_periodic.c`` and ``test_nonperiodic.c``
  under ``theory/tests/`` for tests on simulation volumes. For mock
  catalogs, look at ``mocks/tests/tests_mocks.c``.
* Add a python extension for the new statistic. This extension should reside in file
  ``theory/python_bindings/_countpairs.c`` or
  ``mocks/python_bindings/_countpairs_mocks.c`` for statistics relevant for
  simulations and mocks respectively. It is preferred to have the extension
  documented but not necessary.
* Add a call to this new *extension* in the
  ``python_bindings/call_correlation_functions*.py`` script.

.. note:: Different from corresponding script in ``Corrfunc/`` directory.

* Add a python wrapper for the previous python extension. This wrapper should
  exist in ``Corrfunc/theory/`` or ``Corrfunc/mocks/``. Wrapper **must** have
  inline API docs.
* Add the new wrapper to ``__all__`` in ``__init__.py`` within the relevant
  directory.
* Add an example call to this *wrapper* in
  ``Corrfunc/call_correlation_functions.py`` or
  ``Corrfunc/call_correlation_functions_mocks.py`` for simulations and mocks
  respectively.
  
.. note:: Different from corresponding script in ``python_bindings`` directory.
          
* Add the new wrapper to the API docs within
  ``ROOT_DIR/docs/source/theory_functions.rst`` or
  ``ROOT_DIR/docs/source/mocks_functions.rst``. 
* Add to the contributors list under
  ``ROOT_DIR/docs/source/development/contributors.rst``.
*  Submit pull request
   

.. note:: Please feel free to email the `author <mailto:manodeep@gmail.com>`_ or
          the `Corrfunc Google Groups <https://groups.google.com/forum/#!forum/corrfunc>`_ if you need help at any stage. 


Corrfunc Design
~~~~~~~~~~~~~~~~
All of the algorithms in Corrfunc have the following components:

* Reading in data. Relevant routines are in the ``io/`` directory with a
  mapping within ``io.c`` to handle the file format
* Creating the 3-D lattice structure. Relevant routines are in the
  ``utils/gridlink_impl.c.src``  and ``utils/gridlink_mocks.c.src``. This
  lattice grids up the particle distribution on cell-sizes of ``rmax`` (the
  maximum search radius).

  
.. note:: The current lattice code duplicates the particle memory. If you
  need a lattice that does not duplicate the particle memory, then please email
  the `author <mailto:manodeep@gmail.com>`_. Relevant code existed in Corrfunc
  but has been removed in the current incarnation.

  
* Setting up the OpenMP sections such that threads have local copies of
  histogram arrays. If OpenMP is not enabled, then this section should not
  produce any compilable code.
* Looping over all cells in the 3-D lattice and then looping over all
  neighbouring cells for each cell.
* For a pair of cells, hand over the two sets of arrays into a specialized
  kernel (``count*kernel.c.src``) for computing pairs.  
* Aggregate the results, if OpenMP was enabled.


Directory and file layout
~~~~~~~~~~~~~~~~~~~~~~~~~~

* Codes that compute statistics on simulation volumes (Cartesian XYZ as input)
  go into a separate directory within ``theory``
* Codes that compute statistics on mock catalogs (RA, DEC [CZ]) go into a
  separate directory within ``mocks``
* Public API in a ``count*.h`` file. Corresponding C file simply dispatches to
  appropriate floating point implementation.
* Floating point implmentation in file ``count*_impl.c.src``. This file is
  processed via ``sed`` to generate both single and double precision
  implementations.
* A kernel named ``count*kernels.c.src`` containing implementations for
  counting pairs on two sets of arrays. This kernel file is also preprocessed
  to produce both the single and double precision kernels.
* Tests go within ``tests`` directory under ``theory`` or ``mocks``, as
  appropriate. For simulation routines, tests with and without periodic
  boundaries go into ``test_periodic.c`` and ``test_nonperiodic.c``
* C code to generate the python extensions goes under ``python_bindings``
  directory into the file ``_countpairs*.c``
* Each python extension has a python wrapper within ``Corrfunc`` directory

Coding Guidelines
~~~~~~~~~~~~~~~~~

C guidelines
============

Code contents
-------------

* **Always** check for error conditions when calling a function
* If an error condition occurs when making an kernel/external library call,
  first call ``perror`` and then return the error status. If calling a wrapper
  from within Corrfunc, assume that ``perror`` has already been called and
  simply return the status. Clean up memory before returning status.
* Declare variables in the smallest possible scope.
* Add ``const`` qualifiers liberally
* There **must** not be any compiler warnings (with ``gcc6.0``) under the given set of Warnings
  already enabled within ``common.mk``. If the warning can not be avoided
  because of logic issues, then suppress the warning but note why that
  suppression is required. Warnings are treated as errors on the continuous integration platform (TRAVIS)
* Valgrind should not report any fixable memory or file leaks (memory
  leaks in OpenMP library, e.g., ``libgomp``, are fine)

Style
------
The coding style is loosely based on `Linux Kernel Guideline
<https://www.kernel.org/doc/Documentation/CodingStyle>`_. These are recommended
but not strictly enforced. However, note that if you do contribute code to
Corrfunc, the style may get converted. 

* Braces
  - Opening braces start at the same line, except for functions
  - Closing braces on new line
  - Even single line conditionals must have opening and closing braces
    
* Comments
  - Explanatory comments on top of code segment enclosed with ``/**/``
  - Inline comments must be single-line on the right 

* Indentation is ``tab:=4 spaces``

* Avoid ``typedef`` for ``structs`` and ``unions``

Python guidelines
=================

* Follow the `astropy python code guide <http://docs.astropy.org/en/stable/development/codeguide_emacs.html>`_
* Docs are in ``numpydocs`` format. Follow any of the wrapper routines in
  ``Corrfunc`` (which are, in turn, taken from `halotools <http://halotools.readthedocs.io/>`_)

