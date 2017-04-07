.. _quickstart:

******************************
Getting started with Corrfunc
******************************

Corrfunc is a set of high-performance routines to measure clustering
statistics. The codes are divided conceptually into two different segments:

* theory - calculates clustering statistics on **simulation** volumes. Input
  positions are expected to be Cartesian X/Y/Z. Periodic boundary conditions
  are supported. Relevant C codes are in directory ``theory/``
  
* mocks - calculates clustering statistics on **observation** volumes. Input
  positions are assumed to be in obverser frame, ``Right Ascension``, ``Declination``
  and ``SpeedofLight*Redshift`` (where required; :math:`\omega(\theta)`
  only needs ``RA`` and ``DEC``). Relevant C codes are in directory ``mocks/``

This getting-started guide assumes you have already followed the
:ref:`step_by_step_install` section of the documentation to get the package 
and its dependencies set up on your machine. 

If you want to compute correlation functions and have installed the python
extensions, then see :ref:`function_usage_examples` for typical
tasks. Otherwise, read on for the various interfaces available within Corrfunc.


Computing Clustering Statistics with Corrfunc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
Corrfunc supports three separate mechanisms to compute the clustering statistics:

* **Via python** (if you have ``python`` and ``numpy`` installed)

  Pros: Fully flexible API to modulate code behaviour at runtime. For instance,
  calculations can be performed in double-precision simply by passing arrays of
  doubles (rather than floats).
  
  Cons: Has fixed python overhead. For low particle numbers, can be as much as
  20% slower compared to the command-line  executables.

  See :ref:`python-interface` for details on how to use the python interface.
    
* **Via static libraries directly in C codes**

  Pros: Fully flexible API to modulate code behaviour at runtime. All features supported by the python extensions are also supported here. 

  Cons: Requires coding in C. See example C codes invoking the ``theory`` and
  ``mocks`` in the directories: ``theory/examples/run_correlations.c`` and ``mocks/examples/run_correlations_mocks.c``.

  See :ref:`staticlibrary-interface` for details on how to use the static library interface.
  
* **Command-line executables**
  
  Pros: Fastest possible implementations of all clustering statistics

  Cons: API is fixed. Any changes require full re-compilation. 

  See :ref:`commandline-interface` for details on how to use the command-line executables.

Available Corrfunc interfaces 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
.. toctree::
   :maxdepth: 1

   python-interface
   staticlibrary-interface
   commandline-interface
   all-interfaces








