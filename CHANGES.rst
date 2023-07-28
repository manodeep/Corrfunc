3.0.0 (future)
=================

**Breaking Changes**
---------------------
- Package will be renamed to ``corrfunc`` from ``Corrfunc``

New features
------------
- conda installable package
- GPU version


2.5.1 (28/07/2023)
==================

Enhancements
------------
- Corrfunc now compiles  and runs on Apple M1/M2 cpus (using the `FALLBACK` kernels) [#295]

Changes
-------
- Python >= 3.7 and numpy >= 1.16 are required for python extensions [#291]

Enhancements
------------
- Warn about loss of precision for float32 calculations involving small ``theta`` in ``DDtheta_mocks`` and large ``mu`` in ``DDsmu_mocks`` [#299]

2.5.0 (2022-12-23)
================

Enhancements
------------
- Allow user to specify periodicity and box size per dimension [#276]
- Allow larger Rmax (up to half the boxsize) [#277]

Changes
-------
- Add Corrfunc/tests.py to CI testing [#260]
- Migrate doctests to Python 3.8 [#261]
- Migrate Python tests to pytest [#265]

Fixes
-----
- Add additional check to tell if it's safe to redirect stdout/err [#270]
- Check and fix ``z`` vs ``cz`` in ``DDrppi_mocks`` and ``DDsmu_mocks`` only if comoving distance flag is not set [#275]
- Update GNU assembler bug detection [#278]
- Fix installation instructions and update README.rst [#285]
- Remove mistaken references to "projected" correlation function in docs [#289]


2.4.0 (2021-09-30)
==================
This release adds the ``boxsize`` parameter to the command line interfaces and
requires the user to specify the box size in the Python interfaces to the periodic
theory functions.  It also contains a number of performance, code-quality, and
user-experience improvements.

**Breaking Changes**
--------------------
- Require user to specify ``boxsize`` rather than automatically detect particle
  extent in periodic theory boxes. Applies to Python, command line, and C API.
  The order of some Python keyword args has also changed. [#199]

Enhancements
------------
- In the theoretical VPF calculation (``theory.vpf``), the total volume of the random spheres can now exceed the volume of the sample  [#238]
- Gridlink (the binning of particles into cells) now uses a parallel algorithm for the theory module [#239]
- Add detection of known-bad Cray hugepages library at NERSC [#246]
- Replace ``np.float`` with ``np.float64`` to fix numpy 1.20 deprecation [#250]
- Test Numpy versions as old as 1.14 and recent as 1.20 [#251]
- Add lscpu and preprocessor defs to CI output [#259]

Bug fixes
---------
- Fix Python reference leak to results struct [#229]
- Fix parsing error when ``periodic=False`` and ``boxsize`` not given in the theory module [#257]


2.3.4 (2019-07-21)
==================
This is a bug-fix release and contains general code quality improvements.


Enhancements
------------
- A new helper routine to find the combination of (RA, DEC) refinements that produces fastest runtime in ``DDtheta_mocks`` [#216]
- Further testing via GitHub Actions [#220]
- Added Ubuntu-Xenial on Travis [#222]

Bug fixes
----------
- Fixing docs build failure on Travis [#215]
- Fixing compile failure on missing 'CC' in environment [#226]


2.3.3 (2019-02-03)
==================
This is a bug-fix release and contains general code quality improvements.


Enhancements
------------



Bug fixes
----------
- Installation does not require python(3)-config anymore [#209, #211]
- Better handling of terminal colours for unknown terminals [#209]
- Prevent incorrect calculations with periodic boundaries for large ratios of (zmax, Rmax) to Lbox [#210]


2.3.2 (2019-12-24)
===================
This is a release for bug-fixes and general code quality improvements. Travis
now also tests for ``python3.7``.

Enhancements
------------
- Improved code quality and adherence to PEP8 [#189]
- Documentation no longer shows duplicate entries [#205]


Bug fixes
----------
- Incorrect calculations for non-native endian data [#191]
- Large Rmax to Lbox ratio now supported for periodic boundaries [#192]
- Workaround for GNU Assembler bug causing incorrect calculations [#196]
- Only report runtime isa support if we also have compiler support [#200]
- Example code to illustrate how to code custom weights with AVX512F [#205]

2.3.1 (2019-06-21)
================

Enhancements
------------
- Reduce memory footprint of the cell pairs [#186]


2.3.0 (2019-05-20)
==================

**Breaking Changes**
--------------------

New features
------------
- AVX512F kernels for all pair-counters, faster code from new optimizations using the minimum separation between pairs of cells, option to use the input particle arrays directly and not make a copy of the particle positions, internal code changes to (hopefully) achieve better OpenMP scaling [#167, #170, #173]

Bug fixes
---------
- Fix segmentation fault in vpf_mocks [#168]
- Fix automatic uniform weights array when only one set of weights (or a scalar) is passed [#180]
- Fix memory leak due to Python reference leak when using weights [#181]


2.2.0 (2018-08-18)
==================

**Breaking Changes**
--------------------
- Drop Python 2.6 support

New features
------------
- Progress bar is displayed in Jupyter notebooks [#158]

Bug fixes
---------
- Fix virtualenv install issue [#159]
- Error messages are displayed in Jupyter notebooks
  instead of the unhelpful "TypeError: 'NoneType' object is not iterable". [#158]


2.1.0 (2018-08-17)
==================

New features
------------
- New pair counter `DD(s, mu)` for theory and mocks (contributed by @nickhand,
  in #130 and #132) [#166]


Enhancements
------------
- GSL version now specified and tested by Travis [#164]
- Now possible to specify the number of Newton-Raphson steps to
  improve accuracy of approximate reciprocals. Available in `DD(rp, pi)` for mocks,
  and `DD(s, mu)` for both theory and mocks


2.0.0 (2017-04-06)
==================

New features
------------

- Library behaviour can now be controlled at runtime
- Calculates with ``doubles`` and ``floats`` transparently
  (passing arrays of ``doubles`` ensures calculation in double
  precision)
- Both the API and ABI should be future proof
- Extensive docs (first version with docs)
- Arbitrary cosmology can be accounted for in the mocks routines  `#71 <https://github.com/manodeep/Corrfunc/issues/71>`_

**Breaking Changes**
---------------------

- API has changed from previous version. Two additional inputs are
  now required for every statistic (`#73 <https://github.com/manodeep/Corrfunc/issues/73>`_)


Enhancements
------------

- Ctrl-C now aborts even within python extensions (cleans up memory too!, `see issue #12 <https://github.com/manodeep/Corrfunc/issues/12>`_)
- Significantly improved installation for python

  - compiler can now be specified within ``python setup.py install CC=yourcompiler``
    `(see issue #31) <https://github.com/manodeep/Corrfunc/issues/31>`_
  - python via an alias is now solved `#52 <https://github.com/manodeep/Corrfunc/issues/52>`_


Bug fixes
----------

- Fixed bug in ``DDrppi_mocks`` where the minimum number of grid cells had to
  be 1 `#70 <https://github.com/manodeep/Corrfunc/issues/70>`_



Outstanding issues
-------------------
- Conda package still is pending (`#49 <https://github.com/manodeep/Corrfunc/issues/49>`_)
- Recursive Makefile needs to be replaced with
  a more monolithic Makefile (`#14 <https://github.com/manodeep/Corrfunc/issues/14>`_)
- Parameter parsing in python extensions can be flaky (`#79 <https://github.com/manodeep/Corrfunc/issues/79>`_)


1.1.0 (2016-06-08)
===================

- SSE kernels for all statistics
- Incorrect normalization in ``xi``. **ALL** previous
  ``xi`` calculations were wrong.


1.0.0 (2016-04-14)
==================

- Improved installation process
- Detecting ``AVX`` capable CPU at compile time
- Double-counting bug fixes in ``wp`` and ``xi``


0.2.3 (2016-03-30)
==================

- Streamlined compilation on MACs
- PyPI version is not verbose by default


0.2.2 (2016-02-09)
==================

- First version on `PyPI <https://pypi.python.org/pypi/Corrfunc>`_


0.2.1 (2016-02-06)
==================

- ``AVX`` enabled by default


0.2.0 (2016-02-05)
==================

- Python 2/3 compatible



0.0.1 (2015-11-11)
==================

- Initial release
