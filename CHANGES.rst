3.0.0 (upcoming)
=================

**Breaking Changes**
---------------------
- Package will be renamed to ``corrfunc`` from ``Corrfunc``

New features
------------
- New pair counter `DD(s, mu)` for theory and mocks
- conda installable package


2.1.0
=======

Enhancements
------------
- GSL version now specified and tested by Travis [#164]


2.0.0
=======

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

- Ctrl-C now aborts even within python extensions (cleans up memory too!, `#12 <https://github.com/manodeep/Corrfunc/issues/12>`_)
- Significantly improved installation for python
  - compiler can now be specified within ``python setup.py install
    CC=yourcompiler`` `#31<https://github.com/manodeep/Corrfunc/issues/31>`_
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


1.1.0 (June 8, 2016)
=====================

- SSE kernels for all statistics
- Incorrect normalization in ``xi``. **ALL** previous
  ``xi`` calculations were wrong.


1.0.0 (Apr 14, 2016)
====================

- Improved installation process  
- Detecting ``AVX`` capable CPU at compile time
- Double-counting bug fixes in ``wp`` and ``xi``
  

0.2.3 (Mar 30, 2016)
=====================

- Streamlined compilation on MACs
- PyPI version is not verbose by default


0.2.2 (Feb 9, 2016)
====================

- First version on `PyPI <https://pypi.python.org/pypi/Corrfunc>`_


0.2.1 (Feb 6, 2016)
====================

- ``AVX`` enabled by default


0.2.0 (Feb 5, 2016)
====================

- Python 2/3 compatible
 


0.0.1 (Nov 11, 2015)
====================

- Initial release

