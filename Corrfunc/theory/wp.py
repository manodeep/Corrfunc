#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the theoretical projected
auto-correlation function, wp(rp), in ``theory/wp``. This python
wrapper is in :py:mod:`Corrfunc.theory.wp`.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('wp',)


def wp(boxsize, pimax, nthreads, binfile, X, Y, Z, verbose=False,
       output_rpavg=False, c_api_timer=False, isa='fastest'):
    """
    Function to compute the projected correlation function in a
    periodic cosmological box. Pairs which are separated by less
    than the ``rp`` bins (specified in ``binfile``) in the
    X-Y plane, and less than ``pimax`` in the Z-dimension are
    counted.

    Note that pairs are double-counted. And if ``rpmin`` is set to
    0.0, then all the self-pairs (i'th particle with itself) are
    added to the first bin => minimum number of pairs in the first bin
    is the total number of particles.

    Parameters
    -----------

    boxsize: double
       A double-precision value for the boxsize of the simulation
       in same units as the particle positions and the ``rp`` bins.

    pimax: double
       A double-precision value for the maximum separation along
       the Z-dimension. Note that only pairs with ``0 <= dz < pimax``
       are counted (no equality).

    nthreads: integer
       Number of threads to use.

    binfile: string or an list/array of floats
       For string input: filename specifying the ``rp`` bins for
       ``DDrppi_mocks``. The file should contain white-space separated values
       of (rpmin, rpmax)  for each ``rp`` wanted. The bins do not need to be
       contiguous but must be in increasing order (smallest bins come first).

       For array-like input: A sequence of ``rp`` values that provides the
       bin-edges. For example,
       ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
       input, specifying 15 (logarithmic) bins between 0.1 and 10.0. This
       array does not need to be sorted.

    X/Y/Z: arraytype, real (float/double)
       Particle positions in the 3 axes. Must be within [0, boxsize]
       and specified in the same units as ``rp_bins`` and boxsize. All
       3 arrays must be of the same floating-point type.

       Calculations will be done in the same precision as these arrays,
       i.e., calculations will be in floating point if XYZ are single
       precision arrays (C float type); or in double-precision if XYZ
       are double precision arrays (C double type).

    verbose: boolean (default false)
       Boolean flag to control output of informational messages

    output_rpavg: boolean (default false)
       Boolean flag to output the average ``rp`` for each bin. Code will
       run slower if you set this flag. Also, note, if you are calculating
       in single-precision, ``rpavg`` will suffer from numerical loss of
       precision and can not be trusted. If you need accurate ``rpavg``
       values, then pass in double precision arrays for ``XYZ``.

    c_api_timer: boolean (default false)
       Boolean flag to measure actual time spent in the C libraries. Here
       to allow for benchmarking and scaling studies.

    isa: string (default ``fastest``)
       Controls the runtime dispatch for the instruction set to use. Possible
       options are: [``fastest``, ``avx``, ``sse42``, ``fallback``]

       Setting isa to ``fastest`` will pick the fastest available instruction
       set on the current computer. However, if you set ``isa`` to, say,
       ``avx`` and ``avx`` is not available on the computer, then the code will
       revert to using ``fallback`` (even though ``sse42`` might be available).

       Unless you are benchmarking the different instruction sets, you should
       always leave ``isa`` to the default value. And if you *are*
       benchmarking, then the string supplied here gets translated into an
       ``enum`` for the instruction set defined in ``utils/defs.h``.

    Returns
    --------

    results: Numpy structured array

       A numpy structured array containing [rpmin, rpmax, rpavg, wp, npairs]
       for each radial specified in the ``binfile``. If ``output_rpavg`` is not
       set then ``rpavg`` will be set to 0.0 for all bins. ``wp`` contains the
       projected correlation function while ``npairs`` contains the number of
       unique pairs in that bin.

       if ``c_api_timer`` is set, then the return value is a tuple containing
       (results, api_time). ``api_time`` measures only the time spent within
       the C library and ignores all python overhead.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.theory.wp import wp
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../theory/tests/", "bins")
    >>> N = 10000
    >>> boxsize = 420.0
    >>> pimax = 40.0
    >>> nthreads = 4
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> X = np.random.uniform(0, boxsize, N)
    >>> Y = np.random.uniform(0, boxsize, N)
    >>> Z = np.random.uniform(0, boxsize, N)
    >>> results = wp(boxsize, pimax, nthreads, binfile, X, Y, Z)
    >>> for r in results: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f} "
    ...                         "{4:10d}".format(r['rmin'], r['rmax'],
    ...                         r['rpavg'], r['wp'], r['npairs']))
    ...                   # doctest: +NORMALIZE_WHITESPACE
    0.167536   0.238755   0.000000  66.702471         18
    0.238755   0.340251   0.000000 -15.792466         16
    0.340251   0.484892   0.000000   2.990170         42
    0.484892   0.691021   0.000000 -15.786307         66
    0.691021   0.984777   0.000000 -11.973531        142
    0.984777   1.403410   0.000000  -9.706936        298
    1.403410   2.000000   0.000000 -11.705601        588
    2.000000   2.850200   0.000000   3.839991       1466
    2.850200   4.061840   0.000000  -0.929359       2808
    4.061840   5.788530   0.000000   0.446806       5802
    5.788530   8.249250   0.000000   1.420201      11926
    8.249250  11.756000   0.000000  -1.075778      23478
    11.756000  16.753600   0.000000  -0.561263      47994
    16.753600  23.875500   0.000000  -0.094424      98042

    """

    try:
        from Corrfunc._countpairs import countpairs_wp as wp_extn
    except ImportError:
        msg = "Could not import the C extension for the projected "\
              "correlation function."
        raise ImportError(msg)

    import numpy as np
    from future.utils import bytes_to_native_str
    from Corrfunc.utils import translate_isa_string_to_enum,\
        return_file_with_rbins

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    extn_results, api_time = wp_extn(boxsize, pimax, nthreads, rbinfile,
                                     X, Y, Z,
                                     verbose=verbose,
                                     output_rpavg=output_rpavg,
                                     c_api_timer=c_api_timer,
                                     isa=integer_isa)
    if extn_results is None:
        msg = "RuntimeError occurred"
        raise RuntimeError(msg)

    if delete_after_use:
        import os
        os.remove(rbinfile)

    results_dtype = np.dtype([(bytes_to_native_str(b'rmin'), np.float),
                              (bytes_to_native_str(b'rmax'), np.float),
                              (bytes_to_native_str(b'rpavg'), np.float),
                              (bytes_to_native_str(b'wp'), np.float),
                              (bytes_to_native_str(b'npairs'), np.uint64)])

    nbin = len(extn_results)
    results = np.zeros(nbin, dtype=results_dtype)

    for ii, r in enumerate(extn_results):
        results['rmin'][ii] = r[0]
        results['rmax'][ii] = r[1]
        results['rpavg'][ii] = r[2]
        results['wp'][ii] = r[3]
        results['npairs'][ii] = r[4]

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import doctest
    doctest.testmod()
