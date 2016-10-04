#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the theoretical 3-D
real-space correlation function, :math:`\\xi(r)`. Corresponding
C routines are in ``theory/xi/``, python interface is
:py:mod:`Corrfunc.theory.xi`.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('xi',)


def xi(boxsize, nthreads, binfile, X, Y, Z, verbose=False,
       output_ravg=False, xbin_refine_factor=2, ybin_refine_factor=2,
       zbin_refine_factor=1, c_api_timer=False, isa='fastest'):
    """
    Function to compute the projected correlation function in a
    periodic cosmological box. Pairs which are separated by less
    than the ``r`` bins (specified in ``binfile``) in 3-D real space.

    Note that pairs are double-counted. And if ``rmin`` is set to
    0.0, then all the self-pairs (i'th particle with itself) are
    added to the first bin => minimum number of pairs in the first bin
    is the total number of particles.

    Parameters
    -----------

    boxsize: double
       A double-precision value for the boxsize of the simulation
       in same units as the particle positions and the ``r`` bins.

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

    output_ravg: boolean (default false)
       Boolean flag to output the average ``r`` for each bin. Code will
       run slower if you set this flag. Also, note, if you are calculating
       in single-precision, ``ravg`` will suffer from numerical loss of
       precision and can not be trusted. If you need accurate ``ravg``
       values, then pass in double precision arrays for ``XYZ``.

    (xyz)bin_refine_factor: integer, default is (2,2,1); typically within [1-3]
       Controls the refinement on the cell sizes. Can have up to a 20% impact
       on runtime.

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

       A numpy structured array containing [rmin, rmax, ravg, xi, npairs] for
       each radial specified in the ``binfile``. If ``output_ravg`` is not
       set then ``ravg`` will be set to 0.0 for all bins. ``xi`` contains the
       correlation function while ``npairs`` contains the number of
       pairs in that bin.

       if ``c_api_timer`` is set, then the return value is a tuple containing
       (results, api_time). ``api_time`` measures only the time spent within
       the C library and ignores all python overhead.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.theory.xi import xi
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../theory/tests/", "bins")
    >>> N = 100000
    >>> boxsize = 420.0
    >>> nthreads = 4
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> X = np.random.uniform(0, boxsize, N)
    >>> Y = np.random.uniform(0, boxsize, N)
    >>> Z = np.random.uniform(0, boxsize, N)
    >>> results = xi(boxsize, nthreads, binfile, X, Y, Z, output_ravg=True)
    >>> for r in results: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f}"
    ...                         " {4:10d}".format(r['rmin'], r['rmax'],
    ...                         r['ravg'], r['xi'], r['npairs']))
    ...                   # doctest: +NORMALIZE_WHITESPACE
    0.167536   0.238755   0.226592  -0.205741          4
    0.238755   0.340251   0.289277  -0.176737         12
    0.340251   0.484892   0.426819  -0.051838         40
    0.484892   0.691021   0.596187  -0.131862        106
    0.691021   0.984777   0.850100  -0.049217        336
    0.984777   1.403410   1.225112   0.028532       1052
    1.403410   2.000000   1.737153   0.011392       2994
    2.000000   2.850200   2.474588   0.005395       8614
    2.850200   4.061840   3.532018  -0.014108      24448
    4.061840   5.788530   5.022241  -0.010794      70996
    5.788530   8.249250   7.160648  -0.001598     207392
    8.249250  11.756000  10.207213  -0.000333     601002
    11.756000  16.753600  14.541171  -0.000003    1740084
    16.753600  23.875500  20.728773  -0.001605    5028058

    """

    try:
        from Corrfunc._countpairs import countpairs_xi as xi_extn
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
    extn_results, api_time = xi_extn(boxsize, nthreads, rbinfile,
                                     X, Y, Z,
                                     verbose=verbose,
                                     output_ravg=output_ravg,
                                     xbin_refine_factor=xbin_refine_factor,
                                     ybin_refine_factor=ybin_refine_factor,
                                     zbin_refine_factor=zbin_refine_factor,
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
                              (bytes_to_native_str(b'ravg'), np.float),
                              (bytes_to_native_str(b'xi'), np.float),
                              (bytes_to_native_str(b'npairs'), np.uint64)])

    nbin = len(extn_results)
    results = np.zeros(nbin, dtype=results_dtype)

    for ii, r in enumerate(extn_results):
        results['rmin'][ii] = r[0]
        results['rmax'][ii] = r[1]
        results['ravg'][ii] = r[2]
        results['xi'][ii] = r[3]
        results['npairs'][ii] = r[4]

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import doctest
    doctest.testmod()
