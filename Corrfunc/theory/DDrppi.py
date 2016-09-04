#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the pair counter in
``xi_theory/xi_rp_pi``. `~Corrfunc.theory.DDrppi`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DDrppi', )


def DDrppi(autocorr, nthreads, pimax, binfile, X1, Y1, Z1,
           periodic=True, X2=None, Y2=None, Z2=None, verbose=False,
           boxsize=0.0, output_rpavg=False, c_api_timer=False, isa='fastest'):
    """
    Calculate the 3-D pair-counts corresponding to the real-space correlation
    function, :math:`\\xi(r_p, \pi)` or :math:`\\wp(r_p)`. Pairs which are
    separated by less than the ``rp`` bins (specified in ``binfile``) in the
    X-Y plane, and less than ``pimax`` in the Z-dimension are
    counted.

    Note, that this module only returns pair counts and not the actual
    correlation function :math:`\\xi(r_p, \pi)`. See the
    ``xi_theory/xi_rp_pi/wprp.c`` for computing :math:`wp(r_p)` from
    the pair counts returned.
     
    Parameters:
    -----------
    autocorr: boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.

    nthreads: integer
        The number of OpenMP threads to use. Has no effect if OpenMP was not
        enabled during library compilation.

    pimax: double
       A double-precision value for the maximum separation along
       the Z-dimension. Note that only pairs with ``0 <= dz < pimax``
       are counted (no equality).

       Distances along the :math:``\\pi`` direction are binned with unit
       depth. For instance, if ``pimax=40``, then 40 bins will be created
       along the ``pi`` direction.

    binfile: string or an list/array of floats
       For string input: filename specifying the ``rp`` bins for ``DDrppi``.
       The file should contain white-space separated values of (rpmin, rpmax)
       for each ``rp`` wanted. The bins do not need to be contiguous but must
       be in increasing order (smallest bins come first).

       For array-like input: A sequence of ``rp`` values that provides the
       bin-edges. For example, ``np.logspace(0.1, 10.0, 15)`` is a valid
       input, specifying 15 (logarithmic) bins between 0.1 and 10.0. This
       array does not need to be sorted.

    X1/Y1/Z1: array-like, real (float/double)
        The array of X/Y/Z positions for the first set of points.
        Calculations are done in the precision of the supplied arrays.

    X2/Y2/Z2: array-like, real (float/double)
    Array of XYZ positions for the second set of points. *Must* be the same
    precision as the X1/Y1/Z1 arrays. Only required when ``autocorr==0``.

    periodic: boolean
        Boolean flag to indicate periodic boundary conditions.

    verbose: boolean (default false)
       Boolean flag to control output of informational messages

    boxsize: double
        Present to facilitate exact calculations for periodic wrapping.
        If boxsize is not supplied, then the wrapping is done based on
        the maximum difference within each dimension of the X/Y/Z arrays.
    
    output_rpavg: boolean (default false)
       Boolean flag to output the average ``rp`` for each bin. Code will
       run slower if you set this flag. Also, note, if you are calculating
       in single-precision, ``rpavg`` will suffer from numerical loss of
       precision and can not be trusted. If you need accurate ``rpavg``
       values, then pass in double precision arrays for the particle positions.

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
       
    Returns:
    --------
    results: Numpy structured array

       A numpy structured array containing [rmin, rmax, ravg, pimax, npairs]
       for each radial bin specified in the ``binfile``. If ``output_ravg``
       is not set, then ``ravg`` will be set to 0.0 for all bins. ``npairs``
       contains the number of pairs in that bin and can be used to compute the
       actual :math:`\\xi(r_p, \pi)` by combining with (DR, RR) counts.

       if ``c_api_timer`` is set, then the return value is a tuple containing
       (results, api_time). ``api_time`` measures only the time spent within
       the C library and ignores all python overhead.

    Example:
    --------

    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.theory import DDrppi
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                        "../xi_theory/tests/", "bins")

    >>> N = 10000
    >>> boxsize = 420.0
    >>> nthreads = 4
    >>> autocorr = 1
    >>> pimax = 40.0
    >>> X = np.random.uniform(0, boxsize, N)
    >>> Y = np.random.uniform(0, boxsize, N)
    >>> Z = np.random.uniform(0, boxsize, N)
    >>> results = DDrppi(autocorr, nthreads, pimax, binfile,
                         X, Y, Z,
                         verbose=True, output_ravg=True)

    """
    try:
        from Corrfunc._countpairs import countpairs_rp_pi as DDrppi_extn
    except ImportError:
        msg = "Could not import the C extension for the 3-D "\
              "real-space pair counter."
        raise ImportError(msg)

    import numpy as np
    from Corrfunc.utils import translate_isa_string_to_enum,\
        return_file_with_rbins
    from future.utils import bytes_to_native_str

    if autocorr == 0:
        if X2 is None or Y2 is None or Z2 is None:
            msg = "Must pass valid arrays for X2/Y2/Z2 for "\
                "computing cross-correlation"
            raise ValueError(msg)

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    if autocorr == 1:
        extn_results, api_time = DDrppi_extn(autocorr, nthreads,
                                             pimax, rbinfile,
                                             X1, Y1, Z1,
                                             periodic=periodic,
                                             output_rpavg=output_rpavg,
                                             verbose=verbose,
                                             boxsize=boxsize,
                                             c_api_timer=c_api_timer,
                                             isa=integer_isa)
    else:
        extn_results, api_time = DDrppi_extn(autocorr, nthreads,
                                             pimax, rbinfile,
                                             X1, Y1, Z1,
                                             X2, Y2, Z2,
                                             periodic=periodic,
                                             verbose=verbose,
                                             boxsize=boxsize,
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
                              (bytes_to_native_str(b'pimax'), np.float),
                              (bytes_to_native_str(b'npairs'), np.uint64)])
    
    nbin = len(extn_results)
    results = np.zeros(nbin, dtype=results_dtype)
    for ii, r in enumerate(extn_results):
        results['rmin'][ii] = r[0]
        results['rmax'][ii] = r[1]
        results['rpavg'][ii] = r[2]
        results['pimax'][ii] = r[3]
        results['npairs'][ii] = r[4]

    if not c_api_timer:
        return results
    else:
        return results, api_time

    
if __name__ == '__main__':
    import numpy as np
    import Corrfunc
    import time
    from os.path import dirname, abspath, join as pjoin
    binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                    "../xi_theory/tests/", "bins")
    
    N = 100000
    boxsize = 420.0
    nthreads = 4
    autocorr = 1
    pimax = 40.0
    seed = 42
    np.random.seed(seed)
    X = np.random.uniform(0, boxsize, N)
    Y = np.random.uniform(0, boxsize, N)
    Z = np.random.uniform(0, boxsize, N)
    t0 = time.time()
    results, api_time = DDrppi(autocorr, nthreads, pimax, binfile,
                               X, Y, Z,
                               verbose=True,
                               c_api_timer=True,
                               output_rpavg=True)

    t1 = time.time()
    print("Results from DDrppi (Npts = {0}): Time taken = {1:0.3f} sec "
          "Python time = {2:0.3f} sec".format(N, api_time, t1-t0))

    # Only print the first ten bins (1st rp bin, with 10 pimax bins)
    for r in results[0:10]:
        print("{0}".format(r))

