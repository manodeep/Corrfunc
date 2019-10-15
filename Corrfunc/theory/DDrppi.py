#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the pair counter in
``theory/DDrppi/``. This wrapper is in :py:mod:`Corrfunc.theory.DDrppi`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DDrppi', )


def DDrppi(autocorr, nthreads, pimax, binfile, X1, Y1, Z1, weights1=None,
           periodic=True, X2=None, Y2=None, Z2=None, weights2=None,
           verbose=False, boxsize=0.0, output_rpavg=False,
           xbin_refine_factor=2, ybin_refine_factor=2,
           zbin_refine_factor=1, max_cells_per_dim=100,
           copy_particles=True, enable_min_sep_opt=True,
           c_api_timer=False, isa=r'fastest', weight_type=None):
    """
    Calculate the 3-D pair-counts corresponding to the real-space correlation
    function, :math:`\\xi(r_p, \pi)` or :math:`\\wp(r_p)`. Pairs which are
    separated by less than the ``rp`` bins (specified in ``binfile``) in the
    X-Y plane, and less than ``pimax`` in the Z-dimension are
    counted.

    If ``weights`` are provided, the mean pair weight is stored in the
    ``"weightavg"`` field of the results array.  The raw pair counts in the
    ``"npairs"`` field are not weighted.  The weighting scheme depends on
    ``weight_type``.

    .. note:: that this module only returns pair counts and not the actual
        correlation function :math:`\\xi(r_p, \pi)` or :math:`wp(r_p)`. See the
        utilities :py:mod:`Corrfunc.utils.convert_3d_counts_to_cf` and
        :py:mod:`Corrfunc.utils.convert_rp_pi_counts_to_wp` for computing
        :math:`\\xi(r_p, \pi)` and :math:`wp(r_p)` respectively from the
        pair counts.


    Parameters
    -----------

    autocorr: boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.

    nthreads: integer
        The number of OpenMP threads to use. Has no effect if OpenMP was not
        enabled during library compilation.

    pimax: double
        A double-precision value for the maximum separation along
        the Z-dimension.

        Distances along the :math:``\\pi`` direction are binned with unit
        depth. For instance, if ``pimax=40``, then 40 bins will be created
        along the ``pi`` direction.

        Note: Only pairs with ``0 <= dz < pimax`` are counted (no equality).

    binfile: string or an list/array of floats
        For string input: filename specifying the ``rp`` bins for
        ``DDrppi``. The file should contain white-space separated values
        of (rpmin, rpmax)  for each ``rp`` wanted. The bins need to be
        contiguous and sorted in increasing order (smallest bins come first).

        For array-like input: A sequence of ``rp`` values that provides the
        bin-edges. For example,
        ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
        input specifying **14** (logarithmic) bins between 0.1 and 10.0. This
        array does not need to be sorted.

    X1/Y1/Z1: array-like, real (float/double)
        The array of X/Y/Z positions for the first set of points.
        Calculations are done in the precision of the supplied arrays.

    weights1: array_like, real (float/double), optional
        A scalar, or an array of weights of shape (n_weights, n_positions) or
        (n_positions,). ``weight_type`` specifies how these weights are used;
        results are returned in the ``weightavg`` field.  If only one of
        weights1 and weights2 is specified, the other will be set to uniform
        weights.

    X2/Y2/Z2: array-like, real (float/double)
        Array of XYZ positions for the second set of points. *Must* be the same
        precision as the X1/Y1/Z1 arrays. Only required when ``autocorr==0``.

    weights2: array-like, real (float/double), optional
        Same as weights1, but for the second set of positions

    periodic: boolean
        Boolean flag to indicate periodic boundary conditions.

    verbose: boolean (default false)
        Boolean flag to control output of informational messages

    boxsize: double
        The side-length of the cube in the cosmological simulation.
        Present to facilitate exact calculations for periodic wrapping.
        If boxsize is not supplied, then the wrapping is done based on
        the maximum difference within each dimension of the X/Y/Z arrays.

    output_rpavg: boolean (default false)
        Boolean flag to output the average ``rp`` for each bin. Code will
        run slower if you set this flag.

        Note: If you are calculating in single-precision, ``rpavg`` will
        suffer from numerical loss of precision and can not be trusted. If
        you need accurate ``rpavg`` values, then pass in double precision
        arrays for the particle positions.

    (xyz)bin_refine_factor: integer, default is (2,2,1); typically within [1-3]
        Controls the refinement on the cell sizes. Can have up to a 20% impact
        on runtime.

    max_cells_per_dim: integer, default is 100, typical values in [50-300]
        Controls the maximum number of cells per dimension. Total number of
        cells can be up to (max_cells_per_dim)^3. Only increase if ``rpmax`` is
        too small relative to the boxsize (and increasing helps the runtime).

    copy_particles: boolean (default True)
        Boolean flag to make a copy of the particle positions
        If set to False, the particles will be re-ordered in-place

        .. versionadded:: 2.3.0

    enable_min_sep_opt: boolean (default true)
        Boolean flag to allow optimizations based on min. separation between
        pairs of cells. Here to allow for comparison studies.

        .. versionadded:: 2.3.0

    c_api_timer: boolean (default false)
        Boolean flag to measure actual time spent in the C libraries. Here
        to allow for benchmarking and scaling studies.

    isa: string (default ``fastest``)
        Controls the runtime dispatch for the instruction set to use. Options
        are: [``fastest``, ``avx512f``, ``avx``, ``sse42``, ``fallback``]

        Setting isa to ``fastest`` will pick the fastest available instruction
        set on the current computer. However, if you set ``isa`` to, say,
        ``avx`` and ``avx`` is not available on the computer, then the code
        will revert to using ``fallback`` (even though ``sse42`` might be
        available).  Unless you are benchmarking the different instruction
        sets, you should always leave ``isa`` to the default value. And if
        you *are* benchmarking, then the string supplied here gets translated
        into an ``enum`` for the instruction set defined in ``utils/defs.h``.

    weight_type: string, optional. Default: None.
        The type of weighting to apply.  One of ["pair_product", None].

    Returns
    --------

    results: Numpy structured array
        A numpy structured array containing [rpmin, rpmax, rpavg, pimax,
        npairs, weightavg] for each radial bin specified in the ``binfile``.
        If ``output_rpavg`` is not set, then ``rpavg`` will be set to 0.0 for
        all bins; similarly for ``weightavg``. ``npairs`` contains the number
        of pairs in that bin and can be used to compute :math:`\\xi(r_p, \pi)`
        by combining with (DR, RR) counts.

    api_time: float, optional
        Only returned if ``c_api_timer`` is set.  ``api_time`` measures only
        the time spent within the C library and ignores all python overhead.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.theory.DDrppi import DDrppi
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../theory/tests/", "bins")
    >>> N = 10000
    >>> boxsize = 420.0
    >>> nthreads = 4
    >>> autocorr = 1
    >>> pimax = 40.0
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> X = np.random.uniform(0, boxsize, N)
    >>> Y = np.random.uniform(0, boxsize, N)
    >>> Z = np.random.uniform(0, boxsize, N)
    >>> weights = np.ones_like(X)
    >>> results = DDrppi(autocorr, nthreads, pimax, binfile,
    ...                  X, Y, Z, weights1=weights, weight_type='pair_product', output_rpavg=True)
    >>> for r in results[519:]: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.1f}"
    ...                               " {4:10d} {5:10.6f}".format(r['rmin'], r['rmax'],
    ...                               r['rpavg'], r['pimax'], r['npairs'], r['weightavg']))
    ...                         # doctest: +NORMALIZE_WHITESPACE
     11.756000  16.753600  14.379250       40.0       1150   1.000000
     16.753600  23.875500  20.449131        1.0       2604   1.000000
     16.753600  23.875500  20.604834        2.0       2370   1.000000
     16.753600  23.875500  20.523989        3.0       2428   1.000000
     16.753600  23.875500  20.475181        4.0       2462   1.000000
     16.753600  23.875500  20.458005        5.0       2532   1.000000
     16.753600  23.875500  20.537162        6.0       2522   1.000000
     16.753600  23.875500  20.443087        7.0       2422   1.000000
     16.753600  23.875500  20.474580        8.0       2360   1.000000
     16.753600  23.875500  20.420360        9.0       2512   1.000000
     16.753600  23.875500  20.478355       10.0       2472   1.000000
     16.753600  23.875500  20.485268       11.0       2406   1.000000
     16.753600  23.875500  20.372985       12.0       2420   1.000000
     16.753600  23.875500  20.647998       13.0       2378   1.000000
     16.753600  23.875500  20.556208       14.0       2420   1.000000
     16.753600  23.875500  20.527992       15.0       2462   1.000000
     16.753600  23.875500  20.581017       16.0       2380   1.000000
     16.753600  23.875500  20.491819       17.0       2346   1.000000
     16.753600  23.875500  20.534440       18.0       2496   1.000000
     16.753600  23.875500  20.529129       19.0       2512   1.000000
     16.753600  23.875500  20.501946       20.0       2500   1.000000
     16.753600  23.875500  20.513349       21.0       2544   1.000000
     16.753600  23.875500  20.471915       22.0       2430   1.000000
     16.753600  23.875500  20.450651       23.0       2354   1.000000
     16.753600  23.875500  20.550753       24.0       2460   1.000000
     16.753600  23.875500  20.540262       25.0       2490   1.000000
     16.753600  23.875500  20.559572       26.0       2350   1.000000
     16.753600  23.875500  20.534245       27.0       2382   1.000000
     16.753600  23.875500  20.511302       28.0       2508   1.000000
     16.753600  23.875500  20.491632       29.0       2456   1.000000
     16.753600  23.875500  20.592493       30.0       2386   1.000000
     16.753600  23.875500  20.506234       31.0       2484   1.000000
     16.753600  23.875500  20.482109       32.0       2538   1.000000
     16.753600  23.875500  20.518463       33.0       2544   1.000000
     16.753600  23.875500  20.482515       34.0       2534   1.000000
     16.753600  23.875500  20.503124       35.0       2382   1.000000
     16.753600  23.875500  20.471307       36.0       2356   1.000000
     16.753600  23.875500  20.384231       37.0       2554   1.000000
     16.753600  23.875500  20.454012       38.0       2458   1.000000
     16.753600  23.875500  20.585543       39.0       2394   1.000000
     16.753600  23.875500  20.504965       40.0       2500   1.000000

    """
    try:
        from Corrfunc._countpairs import countpairs_rp_pi as DDrppi_extn
    except ImportError:
        msg = "Could not import the C extension for the 3-D "\
              "real-space pair counter."
        raise ImportError(msg)

    import numpy as np
    from Corrfunc.utils import translate_isa_string_to_enum,\
        return_file_with_rbins, convert_to_native_endian,\
        sys_pipes, process_weights
    from future.utils import bytes_to_native_str

    if not autocorr:
        if X2 is None or Y2 is None or Z2 is None:
            msg = "Must pass valid arrays for X2/Y2/Z2 for "\
                "computing cross-correlation"
            raise ValueError(msg)
    else:
        # TODO: is this needed?
        X2 = np.empty(1)
        Y2 = np.empty(1)
        Z2 = np.empty(1)

    weights1, weights2 = process_weights(weights1, weights2, X1, X2, weight_type, autocorr)

    # Ensure all input arrays are native endian
    X1, Y1, Z1, weights1, X2, Y2, Z2, weights2 = [
            convert_to_native_endian(arr, warn=True) for arr in
            [X1, Y1, Z1, weights1, X2, Y2, Z2, weights2]]

    # Passing None parameters breaks the parsing code, so avoid this
    kwargs = {}
    for k in ['weights1', 'weights2', 'weight_type', 'X2', 'Y2', 'Z2']:
        v = locals()[k]
        if v is not None:
            kwargs[k] = v

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)

    with sys_pipes():
      extn_results = DDrppi_extn(autocorr, nthreads,
                                 pimax, rbinfile,
                                 X1, Y1, Z1,
                                 periodic=periodic,
                                 verbose=verbose,
                                 boxsize=boxsize,
                                 output_rpavg=output_rpavg,
                                 xbin_refine_factor=xbin_refine_factor,
                                 ybin_refine_factor=ybin_refine_factor,
                                 zbin_refine_factor=zbin_refine_factor,
                                 max_cells_per_dim=max_cells_per_dim,
                                 copy_particles=copy_particles,
                                 enable_min_sep_opt=enable_min_sep_opt,
                                 c_api_timer=c_api_timer,
                                 isa=integer_isa, **kwargs)
    if extn_results is None:
        msg = "RuntimeError occurred"
        raise RuntimeError(msg)
    else:
        extn_results, api_time = extn_results

    if delete_after_use:
        import os
        os.remove(rbinfile)

    results_dtype = np.dtype([(bytes_to_native_str(b'rmin'), np.float),
                              (bytes_to_native_str(b'rmax'), np.float),
                              (bytes_to_native_str(b'rpavg'), np.float),
                              (bytes_to_native_str(b'pimax'), np.float),
                              (bytes_to_native_str(b'npairs'), np.uint64),
                              (bytes_to_native_str(b'weightavg'), np.float),])
    results = np.array(extn_results, dtype=results_dtype)

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import doctest
    doctest.testmod()
