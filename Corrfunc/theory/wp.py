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
__all__ = ('wp', 'find_fastest_wp_bin_refs', )


def find_fastest_wp_bin_refs(boxsize, pimax, nthreads, binfile, X, Y, Z,
                             verbose=False, output_rpavg=False,
                             max_cells_per_dim=100,
                             isa=r'fastest',
                             maxbinref=3, nrepeats=3,
                             return_runtimes=False):

    """
    Finds the combination of ``bin refine factors`` that produces the
    fastest computation for the given dataset and ``rp`` limits.

    Parameters
    -----------

    boxsize: double
        A double-precision value for the boxsize of the simulation
        in same units as the particle positions and the ``rp`` bins.

    pimax: double
        A double-precision value for the maximum separation along
        the Z-dimension.

        Note: Only pairs with ``0 <= dz < pimax`` are counted (no equality).

    nthreads: integer
        Number of threads to use.

    binfile: string or an list/array of floats
        For string input: filename specifying the ``rp`` bins for
        ``wp``. The file should contain white-space separated values
        of (rpmin, rpmax)  for each ``rp`` wanted. The bins need to be
        contiguous and sorted in increasing order (smallest bins come first).

        For array-like input: A sequence of ``rp`` values that provides the
        bin-edges. For example,
        ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
        input specifying **14** (logarithmic) bins between 0.1 and 10.0. This
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
        run slower if you set this flag.

        Note: If you are calculating in single-precision, ``rpavg`` will
        suffer from numerical loss of precision and can not be trusted. If
        you need accurate ``rpavg`` values, then pass in double precision
        arrays for the particle positions.

    max_cells_per_dim: integer, default is 100, typical values in [50-300]
        Controls the maximum number of cells per dimension. Total number of
        cells can be up to (max_cells_per_dim)^3. Only increase if ``rpmax`` is
        too small relative to the boxsize (and increasing helps the runtime).

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

    maxbinref: integer (default 3)
        The maximum ``bin refine factor`` to use along each dimension. From
        experience, values larger than 3 do not improve ``wp`` runtime.

        Runtime of module scales as ``maxbinref^3``, so change the value of
        ``maxbinref`` with caution.

    nrepeats: integer (default 3)
        Number of times to repeat the timing for an individual run. Accounts
        for the dispersion in runtimes on computers with multiple user
        processes.

    return_runtimes: boolean (default ``false``)
        If set, also returns the array of runtimes.

    Returns
    --------
    (nx, ny, nz) : tuple of integers
        The combination of ``bin refine factors`` along each dimension that
        produces the fastest code.

    runtimes : numpy structured array

        Only returned if ``return_runtimes`` is set, then the return value is a tuple
        containing ((nx, ny, nz), runtimes). ``runtimes`` is a ``numpy``
        structured array containing the fields, [``nx``, ``ny``, ``nz``,
        ``avg_runtime``, ``sigma_time``]. Here, ``avg_runtime`` is the
        average time, measured over ``nrepeats`` invocations, spent in
        the python extension. ``sigma_time`` is the dispersion of the
        run times across those ``nrepeats`` invocations.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.io import read_catalog
    >>> from Corrfunc.theory.wp import find_fastest_wp_bin_refs
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../theory/tests/", "bins")
    >>> X, Y, Z = read_catalog(return_dtype=np.float32)
    >>> boxsize = 420.0
    >>> pimax = 40.0
    >>> nthreads = 4
    >>> verbose = 1
    >>> best, _ = find_fastest_wp_bin_refs(boxsize, pimax, nthreads, binfile,
    ...                                    X, Y, Z, maxbinref=2, nrepeats=3,
    ...                                    verbose=verbose,
    ...                                    return_runtimes=True)
    >>> print(best) # doctest:+SKIP
    (2, 2, 1)

    .. note:: Since the result might change depending on the computer, doctest
        is skipped for this function.


    """
    try:
        from Corrfunc._countpairs import countpairs_wp as wp_extn
    except ImportError:
        msg = "Could not import the C extension for the projected "\
              "correlation function."
        raise ImportError(msg)

    from Corrfunc.utils import translate_isa_string_to_enum,\
        return_file_with_rbins

    import itertools
    import numpy as np
    from future.utils import bytes_to_native_str
    import time

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    bin_refs = np.arange(1, maxbinref + 1)
    bin_ref_perms = itertools.product(bin_refs, bin_refs, bin_refs)
    dtype = np.dtype([(bytes_to_native_str(b'nx'), np.int),
                      (bytes_to_native_str(b'ny'), np.int),
                      (bytes_to_native_str(b'nz'), np.int),
                      (bytes_to_native_str(b'avg_time'), np.float64),
                      (bytes_to_native_str(b'sigma_time'), np.float64)])
    all_runtimes = np.zeros(maxbinref**3, dtype=dtype)
    all_runtimes[:] = np.inf

    for ii, (nx, ny, nz) in enumerate(bin_ref_perms):
        total_runtime = 0.0
        total_sqr_runtime = 0.0

        for _ in range(nrepeats):
            t0 = time.time()
            extn_results, _, _ = wp_extn(boxsize, pimax, nthreads,
                                         rbinfile,
                                         X, Y, Z,
                                         verbose=verbose,
                                         output_rpavg=output_rpavg,
                                         xbin_refine_factor=nx,
                                         ybin_refine_factor=ny,
                                         zbin_refine_factor=nz,
                                         max_cells_per_dim=max_cells_per_dim,
                                         isa=integer_isa)
            t1 = time.time()

            if extn_results is None:
                msg = "RuntimeError occurred with perms = ({0}, {1}, {2})".\
                                                          format(nx, ny, nz)
                print(msg)
                print("Continuing...")
                continue

            dt = (t1 - t0)
            total_runtime += dt
            total_sqr_runtime += dt*dt

        avg_runtime = total_runtime/nrepeats

        # variance = E(X^2) - E^2(X)
        # disp = sqrt(variance)
        runtime_disp = np.sqrt(total_sqr_runtime/nrepeats -
                               avg_runtime*avg_runtime)

        all_runtimes[ii]['nx'] = nx
        all_runtimes[ii]['ny'] = ny
        all_runtimes[ii]['nz'] = nz
        all_runtimes[ii]['avg_time'] = avg_runtime
        all_runtimes[ii]['sigma_time'] = runtime_disp

    if delete_after_use:
        import os
        os.remove(rbinfile)

    all_runtimes.sort(order=('avg_time', 'sigma_time'))
    results = (all_runtimes[0]['nx'],
               all_runtimes[0]['ny'],
               all_runtimes[0]['nz'])

    optional_returns = return_runtimes
    if not optional_returns:
        ret = results
    else:
        ret = (results, )
        if return_runtimes:
            ret += (all_runtimes, )

    return ret


def _convert_cell_timer(cell_time_lst):
    """
    Converts a the cell timings list returned by the python extensions
    into a more user-friendly numpy structured array.

    The fields correspond to the C ``struct api_cell_timings`` defined
    in ``utils/defs.h``.

    Returns:
    --------

    cell_times : numpy structured array
        The following fields are present in the ``cell_times``:

        N1 -> number of particles in cell 1
        N2 -> number of particles in cell 2
        time_in_ns -> time taken to compute all pairs between two cells
                      (cellidx1, cellidx2)
        cellidx1, cellidx2 -> the 1-D index for the two cells
        tid -> thread-id that computed the pairs (identically 0 for
               serial/single-threaded runs)


    """

    import numpy as np
    from future.utils import bytes_to_native_str

    dtype = np.dtype([(bytes_to_native_str(b'N1'), np.int64),
                      (bytes_to_native_str(b'N2'), np.int64),
                      (bytes_to_native_str(b'time_in_ns'), np.int64),
                      (bytes_to_native_str(b'cellidx1'), np.int32),
                      (bytes_to_native_str(b'cellidx2'), np.int32),
                      (bytes_to_native_str(b'tid'), np.int32)])

    cell_times = np.array(cell_time_lst, dtype=dtype)
    return cell_times


def wp(boxsize, pimax, nthreads, binfile, X, Y, Z,
       weights=None, weight_type=None, verbose=False, output_rpavg=False,
       xbin_refine_factor=2, ybin_refine_factor=2,
       zbin_refine_factor=1, max_cells_per_dim=100,
       copy_particles=True, enable_min_sep_opt=True,
       c_api_timer=False, c_cell_timer=False, isa='fastest'):
    """
    Function to compute the projected correlation function in a
    periodic cosmological box. Pairs which are separated by less
    than the ``rp`` bins (specified in ``binfile``) in the
    X-Y plane, and less than ``pimax`` in the Z-dimension are
    counted.

    If ``weights`` are provided, the resulting correlation function
    is weighted.  The weighting scheme depends on ``weight_type``.


    .. note:: Pairs are double-counted. And if ``rpmin`` is set to
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
        the Z-dimension.

        Note: Only pairs with ``0 <= dz < pimax`` are counted (no equality).

    nthreads: integer
        Number of threads to use.

    binfile: string or an list/array of floats
        For string input: filename specifying the ``rp`` bins for
        ``wp``. The file should contain white-space separated values
        of (rpmin, rpmax)  for each ``rp`` wanted. The bins need to be
        contiguous and sorted in increasing order (smallest bins come first).

        For array-like input: A sequence of ``rp`` values that provides the
        bin-edges. For example,
        ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
        input specifying **14** (logarithmic) bins between 0.1 and 10.0. This
        array does not need to be sorted.

    X/Y/Z: arraytype, real (float/double)
        Particle positions in the 3 axes. Must be within [0, boxsize]
        and specified in the same units as ``rp_bins`` and boxsize. All
        3 arrays must be of the same floating-point type.

        Calculations will be done in the same precision as these arrays,
        i.e., calculations will be in floating point if XYZ are single
        precision arrays (C float type); or in double-precision if XYZ
        are double precision arrays (C double type).

    weights: array_like, real (float/double), optional
        A scalar, or an array of weights of shape (n_weights, n_positions) or
        (n_positions,). ``weight_type`` specifies how these weights are used;
        results are returned in the ``weightavg`` field.

    verbose: boolean (default false)
        Boolean flag to control output of informational messages

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

    c_cell_timer : boolean (default false)
        Boolean flag to measure actual time spent **per cell-pair** within the
        C libraries. A very detailed timer that stores information about the
        number of particles in each cell, the thread id that processed that
        cell-pair and the amount of time in nano-seconds taken to process that
        cell pair. This timer can be used to study the instruction set
        efficiency, and load-balancing of the code.

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

    weight_type: string, optional.  Default: None.
        The type of weighting to apply.  One of ["pair_product", None].

    Returns
    --------

    results: Numpy structured array
        A numpy structured array containing [rpmin, rpmax, rpavg, wp, npairs,
        weightavg] for each radial specified in the ``binfile``. If
        ``output_rpavg`` is not set then ``rpavg`` will be set to 0.0 for all
        bins; similarly for ``weightavg``. ``wp`` contains the projected
        correlation function while ``npairs`` contains the number of unique
        pairs in that bin.  If using weights, ``wp`` will be weighted while
        ``npairs`` will not be.

    api_time: float, optional
        Only returned if ``c_api_timer`` is set.  ``api_time`` measures only
        the time spent within the C library and ignores all python overhead.

    cell_time: list, optional
        Only returned if ``c_cell_timer`` is set. Contains
        detailed stats about each cell-pair visited during pair-counting,
        viz., number of particles in each of the cells in the pair, 1-D
        cell-indices for each cell in the pair, time (in nano-seconds) to
        process the pair and the thread-id for the thread that processed that
        cell-pair.

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
    >>> results = wp(boxsize, pimax, nthreads, binfile, X, Y, Z, weights=np.ones_like(X), weight_type='pair_product')
    >>> for r in results:
    ...     print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f} {4:10d} {5:10.6f}".
    ...           format(r['rmin'], r['rmax'],
    ...           r['rpavg'], r['wp'], r['npairs'], r['weightavg']))
    ...           # doctest: +NORMALIZE_WHITESPACE
      0.167536   0.238755   0.000000  66.717143         18   1.000000
      0.238755   0.340251   0.000000 -15.786045         16   1.000000
      0.340251   0.484892   0.000000   2.998470         42   1.000000
      0.484892   0.691021   0.000000 -15.779885         66   1.000000
      0.691021   0.984777   0.000000 -11.966728        142   1.000000
      0.984777   1.403410   0.000000  -9.699906        298   1.000000
      1.403410   2.000000   0.000000 -11.698771        588   1.000000
      2.000000   2.850200   0.000000   3.848375       1466   1.000000
      2.850200   4.061840   0.000000  -0.921452       2808   1.000000
      4.061840   5.788530   0.000000   0.454851       5802   1.000000
      5.788530   8.249250   0.000000   1.428344      11926   1.000000
      8.249250  11.756000   0.000000  -1.067885      23478   1.000000
     11.756000  16.753600   0.000000  -0.553319      47994   1.000000
     16.753600  23.875500   0.000000  -0.086433      98042   1.000000
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
        return_file_with_rbins, convert_to_native_endian,\
        sys_pipes, process_weights

    weights, _ = process_weights(weights, None, X, None, weight_type, autocorr=True)

    # Ensure all input arrays are native endian
    X, Y, Z, weights = [convert_to_native_endian(arr, warn=True)
                        for arr in [X, Y, Z, weights]]

    # Passing None parameters breaks the parsing code, so avoid this
    kwargs = {}
    for k in ['weights', 'weight_type']:
        v = locals()[k]
        if v is not None:
            kwargs[k] = v

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    with sys_pipes():
      extn_results = wp_extn(boxsize, pimax, nthreads,
                             rbinfile,
                             X, Y, Z,
                             verbose=verbose,
                             output_rpavg=output_rpavg,
                             xbin_refine_factor=xbin_refine_factor,
                             ybin_refine_factor=ybin_refine_factor,
                             zbin_refine_factor=zbin_refine_factor,
                             max_cells_per_dim=max_cells_per_dim,
                             copy_particles=copy_particles,
                             enable_min_sep_opt=enable_min_sep_opt,
                             c_api_timer=c_api_timer,
                             c_cell_timer=c_cell_timer,
                             isa=integer_isa, **kwargs)
    if extn_results is None:
        msg = "RuntimeError occurred"
        raise RuntimeError(msg)
    else:
        extn_results, api_time, cell_time = extn_results

    if delete_after_use:
        import os
        os.remove(rbinfile)

    results_dtype = np.dtype([(bytes_to_native_str(b'rmin'), np.float64),
                              (bytes_to_native_str(b'rmax'), np.float64),
                              (bytes_to_native_str(b'rpavg'), np.float64),
                              (bytes_to_native_str(b'wp'), np.float64),
                              (bytes_to_native_str(b'npairs'), np.uint64),
                              (bytes_to_native_str(b'weightavg'), np.float64)])
    results = np.array(extn_results, dtype=results_dtype)

    # A better solution for returning multiple values based on
    # input parameter. Lifted straight from numpy.unique -- MS 10/26/2016
    optional_returns = c_api_timer or c_cell_timer
    if not optional_returns:
        ret = results
    else:
        ret = (results, )

        if c_api_timer:
            ret += (api_time, )

        if c_cell_timer:
            # Convert to numpy structured array
            np_cell_time = _convert_cell_timer(cell_time)
            ret += (np_cell_time, )

    return ret


if __name__ == '__main__':
    import doctest
    doctest.testmod()
