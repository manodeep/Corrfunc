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


def xi(boxsize, nthreads, binfile, X, Y, Z,
       weights=None, weight_type=None, verbose=False, output_ravg=False,
       xbin_refine_factor=2, ybin_refine_factor=2,
       zbin_refine_factor=1, max_cells_per_dim=100,
       copy_particles=True, enable_min_sep_opt=True,
       c_api_timer=False, isa=r'fastest'):
    """
    Function to compute the correlation function in a
    periodic cosmological box. Pairs which are separated by less
    than the ``r`` bins (specified in ``binfile``) in 3-D real space.

    If ``weights`` are provided, the resulting correlation function
    is weighted.  The weighting scheme depends on ``weight_type``.


    .. note:: Pairs are double-counted. And if ``rmin`` is set to
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
        For string input: filename specifying the ``r`` bins for
        ``xi``. The file should contain white-space separated values
        of (rmin, rmax)  for each ``r`` wanted. The bins need to be
        contiguous and sorted in increasing order (smallest bins come first).

        For array-like input: A sequence of ``r`` values that provides the
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

    output_ravg: boolean (default false)
        Boolean flag to output the average ``r`` for each bin. Code will
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
        cells can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is
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

    weight_type: string, optional, Default: None.
        The type of weighting to apply.  One of ["pair_product", None].


    Returns
    --------

    results: Numpy structured array
        A numpy structured array containing [rmin, rmax, ravg, xi, npairs,
        weightavg] for each radial specified in the ``binfile``. If
        ``output_ravg`` is not set then ``ravg`` will be set to 0.0 for all
        bins; similarly for ``weightavg``. ``xi`` contains the correlation
        function while ``npairs`` contains the number of pairs in that bin.
        If using weights, ``xi`` will be weighted while ``npairs`` will not be.

    api_time: float, optional
        Only returned if ``c_api_timer`` is set.  ``api_time`` measures only
        the time spent within the C library and ignores all python overhead.

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
    >>> weights = np.ones_like(X)
    >>> results = xi(boxsize, nthreads, binfile, X, Y, Z, weights=weights, weight_type='pair_product', output_ravg=True)
    >>> for r in results: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f} {4:10d} {5:10.6f}"
    ...                         .format(r['rmin'], r['rmax'],
    ...                         r['ravg'], r['xi'], r['npairs'], r['weightavg']))
    ...                   # doctest: +NORMALIZE_WHITESPACE
      0.167536   0.238755   0.226592  -0.205733          4   1.000000
      0.238755   0.340251   0.289277  -0.176729         12   1.000000
      0.340251   0.484892   0.426819  -0.051829         40   1.000000
      0.484892   0.691021   0.596187  -0.131853        106   1.000000
      0.691021   0.984777   0.850100  -0.049207        336   1.000000
      0.984777   1.403410   1.225112   0.028543       1052   1.000000
      1.403410   2.000000   1.737153   0.011403       2994   1.000000
      2.000000   2.850200   2.474588   0.005405       8614   1.000000
      2.850200   4.061840   3.532018  -0.014098      24448   1.000000
      4.061840   5.788530   5.022241  -0.010784      70996   1.000000
      5.788530   8.249250   7.160648  -0.001588     207392   1.000000
      8.249250  11.756000  10.207213  -0.000323     601002   1.000000
     11.756000  16.753600  14.541171   0.000007    1740084   1.000000
     16.753600  23.875500  20.728773  -0.001595    5028058   1.000000

    """

    try:
        from Corrfunc._countpairs import countpairs_xi as xi_extn
    except ImportError:
        msg = "Could not import the C extension for the "\
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
      extn_results = xi_extn(boxsize, nthreads, rbinfile,
                             X, Y, Z,
                             verbose=verbose,
                             output_ravg=output_ravg,
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

    results_dtype = np.dtype([(bytes_to_native_str(b'rmin'), np.float64),
                              (bytes_to_native_str(b'rmax'), np.float64),
                              (bytes_to_native_str(b'ravg'), np.float64),
                              (bytes_to_native_str(b'xi'), np.float64),
                              (bytes_to_native_str(b'npairs'), np.uint64),
                              (bytes_to_native_str(b'weightavg'), np.float64)])
    results = np.array(extn_results, dtype=results_dtype)

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import doctest
    doctest.testmod()
