#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the pair counter in
``theory/DDsmu/``. This wrapper is in :py:mod:`Corrfunc.theory.DDsmu`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha', 'Nick Hand')
__all__ = ('DDsmu', )


def DDsmu(autocorr, nthreads, binfile, mu_max, nmu_bins, X1, Y1, Z1, weights1=None,
          periodic=True, X2=None, Y2=None, Z2=None, weights2=None,
          verbose=False, boxsize=0.0, output_savg=False,
          fast_divide_and_NR_steps=0,
          xbin_refine_factor=2, ybin_refine_factor=2,
          zbin_refine_factor=1, max_cells_per_dim=100,
          c_api_timer=False, isa=r'fastest', weight_type=None):
    """
    Calculate the 2-D pair-counts corresponding to the redshift-space 
    correlation function, :math:`\\xi(s, \mu)` Pairs which are separated
    by less than the ``s`` bins (specified in ``binfile``) in 3-D, and
    less than ``s*mu_max`` in the Z-dimension are counted.

    If ``weights`` are provided, the resulting pair counts are weighted.  The
    weighting scheme depends on ``weight_type``.


    .. note:: This module only returns pair counts and not the actual
       correlation function :math:`\\xi(s, \mu)`. See the
       utilities :py:mod:`Corrfunc.utils.convert_3d_counts_to_cf` 
       for computing :math:`\\xi(s, \mu)` from the pair counts.

    .. versionadded:: 2.1.0

    Parameters
    ----------

    autocorr: boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.

    nthreads: integer
        The number of OpenMP threads to use. Has no effect if OpenMP was not
        enabled during library compilation.

    binfile: string or an list/array of floats
        For string input: filename specifying the ``s`` bins for
        ``DDsmu_mocks``. The file should contain white-space separated values
        of (smin, smax) specifying each ``s`` bin wanted. The bins
        need to be contiguous and sorted in increasing order (smallest bins
        come first).

        For array-like input: A sequence of ``s`` values that provides the
        bin-edges. For example,
        ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
        input specifying **14** (logarithmic) bins between 0.1 and 10.0. This
        array does not need to be sorted.         

    mu_max: double. Must be in range (0.0, 1.0]
        A double-precision value for the maximum cosine of the angular 
        separation from the line of sight (LOS). Here, LOS is taken to be 
        along the Z direction. 

        Note: Only pairs with :math:`0 <= \cos(\\theta_{LOS}) < \mu_{max}`
        are counted (no equality).

    nmu_bins: int
        The number of linear ``mu`` bins, with the bins ranging from
        from (0, :math:`\mu_{max}`)

    X1/Y1/Z1 : array-like, real (float/double)
        The array of X/Y/Z positions for the first set of points.
        Calculations are done in the precision of the supplied arrays.

    weights1 : array-like, real (float/double), shape (n_particles,) or \
        (n_weights_per_particle,n_particles), optional
        Weights for computing a weighted pair count.

    weight_type : str, optional
        The type of pair weighting to apply.
        Options: "pair_product", None; Default: None.

    periodic : boolean
        Boolean flag to indicate periodic boundary conditions.

    X2/Y2/Z2 : array-like, real (float/double)
        Array of XYZ positions for the second set of points. *Must* be the same
        precision as the X1/Y1/Z1 arrays. Only required when ``autocorr==0``.

    weights2 : array-like, real (float/double), shape (n_particles,) or \
        (n_weights_per_particle,n_particles), optional
        Weights for computing a weighted pair count.

    verbose : boolean (default false)
        Boolean flag to control output of informational messages

    boxsize : double
        The side-length of the cube in the cosmological simulation.
        Present to facilitate exact calculations for periodic wrapping.
        If boxsize is not supplied, then the wrapping is done based on
        the maximum difference within each dimension of the X/Y/Z arrays.

    output_savg : boolean (default false)
        Boolean flag to output the average ``s`` for each bin. Code will
        run slower if you set this flag. Also, note, if you are calculating
        in single-precision, ``s`` will suffer from numerical loss of
        precision and can not be trusted. If you need accurate ``s``
        values, then pass in double precision arrays for the particle positions.

    fast_divide_and_NR_steps: integer (default 0)
        Replaces the division in ``AVX`` implementation with an approximate
        reciprocal, followed by ``fast_divide_and_NR_steps`` of Newton-Raphson.
        Can improve runtime by ~15-20% on older computers. Value of 0 uses
        the standard division operation.
    
    (xyz)bin_refine_factor: integer (default (2,2,1) typical values in [1-3])
        Controls the refinement on the cell sizes. Can have up to a 20% impact
        on runtime.

    max_cells_per_dim: integer (default 100, typical values in [50-300])
        Controls the maximum number of cells per dimension. Total number of 
        cells can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is
        too small relative to the boxsize (and increasing helps the runtime).

    c_api_timer : boolean (default false)
        Boolean flag to measure actual time spent in the C libraries. Here
        to allow for benchmarking and scaling studies.

    isa : integer (default -1)
      Controls the runtime dispatch for the instruction set to use. Possible
      options are: [-1, AVX, SSE42, FALLBACK]

      Setting isa to -1 will pick the fastest available instruction
      set on the current computer. However, if you set ``isa`` to, say,
      ``AVX`` and ``AVX`` is not available on the computer, then the code will
      revert to using ``FALLBACK`` (even though ``SSE42`` might be available).

      Unless you are benchmarking the different instruction sets, you should
      always leave ``isa`` to the default value. And if you *are* benchmarking,
      then the integer values correspond to the ``enum`` for the instruction set
      defined in ``utils/defs.h``.

    Returns
    --------
    results : A python list
        A python list containing ``nmu_bins`` of [smin, smax, savg, mu_max, npairs, weightavg]
        for each spatial bin specified in the ``binfile``. There will be a total of ``nmu_bins``
        ranging from [0, ``mu_max``) *per* spatial bin. If ``output_savg`` is not set, then ``savg``
        will be set to 0.0 for all bins; similarly for ``weight_avg``. ``npairs``
        contains the number of pairs in that bin.

    time : if ``c_api_timer`` is set, then the return value contains the time spent
        in the API; otherwise time is set to 0.0

    Example
    -------
    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.theory.DDsmu import DDsmu
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../theory/tests/", "bins")
    >>> N = 10000
    >>> boxsize = 420.0
    >>> nthreads = 4
    >>> autocorr = 1
    >>> mu_max = 1.0
    >>> seed = 42
    >>> nmu_bins = 10
    >>> np.random.seed(seed)
    >>> X = np.random.uniform(0, boxsize, N)
    >>> Y = np.random.uniform(0, boxsize, N)
    >>> Z = np.random.uniform(0, boxsize, N)
    >>> weights = np.ones_like(X)
    >>> results = DDsmu(autocorr, nthreads, binfile, mu_max, nmu_bins,
    ...                  X, Y, Z, weights1=weights, weight_type='pair_product', output_savg=True)
    >>> for r in results[100:]: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.1f}"
    ...                               " {4:10d} {5:10.6f}".format(r['smin'], r['smax'],
    ...                               r['savg'], r['mu_max'], r['npairs'], r['weightavg']))
    ...                         # doctest: +NORMALIZE_WHITESPACE
     5.788530   8.249250   7.148213        0.1        230   1.000000
     5.788530   8.249250   7.157218        0.2        236   1.000000
     5.788530   8.249250   7.165338        0.3        208   1.000000
     5.788530   8.249250   7.079905        0.4        252   1.000000
     5.788530   8.249250   7.251661        0.5        184   1.000000
     5.788530   8.249250   7.118536        0.6        222   1.000000
     5.788530   8.249250   7.083466        0.7        238   1.000000
     5.788530   8.249250   7.198184        0.8        170   1.000000
     5.788530   8.249250   7.127409        0.9        208   1.000000
     5.788530   8.249250   6.973090        1.0        206   1.000000
     8.249250  11.756000  10.149183        0.1        592   1.000000
     8.249250  11.756000  10.213009        0.2        634   1.000000
     8.249250  11.756000  10.192220        0.3        532   1.000000
     8.249250  11.756000  10.246931        0.4        544   1.000000
     8.249250  11.756000  10.102675        0.5        530   1.000000
     8.249250  11.756000  10.276180        0.6        644   1.000000
     8.249250  11.756000  10.251264        0.7        666   1.000000
     8.249250  11.756000  10.138399        0.8        680   1.000000
     8.249250  11.756000  10.191916        0.9        566   1.000000
     8.249250  11.756000  10.243229        1.0        608   1.000000
    11.756000  16.753600  14.552776        0.1       1734   1.000000
    11.756000  16.753600  14.579991        0.2       1806   1.000000
    11.756000  16.753600  14.599611        0.3       1802   1.000000
    11.756000  16.753600  14.471100        0.4       1820   1.000000
    11.756000  16.753600  14.480192        0.5       1740   1.000000
    11.756000  16.753600  14.493679        0.6       1746   1.000000
    11.756000  16.753600  14.547713        0.7       1722   1.000000
    11.756000  16.753600  14.465390        0.8       1750   1.000000
    11.756000  16.753600  14.547465        0.9       1798   1.000000
    11.756000  16.753600  14.440975        1.0       1828   1.000000
    16.753600  23.875500  20.720406        0.1       5094   1.000000
    16.753600  23.875500  20.735403        0.2       5004   1.000000
    16.753600  23.875500  20.721069        0.3       5172   1.000000
    16.753600  23.875500  20.723648        0.4       5014   1.000000
    16.753600  23.875500  20.650621        0.5       5094   1.000000
    16.753600  23.875500  20.688135        0.6       5076   1.000000
    16.753600  23.875500  20.735691        0.7       4910   1.000000
    16.753600  23.875500  20.714097        0.8       4864   1.000000
    16.753600  23.875500  20.751836        0.9       4954   1.000000
    16.753600  23.875500  20.721183        1.0       5070   1.000000
    """
    try:
        from Corrfunc._countpairs import countpairs_s_mu as DDsmu_extn
    except ImportError:
        msg = "Could not import the C extension for the 3-D "\
              "redshift-space pair counter."
        raise ImportError(msg)

    import numpy as np
    from Corrfunc.utils import translate_isa_string_to_enum,\
        return_file_with_rbins, sys_pipes
    from future.utils import bytes_to_native_str

    # Broadcast scalar weights to arrays
    if weights1 is not None:
        weights1 = np.atleast_1d(weights1)
    if weights2 is not None:
        weights2 = np.atleast_1d(weights2)

    # Check if mu_max is scalar
    if not np.isscalar(mu_max):
        msg = "The parameter `mu_max` = {0}, has size = {1}. "\
              "The code is expecting a scalar quantity (and not "\
              "not a list, array)".format(mu_max, np.size(mu_max))
        raise TypeError(msg)

    # Check that mu_max is within (0.0, 1.0]
    if mu_max <= 0.0 or mu_max > 1.0:
        msg = "The parameter `mu_max` = {0}, is the max. of cosine of an "
        "angle and should be within (0.0, 1.0]".format(mu_max)
        raise ValueError(msg)
        
    if not autocorr:
        if X2 is None or Y2 is None or Z2 is None:
            msg = "Must pass valid arrays for X2/Y2/Z2 for "\
                "computing cross-correlation"
            raise ValueError(msg)

        # If only one set of points has weights, set the other to uniform weights
        if weights1 is None and weights2 is not None:
            weights1 = np.ones_like(weights2)
        if weights2 is None and weights1 is not None:
            weights2 = np.ones_like(weights1)

    else:
        X2 = np.empty(1)
        Y2 = np.empty(1)
        Z2 = np.empty(1)

    # Passing None parameters breaks the parsing code, so avoid this
    kwargs = {}
    for k in ['weights1', 'weights2', 'weight_type', 'X2', 'Y2', 'Z2']:
        v = locals()[k]
        if v is not None:
            kwargs[k] = v

    integer_isa = translate_isa_string_to_enum(isa)
    sbinfile, delete_after_use = return_file_with_rbins(binfile)
    with sys_pipes():
        extn_results = DDsmu_extn(autocorr, nthreads,
                                  sbinfile,
                                  mu_max, nmu_bins,
                                  X1, Y1, Z1,
                                  periodic=periodic,
                                  verbose=verbose,
                                  boxsize=boxsize,
                                  output_savg=output_savg,
                                  fast_divide_and_NR_steps=fast_divide_and_NR_steps,
                                  xbin_refine_factor=xbin_refine_factor,
                                  ybin_refine_factor=ybin_refine_factor,
                                  zbin_refine_factor=zbin_refine_factor,
                                  max_cells_per_dim=max_cells_per_dim,
                                  c_api_timer=c_api_timer,
                                  isa=integer_isa, **kwargs)
    if extn_results is None:
        msg = "RuntimeError occurred"
        raise RuntimeError(msg)
    else:
        extn_results, api_time = extn_results

    if delete_after_use:
        import os
        os.remove(sbinfile)

    results_dtype = np.dtype([(bytes_to_native_str(b'smin'), np.float),
                              (bytes_to_native_str(b'smax'), np.float),
                              (bytes_to_native_str(b'savg'), np.float),
                              (bytes_to_native_str(b'mu_max'), np.float),
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
