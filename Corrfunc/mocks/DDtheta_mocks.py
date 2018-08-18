#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the angular correlation function
:math:`\\omega(\\theta)`. Corresponding C routines are in 
``mocks/DDtheta_mocks/``, while the python interface is 
:py:mod:`Corrfunc.mocks.DDtheta_mocks`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DDtheta_mocks',)


def DDtheta_mocks(autocorr, nthreads, binfile,
                  RA1, DEC1, weights1=None,
                  RA2=None, DEC2=None, weights2=None,
                  link_in_dec=True, link_in_ra=True,
                  verbose=False, output_thetaavg=False,
                  fast_acos=False, ra_refine_factor=2,
                  dec_refine_factor=2, max_cells_per_dim=100,
                  c_api_timer=False, isa=r'fastest', weight_type=None):
    """
    Function to compute the angular correlation function for points on
    the sky (i.e., mock catalogs or observed galaxies).

    Returns a numpy structured array containing the pair counts for the
    specified angular bins.
    
    If ``weights`` are provided, the resulting pair counts are weighted.  The
    weighting scheme depends on ``weight_type``.


    .. note:: This module only returns pair counts and not the actual
       correlation function :math:`\\omega(\theta)`. See 
       :py:mod:`Corrfunc.utils.convert_3d_counts_to_cf` for computing 
       :math:`\\omega(\theta)` from the pair counts returned.


    Parameters
    -----------

    autocorr : boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.

    nthreads : integer
       Number of threads to use.

    binfile: string or an list/array of floats. Units: degrees.
        For string input: filename specifying the ``theta`` bins for
        ``DDtheta_mocks``. The file should contain white-space separated values
        of (thetamin, thetamax)  for each ``theta`` wanted. The bins need to be
        contiguous and sorted in increasing order (smallest bins come first).

        For array-like input: A sequence of ``theta`` values that provides the
        bin-edges. For example,
        ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
        input specifying **14** (logarithmic) bins between 0.1 and 10.0
        degrees. This array does not need to be sorted.         

    RA1 : array-like, real (float/double)
        The array of Right Ascensions for the first set of points. RA's
        are expected to be in [0.0, 360.0], but the code will try to fix cases
        where the RA's are in [-180, 180.0]. For peace of mind, always supply
        RA's in [0.0, 360.0].

        Calculations are done in the precision of the supplied arrays.

    DEC1 : array-like, real (float/double)
       Array of Declinations for the first set of points. DEC's are expected
       to be in the [-90.0, 90.0], but the code will try to fix cases where
       the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
       DEC's in [-90.0, 90.0].
       Must be of same precision type as RA1.
       
    weights1 : array_like, real (float/double), optional
       A scalar, or an array of weights of shape (n_weights, n_positions) or 
       (n_positions,). `weight_type` specifies how these weights are used; 
       results are returned in the `weightavg` field.  If only one of weights1 
       and weights2 is specified, the other will be set to uniform weights.

    RA2 : array-like, real (float/double)
       The array of Right Ascensions for the second set of points. RA's
       are expected to be in [0.0, 360.0], but the code will try to fix cases
       where the RA's are in [-180, 180.0]. For peace of mind, always supply
       RA's in [0.0, 360.0].
       Must be of same precision type as RA1/DEC1.

    DEC2 : array-like, real (float/double)
       Array of Declinations for the second set of points. DEC's are expected
       to be in the [-90.0, 90.0], but the code will try to fix cases where
       the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
       DEC's in [-90.0, 90.0].
       Must be of same precision type as RA1/DEC1.
       
    weights2 : array-like, real (float/double), optional
       Same as weights1, but for the second set of positions

    link_in_dec : boolean (default True)
       Boolean flag to create lattice in Declination. Code runs faster with
       this option. However, if the angular separations are too small, then
       linking in declination might produce incorrect results. When running
       for the first time, check your results by comparing with the output
       of the code for ``link_in_dec=False`` and ``link_in_ra=False``.

    link_in_ra : boolean (default True)
       Boolean flag to create lattice in Right Ascension. Setting this option
       implies ``link_in_dec=True``. Similar considerations as ``link_in_dec``
       described above.

       If you disable both ``link_in_dec`` and ``link_in_ra``, then
       the code reduces to a brute-force pair counter. No lattices are created
       at all. For very small angular separations, the brute-force method 
       might be the most numerically stable method.

    verbose : boolean (default false)
       Boolean flag to control output of informational messages

    output_thetaavg : boolean (default false)
       Boolean flag to output the average ``\theta`` for each bin. Code will
       run slower if you set this flag. 

       If you are calculating in single-precision, ``thetaavg`` will 
       suffer from numerical loss of precision and can not be trusted. If you 
       need accurate ``thetaavg`` values, then pass in double precision arrays 
       for ``RA/DEC``.

       Code will run significantly slower if you enable this option.
       Use the keyword ``fast_acos`` if you can tolerate some loss of 
       precision.

    fast_acos : boolean (default false)
       Flag to use numerical approximation for the ``arccos`` - gives better
       performance at the expense of some precision. Relevant only if
       ``output_thetaavg==True``.

       Developers: Two versions already coded up in ``utils/fast_acos.h``,
       so you can choose the version you want. There are also notes on how
       to implement faster (and less accurate) functions, particularly relevant
       if you know your ``theta`` range is limited. If you implement a new
       version, then you will have to reinstall the entire Corrfunc package.

       Note: Tests will fail if you run the tests with``fast_acos=True``.

    (radec)_refine_factor : integer, default is (2,2); typically within [1-3]
       Controls the refinement on the cell sizes. Can have up to a 20% impact
       on runtime. 

       Only two refine factors are to be specified and these
       correspond to ``ra`` and ``dec`` (rather, than the usual three of
       ``(xyz)bin_refine_factor`` for all other correlation functions).

    max_cells_per_dim : integer, default is 100, typical values in [50-300]
       Controls the maximum number of cells per dimension. Total number of
       cells can be up to (max_cells_per_dim)^3. Only increase if ``thetamax``
       is too small relative to the boxsize (and increasing helps the runtime).

    c_api_timer : boolean (default false)
       Boolean flag to measure actual time spent in the C libraries. Here
       to allow for benchmarking and scaling studies.

    isa : string (default ``fastest``)
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

    results : Numpy structured array
       A numpy structured array containing [thetamin, thetamax, thetaavg,
       npairs, weightavg] for each angular bin specified in the ``binfile``. If
       ``output_thetaavg`` is not set then ``thetavg`` will be set to 0.0 for
       all bins; similarly for
       ``weightavg``. ``npairs`` contains the number of pairs in that bin.

    api_time : float, optional
       Only returned if ``c_api_timer`` is set.  ``api_time`` measures only the time
       spent within the C library and ignores all python overhead.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> import time
    >>> from math import pi
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../mocks/tests/", "angular_bins")
    >>> N = 100000
    >>> nthreads = 4
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> RA = np.random.uniform(0.0, 2.0*pi, N)*180.0/pi
    >>> cos_theta = np.random.uniform(-1.0, 1.0, N)
    >>> DEC = 90.0 - np.arccos(cos_theta)*180.0/pi
    >>> weights = np.ones_like(RA)
    >>> autocorr = 1
    >>> for isa in ['AVX', 'SSE42', 'FALLBACK']:
    ...     for link_in_dec in [False, True]:
    ...         for link_in_ra in [False, True]:
    ...             results = DDtheta_mocks(autocorr, nthreads, binfile,
    ...                         RA, DEC, output_thetaavg=True, 
    ...                         weights1=weights, weight_type='pair_product',
    ...                         link_in_dec=link_in_dec, link_in_ra=link_in_ra, 
    ...                         isa=isa, verbose=True)
    >>> for r in results: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10d} {4:10.6f}".
    ...                         format(r['thetamin'], r['thetamax'],
    ...                         r['thetaavg'], r['npairs'], r['weightavg']))
    ...                         # doctest: +NORMALIZE_WHITESPACE
      0.010000   0.014125   0.012272         62   1.000000
      0.014125   0.019953   0.016978        172   1.000000
      0.019953   0.028184   0.024380        298   1.000000
      0.028184   0.039811   0.034321        598   1.000000
      0.039811   0.056234   0.048535       1164   1.000000
      0.056234   0.079433   0.068385       2438   1.000000
      0.079433   0.112202   0.096631       4658   1.000000
      0.112202   0.158489   0.136834       9414   1.000000
      0.158489   0.223872   0.192967      19098   1.000000
      0.223872   0.316228   0.272673      37848   1.000000
      0.316228   0.446684   0.385344      75520   1.000000
      0.446684   0.630957   0.543973     150938   1.000000
      0.630957   0.891251   0.768406     301854   1.000000
      0.891251   1.258925   1.085273     599896   1.000000
      1.258925   1.778279   1.533461    1200238   1.000000
      1.778279   2.511886   2.166009    2396338   1.000000
      2.511886   3.548134   3.059159    4775162   1.000000
      3.548134   5.011872   4.321445    9532582   1.000000
      5.011872   7.079458   6.104214   19001930   1.000000
      7.079458  10.000000   8.622400   37842502   1.000000

    """

    try:
        from Corrfunc._countpairs_mocks import countpairs_theta_mocks as\
            DDtheta_mocks_extn
    except ImportError:
        msg = "Could not import the C extension for the angular "\
              "correlation function for mocks."
        raise ImportError(msg)

    import numpy as np
    from warnings import warn
    from Corrfunc.utils import translate_isa_string_to_enum, fix_ra_dec,\
        return_file_with_rbins, convert_to_native_endian,\
        is_native_endian, sys_pipes
    from future.utils import bytes_to_native_str
    
    # Broadcast scalar weights to arrays
    if weights1 is not None:
        weights1 = np.atleast_1d(weights1)
    if weights2 is not None:
        weights2 = np.atleast_1d(weights2)

    if autocorr == 0:
        if RA2 is None or DEC2 is None:
            msg = "Must pass valid arrays for RA2/DEC2 for "\
                  "computing cross-correlation"
            raise ValueError(msg)
            
        # If only one set of points has weights, set the other to uniform weights
        if weights1 is None and weights2 is not None:
            weights1 = np.ones_like(weights2)
        if weights2 is None and weights1 is not None:
            weights2 = np.ones_like(weights1)
    else:
        RA2 = np.empty(1)
        DEC2 = np.empty(1)
        
    # Warn about non-native endian arrays
    if not all(is_native_endian(arr) for arr in [RA1, DEC1, weights1, RA2, DEC2, weights2]):
        warn('One or more input array has non-native endianness!  A copy will be made with the correct endianness.')
    RA1, DEC1, weights1, RA2, DEC2, weights2 = [convert_to_native_endian(arr) for arr in [RA1, DEC1, weights1, RA2, DEC2, weights2]]

    fix_ra_dec(RA1, DEC1)
    if autocorr == 0:
        fix_ra_dec(RA2, DEC2)

    if link_in_ra is True:
        link_in_dec = True

        
    # Passing None parameters breaks the parsing code, so avoid this
    kwargs = {}
    for k in ['weights1', 'weights2', 'weight_type', 'RA2', 'DEC2']:
        v = locals()[k]
        if v is not None:
            kwargs[k] = v

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    with sys_pipes():
      extn_results = DDtheta_mocks_extn(autocorr, nthreads, rbinfile,
                                        RA1, DEC1,
                                        verbose=verbose,
                                        link_in_dec=link_in_dec,
                                        link_in_ra=link_in_ra,
                                        output_thetaavg=output_thetaavg,
                                        fast_acos=fast_acos,
                                        ra_refine_factor=ra_refine_factor,
                                        dec_refine_factor=dec_refine_factor,
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
        os.remove(rbinfile)

    results_dtype = np.dtype([(bytes_to_native_str(b'thetamin'), np.float),
                              (bytes_to_native_str(b'thetamax'), np.float),
                              (bytes_to_native_str(b'thetaavg'), np.float),
                              (bytes_to_native_str(b'npairs'), np.uint64),
                              (bytes_to_native_str(b'weightavg'), np.float)])
    results = np.array(extn_results, dtype=results_dtype)

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
