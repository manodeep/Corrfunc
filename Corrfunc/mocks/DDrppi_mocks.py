#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the pair counter in
``mocks/DDrppi_mocks/``. This python wrapper is
:py:mod:`Corrfunc.mocks.DDrppi_mocks`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DDrppi_mocks', )


def DDrppi_mocks(autocorr, cosmology, nthreads, pimax, binfile,
                 RA1, DEC1, CZ1, weights1=None,
                 RA2=None, DEC2=None, CZ2=None, weights2=None,
                 is_comoving_dist=False,
                 verbose=False, output_rpavg=False,
                 fast_divide_and_NR_steps=0,
                 xbin_refine_factor=2, ybin_refine_factor=2,
                 zbin_refine_factor=1, max_cells_per_dim=100,
                 copy_particles=True, enable_min_sep_opt=True,
                 c_api_timer=False, isa=r'fastest', weight_type=None):
    """
    Calculate the pair-counts corresponding to the 2-D correlation
    function, :math:`\\xi(r_p, \\pi)`. Pairs which are separated by less
    than the ``rp`` bins (specified in ``binfile``) in the
    X-Y plane, and less than ``pimax`` in the Z-dimension are
    counted. The input positions are expected to be on-sky co-ordinates.
    This module is suitable for calculating correlation functions for mock
    catalogs.

    If ``weights`` are provided, the resulting pair counts are weighted.  The
    weighting scheme depends on ``weight_type``.

    Returns a numpy structured array containing the pair counts for the
    specified bins.


    .. note:: that this module only returns pair counts and not the actual
       correlation function :math:`\\xi(r_p, \\pi)` or :math:`wp(r_p)`. See the
       utilities :py:mod:`Corrfunc.utils.convert_3d_counts_to_cf` and
       :py:mod:`Corrfunc.utils.convert_rp_pi_counts_to_wp` for computing
       :math:`\\xi(r_p, \\pi)` and :math:`wp(r_p)` respectively from the
       pair counts.


    Parameters
    -----------

    autocorr : boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.

    cosmology : integer, required
        Integer choice for setting cosmology. Valid values are 1->LasDamas
        cosmology and 2->Planck cosmology. If you need arbitrary cosmology,
        easiest way is to convert the ``CZ`` values into co-moving distance,
        based on your preferred cosmology. Set ``is_comoving_dist=True``, to
        indicate that the co-moving distance conversion has already been done.

        Choices:
                 1. LasDamas cosmology. :math:`\\Omega_m=0.25`, :math:`\\Omega_\\Lambda=0.75`
                 2. Planck   cosmology. :math:`\\Omega_m=0.302`, :math:`\\Omega_\\Lambda=0.698`

        To setup a new cosmology, add an entry to the function,
        ``init_cosmology`` in ``ROOT/utils/cosmology_params.c`` and re-install
        the entire package.

    nthreads : integer
        The number of OpenMP threads to use. Has no effect if OpenMP was not
        enabled during library compilation.

    pimax : double
        A double-precision value for the maximum separation along
        the Z-dimension.

        Distances along the :math:`\\pi` direction are binned with unit
        depth. For instance, if ``pimax=40``, then 40 bins will be created
        along the ``pi`` direction. Only pairs with ``0 <= dz < pimax``
        are counted (no equality).

    binfile: string or an list/array of floats
        For string input: filename specifying the ``rp`` bins for
        ``DDrppi_mocks``. The file should contain white-space separated values
        of (rpmin, rpmax)  for each ``rp`` wanted. The bins need to be
        contiguous and sorted in increasing order (smallest bins come first).

        For array-like input: A sequence of ``rp`` values that provides the
        bin-edges. For example,
        ``np.logspace(np.log10(0.1), np.log10(10.0), 15)`` is a valid
        input specifying **14** (logarithmic) bins between 0.1 and 10.0. This
        array does not need to be sorted.

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

    CZ1 : array-like, real (float/double)
        Array of (Speed Of Light * Redshift) values for the first set of
        points. Code will try to detect cases where ``redshifts`` have been
        passed and multiply the entire array with the ``speed of light``.

        If is_comoving_dist is set, then ``CZ1`` is interpreted as the
        co-moving distance, rather than `cz`.

    weights1 : array_like, real (float/double), optional
        A scalar, or an array of weights of shape (n_weights, n_positions) or (n_positions,).
        `weight_type` specifies how these weights are used; results are returned
        in the `weightavg` field.  If only one of weights1 and weights2 is
        specified, the other will be set to uniform weights.

    RA2 : array-like, real (float/double)
        The array of Right Ascensions for the second set of points. RA's
        are expected to be in [0.0, 360.0], but the code will try to fix cases
        where the RA's are in [-180, 180.0]. For peace of mind, always supply
        RA's in [0.0, 360.0].

        Must be of same precision type as RA1/DEC1/CZ1.

    DEC2 : array-like, real (float/double)
        Array of Declinations for the second set of points. DEC's are expected
        to be in the [-90.0, 90.0], but the code will try to fix cases where
        the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
        DEC's in [-90.0, 90.0].

        Must be of same precision type as RA1/DEC1/CZ1.

    CZ2 : array-like, real (float/double)
        Array of (Speed Of Light * Redshift) values for the second set of
        points. Code will try to detect cases where ``redshifts`` have been
        passed and multiply the entire array with the ``speed of light``.

        If is_comoving_dist is set, then ``CZ2`` is interpreted as the
        co-moving distance, rather than `cz`.

        Must be of same precision type as RA1/DEC1/CZ1.

    weights2 : array-like, real (float/double), optional
        Same as weights1, but for the second set of positions

    is_comoving_dist : boolean (default false)
        Boolean flag to indicate that ``cz`` values have already been
        converted into co-moving distances. This flag allows arbitrary
        cosmologies to be used in ``Corrfunc``.

    verbose : boolean (default false)
        Boolean flag to control output of informational messages

    output_rpavg : boolean (default false)
        Boolean flag to output the average ``rp`` for each bin. Code will
        run slower if you set this flag.

        If you are calculating in single-precision, ``rpavg`` will suffer
        suffer from numerical loss of precision and can not be trusted. If
        you need accurate ``rpavg`` values, then pass in double precision
        arrays for the particle positions.

    fast_divide_and_NR_steps: integer (default 0)
        Replaces the division in ``AVX`` implementation with an approximate
        reciprocal, followed by ``fast_divide_and_NR_steps`` of Newton-Raphson.
        Can improve runtime by ~15-20% on older computers. Value of 0 uses
        the standard division operation.

    (xyz)bin_refine_factor : integer, default is (2,2,1); typically within [1-3]
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

    c_api_timer : boolean (default false)
        Boolean flag to measure actual time spent in the C libraries. Here
        to allow for benchmarking and scaling studies.

    isa: string, case-insensitive (default ``fastest``)
       Controls the runtime dispatch for the instruction set to use. Possible
       options are: [``fastest``, ``avx512f``, ``avx``, ``sse42``, ``fallback``]

       Setting isa to ``fastest`` will pick the fastest available instruction
       set on the current computer. However, if you set ``isa`` to, say,
       ``avx`` and ``avx`` is not available on the computer, then the code will
       revert to using ``fallback`` (even though ``sse42`` might be available).

       Unless you are benchmarking the different instruction sets, you should
       always leave ``isa`` to the default value. And if you *are*
       benchmarking, then the string supplied here gets translated into an
       ``enum`` for the instruction set defined in ``utils/defs.h``.

    weight_type : string, optional (default None)
        The type of weighting to apply.  One of ["pair_product", None].

    Returns
    --------

    results : Numpy structured array

        A numpy structured array containing [rpmin, rpmax, rpavg, pimax,
        npairs, weightavg] for each radial bin specified in the ``binfile``.
        If ``output_ravg`` is not set, then ``rpavg`` will be set to 0.0 for
        all bins; similarly for ``weightavg``. ``npairs`` contains the number
        of pairs in that bin and can be used to compute the actual
        :math:`\\xi(r_p, \\pi)` or :math:`wp(rp)` by combining with
        (DR, RR) counts.

    api_time : float, optional
        Only returned if ``c_api_timer`` is set.  ``api_time`` measures only
        the time spent within the C library and ignores all python overhead.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
    >>> import math
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../mocks/tests/", "bins")
    >>> N = 100000
    >>> boxsize = 420.0
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> X = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    >>> Y = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    >>> Z = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    >>> weights = np.ones_like(X)
    >>> CZ = np.sqrt(X*X + Y*Y + Z*Z)
    >>> inv_cz = 1.0/CZ
    >>> X *= inv_cz
    >>> Y *= inv_cz
    >>> Z *= inv_cz
    >>> DEC = 90.0 - np.arccos(Z)*180.0/math.pi
    >>> RA = (np.arctan2(Y, X)*180.0/math.pi) + 180.0
    >>> autocorr = 1
    >>> cosmology = 1
    >>> nthreads = 2
    >>> pimax = 40.0
    >>> results = DDrppi_mocks(autocorr, cosmology, nthreads,
    ...                        pimax, binfile, RA, DEC, CZ,
    ...                        weights1=weights, weight_type='pair_product',
    ...                        output_rpavg=True, is_comoving_dist=True)
    >>> for r in results[519:]: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.1f}"
    ...                               " {4:10d} {5:10.6f}".format(r['rmin'], r['rmax'],
    ...                               r['rpavg'], r['pimax'], r['npairs'], r['weightavg']))
    ...                         # doctest: +NORMALIZE_WHITESPACE
     11.359969  16.852277  14.285169       40.0     104850   1.000000
     16.852277  25.000000  21.181246        1.0     274144   1.000000
     16.852277  25.000000  21.190844        2.0     272876   1.000000
     16.852277  25.000000  21.183321        3.0     272294   1.000000
     16.852277  25.000000  21.188486        4.0     272506   1.000000
     16.852277  25.000000  21.170832        5.0     272100   1.000000
     16.852277  25.000000  21.165379        6.0     271788   1.000000
     16.852277  25.000000  21.175246        7.0     270040   1.000000
     16.852277  25.000000  21.187417        8.0     269492   1.000000
     16.852277  25.000000  21.172066        9.0     269682   1.000000
     16.852277  25.000000  21.182460       10.0     268266   1.000000
     16.852277  25.000000  21.170594       11.0     268744   1.000000
     16.852277  25.000000  21.178608       12.0     266820   1.000000
     16.852277  25.000000  21.187184       13.0     266510   1.000000
     16.852277  25.000000  21.184937       14.0     265484   1.000000
     16.852277  25.000000  21.180184       15.0     265258   1.000000
     16.852277  25.000000  21.191504       16.0     262952   1.000000
     16.852277  25.000000  21.187746       17.0     262602   1.000000
     16.852277  25.000000  21.189778       18.0     260206   1.000000
     16.852277  25.000000  21.188882       19.0     259410   1.000000
     16.852277  25.000000  21.185684       20.0     256806   1.000000
     16.852277  25.000000  21.194036       21.0     255574   1.000000
     16.852277  25.000000  21.184115       22.0     255406   1.000000
     16.852277  25.000000  21.178255       23.0     252394   1.000000
     16.852277  25.000000  21.184644       24.0     252220   1.000000
     16.852277  25.000000  21.187020       25.0     251668   1.000000
     16.852277  25.000000  21.183827       26.0     249648   1.000000
     16.852277  25.000000  21.183121       27.0     247160   1.000000
     16.852277  25.000000  21.180872       28.0     246238   1.000000
     16.852277  25.000000  21.185251       29.0     246030   1.000000
     16.852277  25.000000  21.183488       30.0     242124   1.000000
     16.852277  25.000000  21.194538       31.0     242426   1.000000
     16.852277  25.000000  21.190702       32.0     239778   1.000000
     16.852277  25.000000  21.188985       33.0     239046   1.000000
     16.852277  25.000000  21.187092       34.0     237640   1.000000
     16.852277  25.000000  21.185515       35.0     236256   1.000000
     16.852277  25.000000  21.190278       36.0     233536   1.000000
     16.852277  25.000000  21.183240       37.0     233274   1.000000
     16.852277  25.000000  21.183796       38.0     231628   1.000000
     16.852277  25.000000  21.200668       39.0     230378   1.000000
     16.852277  25.000000  21.181153       40.0     229006   1.000000

    """
    try:
        from Corrfunc._countpairs_mocks import countpairs_rp_pi_mocks as\
            DDrppi_extn
    except ImportError:
        msg = "Could not import the C extension for the on-sky"\
              "pair counter."
        raise ImportError(msg)

    import numpy as np
    from Corrfunc.utils import translate_isa_string_to_enum, fix_ra_dec,\
        return_file_with_rbins, convert_to_native_endian,\
        sys_pipes, process_weights
    from future.utils import bytes_to_native_str

    if not autocorr:
        if RA2 is None or DEC2 is None or CZ2 is None:
            msg = "Must pass valid arrays for RA2/DEC2/CZ2 for "\
                  "computing cross-correlation"
            raise ValueError(msg)
    else:
        RA2 = np.empty(1)
        DEC2 = np.empty(1)
        CZ2 = np.empty(1)

    weights1, weights2 = process_weights(weights1, weights2, RA1, RA2, weight_type, autocorr)

    # Ensure all input arrays are native endian
    RA1, DEC1, CZ1, weights1, RA2, DEC2, CZ2, weights2 = [
            convert_to_native_endian(arr, warn=True) for arr in
            [RA1, DEC1, CZ1, weights1, RA2, DEC2, CZ2, weights2]]

    fix_ra_dec(RA1, DEC1)
    if autocorr == 0:
        fix_ra_dec(RA2, DEC2)

    # Passing None parameters breaks the parsing code, so avoid this
    kwargs = {}
    for k in ['weights1', 'weights2', 'weight_type', 'RA2', 'DEC2', 'CZ2']:
        v = locals()[k]
        if v is not None:
            kwargs[k] = v

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    with sys_pipes():
        extn_results = DDrppi_extn(autocorr, cosmology, nthreads,
                                   pimax, rbinfile,
                                   RA1, DEC1, CZ1,
                                   is_comoving_dist=is_comoving_dist,
                                   verbose=verbose,
                                   output_rpavg=output_rpavg,
                                   fast_divide_and_NR_steps=fast_divide_and_NR_steps,
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
                              (bytes_to_native_str(b'rpavg'), np.float64),
                              (bytes_to_native_str(b'pimax'), np.float64),
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
