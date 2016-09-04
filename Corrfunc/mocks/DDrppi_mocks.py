#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the pair counter in
``xi_mocks/DDrppi``. This python wrapper is in `~Corrfunc.mocks.DDrppi_mocks`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DDrppi_mocks', )


def DDrppi_mocks(autocorr, cosmology, nthreads, pimax, binfile,
                 RA1, DEC1, CZ1,
                 RA2=None, DEC2=None, CZ2=None,
                 is_comoving_dist=False,
                 verbose=False, output_rpavg=False,
                 fast_divide=False, c_api_timer=False, isa='fastest'):
    """
    Calculate the 2-D pair-counts corresponding to the projected correlation
    function, :math:`\\xi(r_p, \pi)`. Pairs which are separated by less
    than the ``rp`` bins (specified in ``binfile``) in the
    X-Y plane, and less than ``pimax`` in the Z-dimension are
    counted. The input positions are expected to be on-sky co-ordinates.
    This module is suitable for calculating correlation functions for mock
    catalogs.


    Returns a numpy structured array containing the pair counts for the
    specified bins.

    Note, that this module only returns pair counts and not the actual
    correlation function :math:`\\xi(r_p, \pi)`. See the
    ``xi_mocks/DDrppi/wprp_mocks.c`` for computing :math:`wp(r_p)` from
    the pair counts returned.
     
    Parameters:
    -----------
    autocorr: boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.

    cosmology: integer, required
        Integer choice for setting cosmology. Valid values are 1->LasDamas
        cosmology and 2->Planck cosmology. If you need arbitrary cosmology,
        easiest way is to convert the ``CZ`` values into co-moving distance,
        based on your preferred cosmology. Set ``is_comoving_dist=True``, to
        indicate that the co-moving distance conversion has already been done.
    
        Choices: 1 -> LasDamas cosmology. Om=0.25,  Ol=0.75
                 2 -> Planck   cosmology. Om=0.302, Ol=0.698
     
        To setup a new cosmology, add an entry to the function,
        ``init_cosmology`` in ``ROOT/utils/cosmology_params.c` and re-install
        the entire package.

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
       For string input: filename specifying the ``rp`` bins for
       ``DDrppi_mocks``. The file should contain white-space separated values
       of (rpmin, rpmax)  for each ``rp`` wanted. The bins do not need to be
       contiguous but must be in increasing order (smallest bins come first).

       For array-like input: A sequence of ``rp`` values that provides the
       bin-edges. For example, ``np.logspace(0.1, 10.0, 15)`` is a valid
       input, specifying 15 (logarithmic) bins between 0.1 and 10.0. This
       array does not need to be sorted.

    RA1: array-like, real (float/double)
        The array of Right Ascensions for the first set of points. RA's
        are expected to be in [0.0, 360.0], but the code will try to fix cases
        where the RA's are in [-180, 180.0]. For peace of mind, always supply
        RA's in [0.0, 360.0].

        Calculations are done in the precision of the supplied arrays.

    DEC1: array-like, real (float/double)
         Array of Declinations for the first set of points. DEC's are expected
         to be in the [-90.0, 90.0], but the code will try to fix cases where
         the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
         DEC's in [-90.0, 90.0].

         Must be of same precision type as RA1.

    CZ1: array-like, real (float/double)
         Array of (Speed Of Light * Redshift) values for the first set of
         points. Code will try to detect cases where ``redshifts`` have been
         passed and multiply the entire array with the ``speed of light``.

         * If is_comoving_dist is set, then ``CZ1`` is interpreted as the
         co-moving distance, rather than `cz`.

    RA2: array-like, real (float/double)
        The array of Right Ascensions for the second set of points. RA's
        are expected to be in [0.0, 360.0], but the code will try to fix cases
        where the RA's are in [-180, 180.0]. For peace of mind, always supply
        RA's in [0.0, 360.0].

        Must be of same precision type as RA1/DEC1/CZ1.

    DEC2: array-like, real (float/double)
         Array of Declinations for the second set of points. DEC's are expected
         to be in the [-90.0, 90.0], but the code will try to fix cases where
         the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
         DEC's in [-90.0, 90.0].

         Must be of same precision type as RA1/DEC1/CZ1.

    CZ2: array-like, real (float/double)
         Array of (Speed Of Light * Redshift) values for the second set of
         points. Code will try to detect cases where ``redshifts`` have been
         passed and multiply the entire array with the ``speed of light``.

         * If is_comoving_dist is set, then ``CZ2`` is interpreted as the
         co-moving distance, rather than `cz`.
        
         Must be of same precision type as RA1/DEC1/CZ1.

    is_comoving_dist: boolean (default false)
         Boolean flag to indicate that ``cz`` values have already been
         converted into co-moving distances. This flag allows arbitrary
         cosmologies to be used in ``Corrfunc``.

    verbose: boolean (default false)
       Boolean flag to control output of informational messages

    output_rpavg: boolean (default false)
       Boolean flag to output the average ``rp`` for each bin. Code will
       run slower if you set this flag. Also, note, if you are calculating
       in single-precision, ``rpavg`` will suffer from numerical loss of
       precision and can not be trusted. If you need accurate ``rpavg``
       values, then pass in double precision arrays for the particle positions.

    fast_divide: boolean (default false)
       Boolean flag to replace the division in ``AVX`` implementation with an
       approximate reciprocal, followed by a Newton-Raphson step. Improves
       runtime by ~15-20%. Loss of precision is at the 5-6th decimal place.

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
       for each radial bin specified in the ``binfile``. If ``output_ravg`` is
       not set, then ``ravg`` will be set to 0.0 for all bins. ``npairs``
       contains the number of pairs in that bin and can be used to compute the
       actual :math:`\\xi(r_p, \pi)` or :math:`\\wp(r_p)` by combining with
       (DR, RR) counts.

       if ``c_api_timer`` is set, then the return value is a tuple containing
       (results, api_time). ``api_time`` measures only the time spent within
       the C library and ignores all python overhead.

    Example:
    --------
    >>> import Corrfunc
    >>> from Corrfunc.mocks import DDrppi_mocks
    >>> import math
    >>> from os.path import dirname, abspath, join as pjoin
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                        "../xi_theory/tests/", "bins")
    
    >>> N = 100000
    >>> boxsize = 420.0
    >>> X = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    >>> Y = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    >>> Z = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    
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
                               pimax, binfile,
                               RA, DEC, CZ,
                               fast_divide=False,
                               verbose=False,
                               is_comoving_dist=True,
                               output_rpavg=True)

    >>> print("Results = {0}".format(results))

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
        return_file_with_rbins
    from future.utils import bytes_to_native_str
    
    if autocorr == 0:
        if RA2 is None or DEC2 is None or CZ2 is None:
            msg = "Must pass valid arrays for RA2/DEC2/CZ2 for "\
                  "computing cross-correlation"
            raise ValueError(msg)
    else:
        RA2 = np.empty(1)
        DEC2 = np.empty(1)
        CZ2 = np.empty(1)

    fix_ra_dec(RA1, DEC1)
    if autocorr == 0:
        fix_ra_dec(RA2, DEC2)

    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    extn_results, api_time = DDrppi_extn(autocorr, cosmology, nthreads,
                                         pimax, rbinfile,
                                         RA1, DEC1, CZ1,
                                         RA2, DEC2, CZ2,
                                         is_comoving_dist=is_comoving_dist,
                                         verbose=verbose,
                                         output_rpavg=output_rpavg,
                                         fast_divide=fast_divide,
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
    from os.path import dirname, abspath, join as pjoin

    print("\nRunning 2-D correlation function for mocks DD(rp,pi)")
    binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                    "../xi_theory/tests/", "bins")
    
    N = 100000
    boxsize = 420.0
    X = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    Y = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    Z = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)
    
    # Convert XYZ into RA/DEC/CZ
    CZ = np.sqrt(X*X + Y*Y + Z*Z)
    inv_cz = 1.0/CZ

    # Convert to unit sphere
    X *= inv_cz
    Y *= inv_cz
    Z *= inv_cz

    import math
    DEC = 90.0 - np.arccos(Z)*180.0/math.pi
    RA = (np.arctan2(Y, X)*180.0/math.pi) + 180.0

    autocorr = 1
    cosmology = 1
    nthreads = 2
    pimax = 40.0

    import time
    t0 = time.time()
    results, api_time = DDrppi_mocks(autocorr, cosmology, nthreads,
                                     pimax, binfile,
                                     RA, DEC, CZ,
                                     fast_divide=False,
                                     verbose=False,
                                     is_comoving_dist=True,
                                     c_api_timer=True,
                                     output_rpavg=True)
    t1 = time.time()
    print("Results from DDrppi_mocks: Time taken = {0:0.3f} sec "
          "Python time = {0:0.3f} sec".format(api_time, t1-t0))
    for r in results[0:10]:
        print("{0}".format(r))

