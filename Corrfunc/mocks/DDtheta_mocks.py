#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the angular correlation function
:math:`\\omega(\theta)`. Corresponding C routines are in mocks/wtheta/,
python interface is `~Corrfunc.mocks.DDtheta_mocks`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DDtheta_mocks',)


def DDtheta_mocks(autocorr, nthreads, binfile,
                  RA1, DEC1,
                  RA2=None, DEC2=None,
                  link_in_dec=True, link_in_ra=True,
                  verbose=False, output_thetaavg=False,
                  fast_acos=False, c_api_timer=False,
                  isa='fastest'):
    """
    Function to compute the angular correlation function for points on
    the sky (i.e., mock catalogs or observed galaxies).

    Returns a numpy structured array containing the pair counts for the
    specified angular bins.

    Note, that this module only returns pair counts and not the actual
    correlation function :math:`\\omega(\theta)`. See the
    ``mocks/wtheta/wtheta`` for computing :math:`\\omega(\theta)` from
    the pair counts returned.

    Parameters
    -----------

    autocorr: boolean, required
        Boolean flag for auto/cross-correlation. If autocorr is set to 1,
        then the second set of particle positions are not required.
    
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

    RA2: array-like, real (float/double)
       The array of Right Ascensions for the second set of points. RA's
       are expected to be in [0.0, 360.0], but the code will try to fix cases
       where the RA's are in [-180, 180.0]. For peace of mind, always supply
       RA's in [0.0, 360.0].

        Must be of same precision type as RA1/DEC1.

    DEC2: array-like, real (float/double)
       Array of Declinations for the second set of points. DEC's are expected
       to be in the [-90.0, 90.0], but the code will try to fix cases where
       the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
       DEC's in [-90.0, 90.0].

       Must be of same precision type as RA1/DEC1.

    link_in_dec: boolean (default True)
       Boolean flag to create lattice in Declination. Code runs faster with
       this option. However, if the angular separations are too small, then
       linking in declination might produce incorrect results. When running
       for the first time, check your results by comparing with the output
       of the code for ``link_in_dec=False`` and ``link_in_ra=False``.

    link_in_ra: boolean (default True)
       Boolean flag to create lattice in Right Ascension. Setting this option
       implies ``link_in_dec=True``. Similar considerations as ``link_in_dec``
       described above.

       *Note*: If you disable both ``link_in_dec`` and ``link_in_ra``, then
       the code reduces to a brute-force pair counter. No lattices are created
       at all. For very small angular separations, that might be the most
       numerically stable method.

    verbose: boolean (default false)
       Boolean flag to control output of informational messages
    
    output_thetaavg: boolean (default false)
       Boolean flag to output the average ``\theta`` for each bin. Code will
       run slower if you set this flag. Also, note, if you are calculating
       in single-precision, ``thetaavg`` will suffer from numerical loss of
       precision and can not be trusted. If you need accurate ``ravg``
       values, then pass in double precision arrays for ``XYZ``.

       **NOTE** Code will run significantly slower if you enable this option.
       Use ``fast_acos`` if you can tolerate some loss of precision.

    fast_acos: boolean (default false)
       Flag to use numerical approximation for the ``arccos`` - gives better
       performance at the expense of some precision. Relevant only if
       ``output_thetaavg==True``.

       Developers: Two versions already coded up in ``utils/fast_acos.h``,
       so you can choose the version you want. There are also notes on how
       to implement faster (and less accurate) functions, particularly relevant
       if you know your ``theta`` range is limited. If you implement a new
       version, then you will have to reinstall the entire Corrfunc package.
       
       Note that tests will fail if you run the tests with``fast_acos=True``.
       
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
       projected correlation function while ``npairs`` contains the number of
       unique pairs in that bin.

    Example
    --------
    
    >>> import numpy as np
    >>> import time
    >>> from math import pi
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                  "../mocks/tests/", "angular_bins")
    >>> N = 100000
    >>> nthreads = 4
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> RA = np.random.uniform(0.0, 2.0*pi, N)*180.0/pi
    >>> cos_theta = np.random.uniform(-1.0, 1.0, N)
    >>> DEC = 90.0 - np.arccos(cos_theta)*180.0/pi
    >>> autocorr = 1
    >>> results = DDtheta_mocks(autocorr, nthreads, binfile,
                                RA, DEC,
                                verbose=True)

    """

    try:
        from Corrfunc._countpairs_mocks import countpairs_theta_mocks as\
            DDtheta_mocks_extn
    except ImportError:
        msg = "Could not import the C extension for the angular "\
              "correlation function for mocks."
        raise ImportError(msg)

    import numpy as np
    from Corrfunc.utils import translate_isa_string_to_enum, fix_ra_dec,\
        return_file_with_rbins
    from future.utils import bytes_to_native_str
    
    if autocorr == 0:
        if RA2 is None or DEC2 is None:
            msg = "Must pass valid arrays for RA2/DEC2 for "\
                  "computing cross-correlation"
            raise ValueError(msg)
    else:
        RA2 = np.empty(1)
        DEC2 = np.empty(1)

    fix_ra_dec(RA1, DEC1)
    if autocorr == 0:
        fix_ra_dec(RA2, DEC2)

    if link_in_ra is True:
        link_in_dec = True
        
    integer_isa = translate_isa_string_to_enum(isa)
    rbinfile, delete_after_use = return_file_with_rbins(binfile)
    extn_results, api_time = DDtheta_mocks_extn(autocorr, nthreads, rbinfile,
                                                RA1, DEC1,
                                                RA2, DEC2,
                                                verbose=verbose,
                                                link_in_dec=link_in_dec,
                                                link_in_ra=link_in_ra,
                                                output_thetaavg=output_thetaavg,
                                                c_api_timer=c_api_timer,
                                                fast_acos=fast_acos,
                                                isa=integer_isa)
        
    if extn_results is None:
        msg = "RuntimeError occurred"
        raise RuntimeError(msg)

    if delete_after_use:
        import os
        os.remove(rbinfile)
    
    results_dtype = np.dtype([(bytes_to_native_str(b'thetamin'), np.float),
                              (bytes_to_native_str(b'thetamax'), np.float),
                              (bytes_to_native_str(b'thetaavg'), np.float),
                              (bytes_to_native_str(b'npairs'), np.uint64)])
    
    nbin = len(extn_results)
    results = np.zeros(nbin, dtype=results_dtype)
    
    for ii, r in enumerate(extn_results):
        results['thetamin'][ii] = r[0]
        results['thetamax'][ii] = r[1]
        results['thetaavg'][ii] = r[2]
        results['npairs'][ii] = r[3]

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import numpy as np
    import time
    from math import pi
    from os.path import dirname, abspath, join as pjoin
    import Corrfunc
    print("\nRunning angular correlation function w(theta)")

    binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                    "../mocks/tests/", "angular_bins")
    
    N = 100000
    nthreads = 4
    t0 = time.time()
    seed = 42
    np.random.seed(seed)
    
    # Faster way of generating random RA's and DEC's
    RA = np.random.uniform(0.0, 2.0*pi, N)*180.0/pi
    cos_theta = np.random.uniform(-1.0, 1.0, N)
    DEC = 90.0 - np.arccos(cos_theta)*180.0/pi

    autocorr = 1
    results, api_time = DDtheta_mocks(autocorr, nthreads, binfile,
                                      RA, DEC,
                                      output_thetaavg=False,
                                      c_api_timer=True,
                                      verbose=True)

    t1 = time.time()
    print("Results from DDtheta_mocks (Npts = {0}): API time = {1:0.3f} sec "
          "python overhead = {2:0.3f} sec".format(N, api_time,
                                                  (t1-t0) - api_time))
    for r in results:
        print("{0}".format(r))

