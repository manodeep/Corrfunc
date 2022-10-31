#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper around the C extension for the counts-in-cells
for positions on the sky. Corresponding C codes are in ``mocks/vpf_mocks/``
while the python wrapper is in :py:mod:`Corrfunc.mocks.vpf_mocks`
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('vpf_mocks', )


def vpf_mocks(rmax, nbins, nspheres, numpN,
              threshold_ngb, centers_file, cosmology,
              RA, DEC, CZ,
              RAND_RA, RAND_DEC, RAND_CZ,
              verbose=False, is_comoving_dist=False,
              xbin_refine_factor=1, ybin_refine_factor=1,
              zbin_refine_factor=1, max_cells_per_dim=100,
              copy_particles=True, c_api_timer=False, isa=r'fastest'):
    """
    Function to compute the counts-in-cells on points on the sky. Suitable
    for mock catalogs and observed galaxies.

    Returns a numpy structured array containing the probability of a
    sphere of radius up to ``rmax`` containing ``0--numpN-1`` galaxies.

    Parameters
    ----------

    rmax : double
        Maximum radius of the sphere to place on the particles

    nbins : integer
        Number of bins in the counts-in-cells. Radius of first shell
        is rmax/nbins

    nspheres : integer (>= 0)
        Number of random spheres to place within the particle distribution.
        For a small number of spheres, the error is larger in the measured
        pN's.

    numpN : integer (>= 1)
        Governs how many unique pN's are to returned. If ``numpN`` is set to 1,
        then only the vpf (p0) is returned. For ``numpN=2``, p0 and p1 are
        returned.

        More explicitly, the columns in the results look like the following:

          ======   ==========================
          numpN    Columns in output
          ======   ==========================
             1      p0
             2      p0      p1
             3      p0      p1     p2
             4      p0      p1     p2     p3
          ======   ==========================

        and so on...

        Note: ``p0`` is the vpf

    threshold_ngb : integer
        Minimum number of random points needed in a ``rmax`` sphere such that
        it is considered to be entirely within the mock footprint. The
        command-line version, ``mocks/vpf/vpf_mocks.c``, assumes that the
        minimum number of randoms can be at most a 1-sigma deviation from
        the expected random number density.

    centers_file : string, filename
        A file containing random sphere centers. If the file does not exist,
        then a list of random centers will be written out. In that case, the
        randoms arrays, ``RAND_RA``, ``RAND_DEC`` and ``RAND_CZ`` are used to
        check that the sphere is entirely within the footprint. If the file
        does exist but either ``rmax`` is too small or there are not enough
        centers then the file will be overwritten.

        Note: If the centers file has to be written, the code will take
        significantly longer to finish. However, subsequent runs can re-use
        that centers file and will be faster.

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

    RA : array-like, real (float/double)
        The array of Right Ascensions for the first set of points. RA's
        are expected to be in [0.0, 360.0], but the code will try to fix cases
        where the RA's are in [-180, 180.0]. For peace of mind, always supply
        RA's in [0.0, 360.0].

        Calculations are done in the precision of the supplied arrays.

    DEC : array-like, real (float/double)
        Array of Declinations for the first set of points. DEC's are expected
        to be in the [-90.0, 90.0], but the code will try to fix cases where
        the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply
        DEC's in [-90.0, 90.0].

        Must be of same precision type as RA.

    CZ : array-like, real (float/double)
        Array of (Speed Of Light * Redshift) values for the first set of
        points. Code will try to detect cases where ``redshifts`` have been
        passed and multiply the entire array with the ``speed of light``.

        If ``is_comoving_dist`` is set, then ``CZ`` is interpreted as the
        co-moving distance, rather than (Speed Of Light * Redshift).

    RAND_RA : array-like, real (float/double)
        The array of Right Ascensions for the randoms. RA's are expected to be
        in [0.0, 360.0], but the code will try to fix cases where the RA's are
        in [-180, 180.0]. For peace of mind, always supply RA's in
        [0.0, 360.0].

        Must be of same precision type as RA/DEC/CZ.

    RAND_DEC : array-like, real (float/double)
        Array of Declinations for the randoms. DEC's are expected to be in the
        [-90.0, 90.0], but the code will try to fix cases where the DEC's are
        in [0.0, 180.0]. Again, for peace of mind, always supply DEC's in
        [-90.0, 90.0].

        Must be of same precision type as RA/DEC/CZ.

    RAND_CZ : array-like, real (float/double)
        Array of (Speed Of Light * Redshift) values for the randoms. Code
        will try to detect cases where ``redshifts`` have been
        passed and multiply the entire array with the ``speed of light``.

        If ``is_comoving_dist`` is set, then ``CZ2`` is interpreted as the
        co-moving distance, rather than ``(Speed Of Light * Redshift)``.

        Note: RAND_RA, RAND_DEC and RAND_CZ are only used when the
           ``centers_file``  needs to be written out. In that case, the
           RAND_RA, RAND_DEC, and RAND_CZ are used as random centers.

    verbose : boolean (default false)
        Boolean flag to control output of informational messages

    is_comoving_dist : boolean (default false)
        Boolean flag to indicate that ``cz`` values have already been
        converted into co-moving distances. This flag allows arbitrary
        cosmologies to be used in ``Corrfunc``.

    (xyz)bin_refine_factor : integer, default is (1, 1, 1); typically in [1-2]
        Controls the refinement on the cell sizes. Higher numbers might have
        a negative impact on runtime.

        Note: Since the counts in spheres calculation is symmetric
        in all 3 dimensions, the defaults are different from the clustering
        routines.

    max_cells_per_dim : integer, default is 100, typical values in [50-300]
        Controls the maximum number of cells per dimension. Total number of
        cells can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is
        too small relative to the boxsize (and increasing helps the runtime).

    copy_particles: boolean (default True)
        Boolean flag to make a copy of the particle positions
        If set to False, the particles will be re-ordered in-place

        .. versionadded:: 2.3.0

    c_api_timer : boolean (default false)
        Boolean flag to measure actual time spent in the C libraries. Here
        to allow for benchmarking and scaling studies.

    isa: string, case-insensitive (default ``fastest``)
        Controls the runtime dispatch for the instruction set to use. Options
        are: [``fastest``, ``avx512f``, ``avx``, ``sse42``, ``fallback``]

        Setting isa to ``fastest`` will pick the fastest available instruction
        set on the current computer. However, if you set ``isa`` to, say,
        ``avx`` and ``avx`` is not available on the computer, then the code
        will revert to using ``fallback`` (even though ``sse42`` might be
        available).

        Unless you are benchmarking the different instruction sets, you should
        always leave ``isa`` to the default value. And if you *are*
        benchmarking, then the string supplied here gets translated into an
        ``enum`` for the instruction set defined in ``utils/defs.h``.

    Returns
    --------

    results : Numpy structured array
        A numpy structured array containing [rmax, pN[numpN]] with ``nbins``
        elements. Each row contains the maximum radius of the sphere and the
        ``numpN`` elements in the ``pN`` array. Each element of this array
        contains the probability that a sphere of radius ``rmax`` contains
        *exactly* ``N`` galaxies. For example, pN[0] (p0, the void probibility
        function) is the probability that a sphere of radius ``rmax`` contains
        0 galaxies.

    api_time : float, optional
        Only returned if ``c_api_timer`` is set.  ``api_time`` measures only
        the time spent within the C library and ignores all python overhead.


    Example
    --------

    >>> from __future__ import print_function
    >>> import math
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import numpy as np
    >>> import Corrfunc
    >>> from Corrfunc.mocks.vpf_mocks import vpf_mocks
    >>> rmax = 10.0
    >>> nbins = 10
    >>> numbins_to_print = nbins
    >>> nspheres = 10000
    >>> numpN = 6
    >>> threshold_ngb = 1  # does not matter since we have the centers
    >>> cosmology = 1  # LasDamas cosmology
    >>> centers_file = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                      "../mocks/tests/data/",
    ...                      "Mr19_centers_xyz_forVPF_rmax_10Mpc.txt")
    >>> N = 1000000
    >>> boxsize = 420.0
    >>> seed = 42
    >>> np.random.seed(seed)
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
    >>> results = vpf_mocks(rmax, nbins, nspheres, numpN, threshold_ngb,
    ...                     centers_file, cosmology,
    ...                     RA, DEC, CZ,
    ...                     RA, DEC, CZ,
    ...                     is_comoving_dist=True)
    >>> for r in results:
    ...     print("{0:10.1f} ".format(r[0]), end="")
    ...     # doctest: +NORMALIZE_WHITESPACE
    ...     for pn in r[1]:
    ...         print("{0:10.3f} ".format(pn), end="")
    ...         # doctest: +NORMALIZE_WHITESPACE
    ...     print("") # doctest: +NORMALIZE_WHITESPACE
       1.0      0.999      0.001      0.000      0.000      0.000      0.000
       2.0      0.992      0.007      0.001      0.000      0.000      0.000
       3.0      0.982      0.009      0.005      0.002      0.001      0.000
       4.0      0.975      0.006      0.006      0.005      0.003      0.003
       5.0      0.971      0.004      0.003      0.003      0.004      0.003
       6.0      0.967      0.003      0.003      0.001      0.003      0.002
       7.0      0.962      0.004      0.002      0.003      0.002      0.001
       8.0      0.958      0.004      0.002      0.003      0.001      0.002
       9.0      0.953      0.003      0.003      0.002      0.003      0.001
      10.0      0.950      0.003      0.002      0.002      0.001      0.002

    """

    try:
        from Corrfunc._countpairs_mocks import countspheres_vpf_mocks\
            as vpf_extn
    except ImportError:
        msg = "Could not import the C extension for the Counts-in-Cells "\
              " (vpf)"
        raise ImportError(msg)

    import numpy as np
    from future.utils import bytes_to_native_str
    from Corrfunc.utils import translate_isa_string_to_enum,\
        convert_to_native_endian, sys_pipes

    # Ensure all input arrays are native endian
    RA, DEC, CZ, RAND_RA, RAND_DEC, RAND_CZ = [
            convert_to_native_endian(arr, warn=True) for arr in
            [RA, DEC, CZ, RAND_RA, RAND_DEC, RAND_CZ]]

    integer_isa = translate_isa_string_to_enum(isa)
    with sys_pipes():
      extn_results = vpf_extn(rmax, nbins, nspheres, numpN,
                              threshold_ngb, centers_file,
                              cosmology,
                              RA, DEC, CZ,
                              RAND_RA, RAND_DEC, RAND_CZ,
                              verbose=verbose,
                              is_comoving_dist=is_comoving_dist,
                              xbin_refine_factor=xbin_refine_factor,
                              ybin_refine_factor=ybin_refine_factor,
                              zbin_refine_factor=zbin_refine_factor,
                              max_cells_per_dim=max_cells_per_dim,
                              copy_particles=copy_particles,
                              c_api_timer=c_api_timer,
                              isa=integer_isa)
    if extn_results is None:
        msg = "RuntimeError occurred"
        raise RuntimeError(msg)
    else:
        extn_results, api_time = extn_results

    results_dtype = np.dtype([(bytes_to_native_str(b'rmax'), np.float64),
                              (bytes_to_native_str(b'pN'),
                               (np.float64, numpN))])
    nbin = len(extn_results)
    results = np.zeros(nbin, dtype=results_dtype)

    for ii, r in enumerate(extn_results):
        results['rmax'][ii] = r[0]
        if numpN == 1:
            results['pN'] = r[1]
        else:
            for j in range(numpN):
                results['pN'][ii][j] = r[1 + j]

    if not c_api_timer:
        return results
    else:
        return results, api_time


if __name__ == '__main__':
    import doctest
    doctest.testmod()
