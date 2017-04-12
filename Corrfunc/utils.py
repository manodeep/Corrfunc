#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A set of utility routines
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
from os.path import exists as file_exists

__all__ = ['convert_3d_counts_to_cf', 'convert_rp_pi_counts_to_wp',
           'translate_isa_string_to_enum', 'return_file_with_rbins',
           'fix_ra_dec', 'fix_cz', 'compute_nbins', 'gridlink_sphere', ]
if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]

try:
    xrange
except NameError:
    xrange = range

try:
    long
except NameError:
    long = int
    
def convert_3d_counts_to_cf(ND1, ND2, NR1, NR2,
                            D1D2, D1R2, D2R1, R1R2,
                            estimator='LS'):
    """
    Converts raw pair counts to a correlation function.

    Parameters
    ----------

    ND1 : integer
       Number of points in the first dataset

    ND2 : integer
        Number of points in the second dataset

    NR1 : integer
        Number of points in the randoms for first dataset

    NR2 : integer
        Number of points in the randoms for second dataset

    D1D2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and D2

    D1R2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and R2

    D2R1 : array-like, integer
        Pair-counts for the cross-correlation between D2 and R1

    R1R2 : array-like, integer
        Pair-counts for the cross-correlation between R1 and R2

    For all of these pair-counts arrays, the corresponding ``numpy``
    struct returned by the theory/mocks modules can also be passed

    estimator: string, default='LS' (Landy-Szalay)
        The kind of estimator to use for computing the correlation
        function. Currently, only supports Landy-Szalay

    Returns
    ---------

    cf : A numpy array
        The correlation function, calculated using the chosen estimator,
        is returned. NAN is returned for the bins where the ``RR`` count
        is 0.


    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from Corrfunc.theory.DD import DD
    >>> from Corrfunc.io import read_catalog
    >>> from Corrfunc.utils import convert_3d_counts_to_cf
    >>> X, Y, Z = read_catalog()
    >>> N = len(X)
    >>> boxsize = 420.0
    >>> rand_N = 3*N
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> rand_X = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Y = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Z = np.random.uniform(0, boxsize, rand_N)
    >>> nthreads = 2
    >>> rmin = 0.1
    >>> rmax = 15.0
    >>> nbins = 10
    >>> bins = np.linspace(rmin, rmax, nbins + 1)
    >>> autocorr = 1
    >>> DD_counts = DD(autocorr, nthreads, bins, X, Y, Z)
    >>> autocorr = 0
    >>> DR_counts = DD(autocorr, nthreads, bins,
    ...                X, Y, Z,
    ...                X2=rand_X, Y2=rand_Y, Z2=rand_Z)
    >>> autocorr = 1
    >>> RR_counts = DD(autocorr, nthreads, bins, rand_X, rand_Y, rand_Z)
    >>> cf = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
    ...                              DD_counts, DR_counts,
    ...                              DR_counts, RR_counts)
    >>> for xi in cf: print("{0:10.6f}".format(xi))
    ...                    # doctest: +NORMALIZE_WHITESPACE
    22.769019
     3.612709
     1.621372
     1.000969
     0.691646
     0.511819
     0.398872
     0.318815
     0.255643
     0.207759

    """

    import numpy as np
    pair_counts = dict()
    fields = ['D1D2', 'D1R2', 'D2R1', 'R1R2']
    arrays = [D1D2, D1R2, D2R1, R1R2]
    for (field, array) in zip(fields, arrays):
        try:
            npairs = array['npairs']
            pair_counts[field] = npairs
        except IndexError:
            pair_counts[field] = array

    nbins = len(pair_counts['D1D2'])
    if (nbins != len(pair_counts['D1R2'])) or \
       (nbins != len(pair_counts['D2R1'])) or \
       (nbins != len(pair_counts['R1R2'])):
        msg = 'Pair counts must have the same number of elements (same bins)'
        raise ValueError(msg)

    nonzero = pair_counts['R1R2'] > 0
    if 'LS' in estimator or 'Landy' in estimator:
        fN1 = np.float(NR1) / np.float(ND1)
        fN2 = np.float(NR2) / np.float(ND2)
        cf = np.zeros(nbins)
        cf[:] = np.nan
        cf[nonzero] = (fN1 * fN2 * pair_counts['D1D2'][nonzero] -
                       fN1 * pair_counts['D1R2'][nonzero] -
                       fN2 * pair_counts['D2R1'][nonzero] +
                       pair_counts['R1R2'][nonzero]) / pair_counts['R1R2'][nonzero]
        if len(cf) != nbins:
            msg = 'Bug in code. Calculated correlation function does not '\
                  'have the same number of bins as input arrays. Input bins '\
                  '={0} bins in (wrong) calculated correlation = {1}'.format(
                      nbins, len(cf))
            raise RuntimeError(msg)
    else:
        msg = "Only the Landy-Szalay estimator is supported. Pass estimator"\
              "='LS'. (Got estimator = {0})".format(estimator)
        raise ValueError(msg)

    return cf


def convert_rp_pi_counts_to_wp(ND1, ND2, NR1, NR2,
                               D1D2, D1R2, D2R1, R1R2,
                               nrpbins, pimax, dpi=1.0,
                               estimator='LS'):
    """
    Converts raw pair counts to a correlation function.

    Parameters
    ----------

    ND1 : integer
       Number of points in the first dataset

    ND2 : integer
        Number of points in the second dataset

    NR1 : integer
        Number of points in the randoms for first dataset

    NR2 : integer
        Number of points in the randoms for second dataset

    D1D2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and D2

    D1R2 : array-like, integer
        Pair-counts for the cross-correlation between D1 and R2

    D2R1 : array-like, integer
        Pair-counts for the cross-correlation between D2 and R1

    R1R2 : array-like, integer
        Pair-counts for the cross-correlation between R1 and R2

    For all of these pair-counts arrays, the corresponding ``numpy``
    struct returned by the theory/mocks modules can also be passed

    nrpbins : integer
        Number of bins in ``rp``

    pimax : float
        Integration distance along the line of sight direction

    dpi : float, default=1.0 Mpc/h
        Binsize in the line of sight direction

    estimator: string, default='LS' (Landy-Szalay)
        The kind of estimator to use for computing the correlation
        function. Currently, only supports Landy-Szalay

    Returns
    ---------

    wp : A numpy array
        The projected correlation function, calculated using the chosen
        estimator, is returned. If *any* of the ``pi`` bins (in an ``rp``
        bin) contains 0 for the ``RR`` counts, then ``NAN`` is returned
        for that ``rp`` bin.

    Example
    --------

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from Corrfunc.theory.DDrppi import DDrppi
    >>> from Corrfunc.io import read_catalog
    >>> from Corrfunc.utils import convert_rp_pi_counts_to_wp
    >>> X, Y, Z = read_catalog()
    >>> N = len(X)
    >>> boxsize = 420.0
    >>> rand_N = 3*N
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> rand_X = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Y = np.random.uniform(0, boxsize, rand_N)
    >>> rand_Z = np.random.uniform(0, boxsize, rand_N)
    >>> nthreads = 4
    >>> pimax = 40.0
    >>> nrpbins = 20
    >>> rpmin = 0.1
    >>> rpmax = 10.0
    >>> bins = np.linspace(rpmin, rpmax, nrpbins + 1)
    >>> autocorr = 1
    >>> DD_counts = DDrppi(autocorr, nthreads, pimax, bins,
    ...                    X, Y, Z)
    >>> autocorr = 0
    >>> DR_counts = DDrppi(autocorr, nthreads, pimax, bins,
    ...                    X, Y, Z,
    ...                    X2=rand_X, Y2=rand_Y, Z2=rand_Z)
    >>> autocorr = 1
    >>> RR_counts = DDrppi(autocorr, nthreads, pimax, bins,
    ...                    rand_X, rand_Y, rand_Z)
    >>> wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
    ...                                 DD_counts, DR_counts,
    ...                                 DR_counts, RR_counts,
    ...                                 nrpbins, pimax)
    >>> for w in wp: print("{0:10.6f}".format(w))
    ...                    # doctest: +NORMALIZE_WHITESPACE
    187.592199
     83.059181
     53.200599
     40.389354
     33.356371
     29.045476
     26.088133
     23.628340
     21.703961
     20.153125
     18.724781
     17.433235
     16.287183
     15.443230
     14.436193
     13.592727
     12.921226
     12.330074
     11.696364
     11.208365

    """
    
    import numpy as np
    if dpi <= 0.0:
        msg = 'Binsize along the line of sight (dpi) = {0}'\
              'must be positive'.format(dpi)
        raise ValueError(msg)

    xirppi = convert_3d_counts_to_cf(ND1, ND2, NR1, NR2,
                                     D1D2, D1R2, D2R1, R1R2,
                                     estimator=estimator)
    wp = np.empty(nrpbins)
    npibins = len(xirppi) // nrpbins
    if ((npibins * nrpbins) != len(xirppi)):
        msg = 'Number of pi bins could not be calculated correctly.'\
              'Expected to find that the total number of bins = {0} '\
              'would be the product of the number of pi bins = {1} '\
              'and the number of rp bins = {2}'.format(len(xirppi),
                                                       npibins,
                                                       nrpbins)
        raise ValueError(msg)

    # Check that dpi/pimax/npibins are consistent
    # Preventing issue #96 (https://github.com/manodeep/Corrfunc/issues/96)
    # where npibins would be calculated incorrectly, and the summation would
    # be wrong.
    if (dpi*npibins != pimax):
        msg = 'Pimax = {0} should be equal to the product of '\
              'npibins = {1} and dpi = {2}. Check your binning scheme.'\
              .format(pimax, npibins, dpi)
        raise ValueError(msg)
    
    for i in range(nrpbins):
        wp[i] = 2.0 * dpi * np.sum(xirppi[i * npibins:(i + 1) * npibins])

    return wp


def return_file_with_rbins(rbins):
    """
    Helper function to ensure that the ``binfile`` required by the Corrfunc
    extensions is a actually a string.

    Checks if the input is a string and file; return if True. If not, and
    the input is an array, then a temporary file is created and the contents
    of rbins is written out.

    Parameters
    -----------
    rbins: string or array-like
       Expected to be a string or an array containing the bins

    Returns
    ---------
    binfile: string, filename
       If the input ``rbins`` was a valid filename, then returns the same
       string. If ``rbins`` was an array, then this function creates a
       temporary file with the contents of the ``rbins`` arrays. This
       temporary filename is returned

    """

    is_string = False
    delete_after_use = False
    try:
        if isinstance(rbins, basestring):
            is_string = True
    except NameError:
        if isinstance(rbins, str):
            is_string = True

    if is_string:
        if file_exists(rbins):
            delete_after_use = False
            return rbins, delete_after_use
        else:
            msg = "Could not find file = `{0}` containing the bins"\
                  .format(rbins)
            raise IOError(msg)

    # For a valid bin specifier, there must be at least 1 bin.
    if len(rbins) >= 1:
        import tempfile
        rbins = sorted(rbins)
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as f:
            for i in range(len(rbins) - 1):
                f.write("{0} {1}\n".format(rbins[i], rbins[i + 1]))

            tmpfilename = f.name

        delete_after_use = True
        return tmpfilename, delete_after_use

    msg = "Input `binfile` was not a valid array (>= 1 element)."\
          "Num elements = {0}".format(len(rbins))
    raise TypeError(msg)


def fix_cz(cz):
    """
    Multiplies the input array by speed of light, if the input values are
    too small.

    Essentially, converts redshift into `cz`, if the user passed
    redshifts instead of `cz`.

    Parameters
    -----------
    cz: array-like, reals
       An array containing ``[Speed of Light *] redshift`` values.

    Returns
    ---------
    cz: array-like
       Actual ``cz`` values, multiplying the input ``cz`` array by the
       ``Speed of Light``, if ``redshift`` values were passed as input ``cz``.

    """

    # if I find that max cz is smaller than this threshold,
    # then I will assume z has been supplied rather than cz
    max_cz_threshold = 10.0
    try:
        input_dtype = cz.dtype
    except:
        msg = "Input cz array must be a numpy array"
        raise TypeError(msg)

    if max(cz) < max_cz_threshold:
        speed_of_light = 299800.0
        cz *= speed_of_light

    return cz.astype(input_dtype)


def fix_ra_dec(ra, dec):
    """
    Wraps input RA and DEC values into range expected by the extensions.

    Parameters
    ------------
    RA: array-like, units must be degrees
       Right Ascension values (astronomical longitude)

    DEC: array-like, units must be degrees
       Declination values (astronomical latitude)

    Returns
    --------
    Tuple (RA, DEC): array-like
         RA is wrapped into range [0.0, 360.0]
         Declination is wrapped into range [-90.0, 90.0]

    """

    try:
        input_dtype = ra.dtype
    except:
        msg = "Input RA array must be a numpy array"
        raise TypeError(msg)

    if ra is None or dec is None:
        msg = "RA or DEC must be valid arrays"
        raise ValueError(msg)

    if min(ra) < 0.0:
        print("Warning: found negative RA values, wrapping into [0.0, 360.0] "
              " range")
        ra += 180.0

    if max(dec) > 90.0:
        print("Warning: found DEC values more than 90.0; wrapping into "
              "[-90.0, 90.0] range")
        dec += 90.0

    return ra.astype(input_dtype), dec.astype(input_dtype)


def translate_isa_string_to_enum(isa):
    """
    Helper function to convert an user-supplied string to the
    underlying enum in the C-API. The extensions only have specific
    implementations for AVX, SSE42 and FALLBACK. Any other value
    will raise a ValueError.

    Parameters
    ------------
    isa: string
       A string containing the desired instruction set. Valid values are
       ['AVX', 'SSE42', 'FALLBACK', 'FASTEST']

    Returns
    --------
    instruction_set: integer
       An integer corresponding to the desired instruction set, as used in the
       underlying C API. The enum used here should be defined *exactly* the
       same way as the enum in ``utils/defs.h``.

    """

    msg = "Input to translate_isa_string_to_enum must be "\
          "of string type. Found type = {0}".format(type(isa))
    try:
        if not isinstance(isa, basestring):
            raise TypeError(msg)
    except NameError:
        if not isinstance(isa, str):
            raise TypeError(msg)
    valid_isa = ['FALLBACK', 'AVX', 'SSE42', 'FASTEST']
    isa_upper = isa.upper()
    if isa_upper not in valid_isa:
        msg = "Desired instruction set = {0} is not in the list of valid "\
              "instruction sets = {1}".format(isa, valid_isa)
        raise ValueError(msg)

    enums = {'FASTEST': -1,
             'FALLBACK': 0,
             'SSE': 1,
             'SSE2': 2,
             'SSE3': 3,
             'SSSE3': 4,
             'SSE4': 5,
             'SSE42': 6,
             'AVX': 7,
             'AVX2': 8,
             'AVX512F': 9
             }
    try:
        return enums[isa_upper]
    except KeyError:
        print("Do not know instruction type = {0}".format(isa))
        print("Valid instructions are {0}".format(enums.keys()))
        raise


def compute_nbins(max_diff, binsize,
                 refine_factor=1,
                 max_nbins=None):
    """
    Helper utility to find the number of bins for
    that satisfies the constraints of (binsize, refine_factor, and max_nbins).
      
    Parameters
    ------------

    max_diff : double
       Max. difference (spatial or angular) to be spanned,
       (i.e., range of allowed domain values)

    binsize : double
       Min. allowed binsize (spatial or angular) 

    refine_factor : integer, default 1
       How many times to refine the bins. The refinements occurs
       after ``nbins`` has already been determined (with ``refine_factor-1``).
       Thus, the number of bins will be **exactly** higher by 
       ``refine_factor`` compared to the base case of ``refine_factor=1``

    max_nbins : integer, default None
       Max number of allowed cells

    Returns 
    ---------

    nbins: integer, >= 1
       Number of bins that satisfies the constraints of 
       bin size >= ``binsize``, the refinement factor 
       and nbins <= ``max_nbins``.

    Example
    ---------
    
    >>> from Corrfunc.utils import compute_nbins
    >>> max_diff = 180
    >>> binsize = 10
    >>> compute_nbins(max_diff, binsize)
    18
    >>> refine_factor=2
    >>> max_nbins = 20
    >>> compute_nbins(max_diff, binsize, refine_factor=refine_factor, 
    ...              max_nbins=max_nbins)
    20

    """

    if max_diff <= 0 or binsize <= 0:
        msg = 'Error: Invalid value for max_diff = {0} or binsize = {1}. '\
              'Both must be positive'.format(max_diff, binsize)
        raise ValueError(msg)
    if max_nbins < 1:
        msg = 'Error: Invalid for the max. number of bins allowed = {0}.'\
              'Max. nbins must be >= 1'.format(max_nbins)
        raise ValueError(msg)

    if refine_factor < 1:
        msg = 'Error: Refine factor must be >=1. Found refine_factor = '\
              '{0}'.format(refine_factor)
        raise ValueError(msg)

    # At least 1 bin
    ngrid = max(1, long(max_diff/binsize))

    # Then refine
    ngrid *= refine_factor

    # But don't exceed max number of bins
    # (if passed as a parameter)
    if max_nbins:
        ngrid = min(max_nbins, ngrid)

    return ngrid             
                     
    
def gridlink_sphere(thetamax,
                    ra_limits=None,
                    dec_limits=None,
                    link_in_ra=True,
                    ra_refine_factor=1, dec_refine_factor=1,
                    max_ra_cells=100, max_dec_cells=200,
                    return_num_ra_cells=False,
                    input_in_degrees=True):
    """
    A method to optimally partition spherical regions such that pairs of 
    points within a certain angular separation, ``thetamax``, can be quickly
    computed. 

    Generates the  binning scheme used in :py:mod:`Corrfunc.mocks.DDtheta_mocks` 
    for a spherical region in Right Ascension (RA), Declination (DEC) 
    and a maximum angular separation. 

    For a given ``thetamax``, regions on the sphere are divided into bands
    in DEC bands, with the width in DEC equal to ``thetamax``. If 
    ``link_in_ra`` is set, then these DEC bands are further sub-divided 
    into RA cells. 

    Parameters
    ----------

    thetamax : double
       Max. angular separation of pairs. Expected to be in degrees 
       unless ``input_in_degrees`` is set to ``False``.

    ra_limits : array of 2 doubles. Default [0.0, 2*pi]
       Range of Righ Ascension (longitude) for the spherical region

    dec_limits : array of 2 doubles. Default [-pi/2, pi/2]
       Range of Declination (latitude) values for the spherical region

    link_in_ra : Boolean. Default True
       Whether linking in RA is done (in addition to linking in DEC)

    ra_refine_factor : integer, >= 1. Default 1
       Controls the sub-division of the RA cells. For a large number of 
       particles, higher `ra_refine_factor` typically results in a faster
       runtime

    dec_refine_factor : integer, >= 1. Default 1
       Controls the sub-division of the DEC cells. For a large number of 
       particles, higher `dec_refine_factor` typically results in a faster
       runtime

    max_ra_cells : integer, >= 1. Default 100
       The max. number of RA cells **per DEC band**.

    max_dec_cells : integer >= 1. Default 200
       The max. number of total DEC bands 

    return_num_ra_cells: bool, default False
       Flag to return the number of RA cells per DEC band

    input_in_degrees : Boolean. Default True
       Flag to show if the input quantities are in degrees. If set to 
       False, all angle inputs will be taken to be in radians.

    Returns
    ---------

    sphere_grid : A numpy compound array 
       A numpy compound array with fields ``dec_limit`` and ``ra_limit`` of 
       size 2 each. These arrays contain the beginning and end of DEC 
       and RA regions for the cell. 
    

    .. note:: If ``link_in_ra=False``, then there is effectively one RA bin
       per DEC band. The  'ra_limit' field will show the range of allowed 
       RA values.

    .. seealso:: :py:mod:`Corrfunc.mocks.DDtheta_mocks`

    Example
    --------

    >>> from Corrfunc.utils import gridlink_sphere
    >>> thetamax=30
    >>> gridlink_sphere(thetamax)
    array([([-1.57079633, -1.04719755], [ 0.        ,  3.14159265]),
       ([-1.57079633, -1.04719755], [ 3.14159265,  6.28318531]),
       ([-1.04719755, -0.52359878], [ 0.        ,  3.14159265]),
       ([-1.04719755, -0.52359878], [ 3.14159265,  6.28318531]),
       ([-0.52359878,  0.        ], [ 0.        ,  1.25663706]),
       ([-0.52359878,  0.        ], [ 1.25663706,  2.51327412]),
       ([-0.52359878,  0.        ], [ 2.51327412,  3.76991118]),
       ([-0.52359878,  0.        ], [ 3.76991118,  5.02654825]),
       ([-0.52359878,  0.        ], [ 5.02654825,  6.28318531]),
       ([ 0.        ,  0.52359878], [ 0.        ,  1.25663706]),
       ([ 0.        ,  0.52359878], [ 1.25663706,  2.51327412]),
       ([ 0.        ,  0.52359878], [ 2.51327412,  3.76991118]),
       ([ 0.        ,  0.52359878], [ 3.76991118,  5.02654825]),
       ([ 0.        ,  0.52359878], [ 5.02654825,  6.28318531]),
       ([ 0.52359878,  1.04719755], [ 0.        ,  3.14159265]),
       ([ 0.52359878,  1.04719755], [ 3.14159265,  6.28318531]),
       ([ 1.04719755,  1.57079633], [ 0.        ,  3.14159265]),
       ([ 1.04719755,  1.57079633], [ 3.14159265,  6.28318531])], 
      dtype=[(u'dec_limit', '<f8', (2,)), (u'ra_limit', '<f8', (2,))])
    >>> gridlink_sphere(60, dec_refine_factor=3, ra_refine_factor=2)
    array([([-1.57079633, -1.22173048], [ 0.        ,  1.57079633]),
           ([-1.57079633, -1.22173048], [ 1.57079633,  3.14159265]),
           ([-1.57079633, -1.22173048], [ 3.14159265,  4.71238898]),
           ([-1.57079633, -1.22173048], [ 4.71238898,  6.28318531]),
           ([-1.22173048, -0.87266463], [ 0.        ,  1.57079633]),
           ([-1.22173048, -0.87266463], [ 1.57079633,  3.14159265]),
           ([-1.22173048, -0.87266463], [ 3.14159265,  4.71238898]),
           ([-1.22173048, -0.87266463], [ 4.71238898,  6.28318531]),
           ([-0.87266463, -0.52359878], [ 0.        ,  1.57079633]),
           ([-0.87266463, -0.52359878], [ 1.57079633,  3.14159265]),
           ([-0.87266463, -0.52359878], [ 3.14159265,  4.71238898]),
           ([-0.87266463, -0.52359878], [ 4.71238898,  6.28318531]),
           ([-0.52359878, -0.17453293], [ 0.        ,  1.57079633]),
           ([-0.52359878, -0.17453293], [ 1.57079633,  3.14159265]),
           ([-0.52359878, -0.17453293], [ 3.14159265,  4.71238898]),
           ([-0.52359878, -0.17453293], [ 4.71238898,  6.28318531]),
           ([-0.17453293,  0.17453293], [ 0.        ,  1.57079633]),
           ([-0.17453293,  0.17453293], [ 1.57079633,  3.14159265]),
           ([-0.17453293,  0.17453293], [ 3.14159265,  4.71238898]),
           ([-0.17453293,  0.17453293], [ 4.71238898,  6.28318531]),
           ([ 0.17453293,  0.52359878], [ 0.        ,  1.57079633]),
           ([ 0.17453293,  0.52359878], [ 1.57079633,  3.14159265]),
           ([ 0.17453293,  0.52359878], [ 3.14159265,  4.71238898]),
           ([ 0.17453293,  0.52359878], [ 4.71238898,  6.28318531]),
           ([ 0.52359878,  0.87266463], [ 0.        ,  1.57079633]),
           ([ 0.52359878,  0.87266463], [ 1.57079633,  3.14159265]),
           ([ 0.52359878,  0.87266463], [ 3.14159265,  4.71238898]),
           ([ 0.52359878,  0.87266463], [ 4.71238898,  6.28318531]),
           ([ 0.87266463,  1.22173048], [ 0.        ,  1.57079633]),
           ([ 0.87266463,  1.22173048], [ 1.57079633,  3.14159265]),
           ([ 0.87266463,  1.22173048], [ 3.14159265,  4.71238898]),
           ([ 0.87266463,  1.22173048], [ 4.71238898,  6.28318531]),
           ([ 1.22173048,  1.57079633], [ 0.        ,  1.57079633]),
           ([ 1.22173048,  1.57079633], [ 1.57079633,  3.14159265]),
           ([ 1.22173048,  1.57079633], [ 3.14159265,  4.71238898]),
           ([ 1.22173048,  1.57079633], [ 4.71238898,  6.28318531])], 
          dtype=[(u'dec_limit', '<f8', (2,)), (u'ra_limit', '<f8', (2,))])

    """

    from math import radians, pi
    import numpy as np
    

    if input_in_degrees:
        thetamax = radians(thetamax)
        for x in [ra_limits, dec_limits]:
            if x:
                x = radians(x)
            
    if not ra_limits:
        ra_limits = [0.0, 2.0*pi]
        
    if not dec_limits:
        dec_limits = [-0.5*pi, 0.5*pi]

    if dec_limits[0] >= dec_limits[1]:
        msg = 'Declination limits should be sorted in increasing '\
              'order. However, dec_limits = [{0}, {1}] is not'.\
              format(dec_limits[0], dec_limits[1])
        raise ValueError(msg)

    if ra_limits[0] >= ra_limits[1]:
        msg = 'Declination limits should be sorted in increasing '\
              'order. However, ra_limits = [{0}, {1}] is not'.\
              format(ra_limits[0], ra_limits[1])
        raise ValueError(msg)

    if dec_limits[0] < -0.5*pi or dec_limits[1] > 0.5*pi:
        msg = 'Valid range of values for declination are [-pi/2, +pi/2] deg. '\
              'However, dec_limits = [{0}, {1}] does not fall within that '\
              'range'.format(dec_limits[0], dec_limits[1])
        raise ValueError(msg)
    
    if ra_limits[0] < 0.0 or ra_limits[1] > 2.0*pi:
        msg = 'Valid range of values for declination are [0.0, 2*pi] deg. '\
              'However, ra_limits = [{0}, {1}] does not fall within that '\
              'range'.format(ra_limits[0], ra_limits[1])
        raise ValueError(msg)
    
    dec_diff = abs(dec_limits[1] - dec_limits[0])
    ngrid_dec = compute_nbins(dec_diff, thetamax,
                             refine_factor=dec_refine_factor,
                             max_nbins=max_dec_cells)
    
    dec_binsize = dec_diff/ngrid_dec

    # Upper and lower limits of the declination bands
    grid_dtype= np.dtype({'names':['dec_limit','ra_limit'],
                          'formats':[(np.float, (2, )), (np.float, (2, ))]
    })
    if not link_in_ra:
        sphere_grid = np.zeros(ngrid_dec, dtype=grid_dtype)
        for i, r in enumerate(sphere_grid['dec_limit']):
            r[0] = dec_limits[0] + i*dec_binsize
            r[1] = dec_limits[0] + (i+1)*dec_binsize

        for r in sphere_grid['ra_limit']:
            r[0] = ra_limits[0]
            r[1] = ra_limits[1]

        return sphere_grid

    # RA linking is requested
    ra_diff = ra_limits[1] - ra_limits[0]
    sin_half_thetamax = np.sin(thetamax)
    costhetamax = np.cos(thetamax)
    max_nmesh_ra = 1

    totncells = 0
    num_ra_cells = np.zeros(ngrid_dec, dtype=np.int64)
    num_ra_cells[:] = ra_refine_factor
    # xrange is replaced by range for python3
    # by using a try/except at the top
    for idec in xrange(ngrid_dec):
        dec_min = dec_limits[0] + idec*dec_binsize
        dec_max = dec_min + dec_binsize

        cos_dec_min = np.cos(dec_min)
        cos_dec_max = np.cos(dec_max)

        if cos_dec_min < cos_dec_max:
            min_cos = cos_dec_min
        else:
            min_cos = cos_dec_max

        if min_cos > 0:
            _tmp = sin_half_thetamax/min_cos
            # clamp to range [0.0, 1.0]
            _tmp = max(min(_tmp, 1.0), 0.0)
            ra_binsize = min(2.0 * np.arcsin(_tmp), ra_diff)
            num_ra_cells[idec] = compute_nbins(ra_diff, ra_binsize,
                                              refine_factor=ra_refine_factor,
                                              max_nbins=max_ra_cells)
            
    totncells = num_ra_cells.sum()
    sphere_grid = np.zeros(totncells, dtype=grid_dtype)
    ra_binsizes = ra_diff/num_ra_cells

    start = 0
    for idec in xrange(ngrid_dec):
        assert start + num_ra_cells[idec] <= totncells
        source_sel = np.s_[start:start+num_ra_cells[idec]]
        for ira, r in enumerate(sphere_grid[source_sel]):
            r['dec_limit'][0] = dec_limits[0] + dec_binsize*idec
            r['dec_limit'][1] = dec_limits[0] + dec_binsize*(idec + 1)
            r['ra_limit'][0] = ra_limits[0] + ra_binsizes[idec] * ira
            r['ra_limit'][1] = ra_limits[0] + ra_binsizes[idec] * (ira + 1)
            
        start += num_ra_cells[idec]

    if return_num_ra_cells:
        return sphere_grid, num_ra_cells
    else:
        return sphere_grid

if __name__ == '__main__':
    import doctest
    doctest.testmod()
