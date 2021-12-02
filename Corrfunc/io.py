#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Routines to read galaxy catalogs from disk.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from os.path import dirname, abspath, splitext, exists as file_exists,\
    join as pjoin
import numpy as np
try:
    import pandas as pd
except ImportError:
    pd = None


__all__ = ('read_fastfood_catalog', 'read_ascii_catalog', 'read_catalog')


def read_fastfood_catalog(filename, return_dtype=None, need_weights=None):
    """
    Read a galaxy catalog from a fast-food binary file.

    Parameters
    -----------
    filename: string
        Filename containing the galaxy positions

    return_dtype: numpy dtype for returned arrays. Default ``numpy.float``
        Specifies the datatype for the returned arrays. Must be in
        {np.float64, np.float32}

    need_weights: boolean, default None.
        Returns weight array in addition to the X/Y/Z arrays.

    Returns
    --------

    X, Y, Z: numpy arrays
        Returns the triplet of X/Y/Z positions as separate numpy arrays.
        

    Example
    --------
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.io import read_fastfood_catalog
    >>> filename = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                  "../theory/tests/data/",
    ...                  "gals_Mr19.ff")
    >>> X, Y, Z = read_fastfood_catalog(filename)
    >>> N = 20
    >>> for x,y,z in zip(X[0:N], Y[0:N], Z[0:]):
    ...     print("{0:10.5f} {1:10.5f} {2:10.5f}".format(x, y, z))
    ...     # doctest: +NORMALIZE_WHITESPACE
    419.94550    1.96340    0.01610
    419.88272    1.79736    0.11960
    0.32880   10.63620    4.16550
    0.15314   10.68723    4.06529
    0.46400    8.91150    6.97090
    6.30690    9.77090    8.61080
    5.87160    9.65870    9.29810
    8.06210    0.42350    4.89410
    11.92830    4.38660    4.54410
    11.95543    4.32622    4.51485
    11.65676    4.34665    4.53181
    11.75739    4.26262    4.31666
    11.81329    4.27530    4.49183
    11.80406    4.54737    4.26824
    12.61570    4.14470    3.70140
    13.23640    4.34750    5.26450
    13.19833    4.33196    5.29435
    13.21249    4.35695    5.37418
    13.06805    4.24275    5.35126
    13.19693    4.37618    5.28772

    """
    if return_dtype is None:
        return_dtype = np.float64

    if return_dtype not in [np.float32, np.float64]:
        msg = "Return data-type must be set and a valid numpy float"
        raise ValueError(msg)

    if not file_exists(filename):
        msg = "Could not find file = {0}".format(filename)
        raise IOError(msg)

    import struct
    try:
        from future.utils import bytes_to_native_str
    except ImportError:
        print("\n\tPlease run python -m pip install Corrfunc before using "
              "the 'Corrfunc' package\n")
        raise

    with open(filename, "rb") as f:
        skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
        idat = struct.unpack(bytes_to_native_str(b'@iiiii'),
                             f.read(20))[0:5]
        skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
        assert skip1 == 20 and skip2 == 20,\
            "fast-food file seems to be incorrect (reading idat)"
        ngal = idat[1]

        fdat_bytes = 4 + 36 + 4
        znow_bytes = 4 + 4 + 4
        # seek over the fdat + znow fields + padding bytes
        # from current position
        f.seek(fdat_bytes + znow_bytes, 1)

        # read the padding bytes for the x-positions
        skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
        assert skip1 == ngal * 4 or skip1 == ngal * 8, \
            "fast-food file seems to be corrupt (padding bytes)"

        # seek back 4 bytes from current position
        f.seek(-4, 1)
        pos = {}
        for field in 'xyz' + ('w' if need_weights else ''):
            skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            assert skip1 == ngal * 4 or skip1 == ngal * 8, \
                "fast-food file seems to be corrupt (padding bytes a)"
            # the next division must be the integer division
            input_dtype = np.float32 if skip1 // ngal == 4 else np.float64
            array = np.fromfile(f, input_dtype, ngal)
            skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            if return_dtype == input_dtype:
                pos[field] = array
            else:
                pos[field] = [return_dtype(a) for a in array]

    toret = [np.array(pos[name]) for name in ['x','y','z']]
    if need_weights:
        toret.append(np.array(pos['w']))

    return toret


def read_ascii_catalog(filename, return_dtype=None):
    """
    Read a galaxy catalog from an ascii file.

    Parameters
    -----------
    filename: string
        Filename containing the galaxy positions

    return_dtype: numpy dtype for returned arrays. Default ``numpy.float``
        Specifies the datatype for the returned arrays. Must be in
        {np.float64, np.float32}

    Returns
    --------

    X, Y, Z: numpy arrays
        Returns the triplet of X/Y/Z positions as separate numpy arrays.

    Example
    --------
    >>> from __future__ import print_function
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.io import read_ascii_catalog
    >>> filename = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../mocks/tests/data/", "Mr19_mock_northonly.rdcz.dat")
    >>> ra, dec, cz = read_ascii_catalog(filename)
    >>> N = 20
    >>> for r,d,c in zip(ra[0:N], dec[0:N], cz[0:]):
    ...     print("{0:10.5f} {1:10.5f} {2:10.5f}".format(r, d, c))
    ...     # doctest: +NORMALIZE_WHITESPACE
    178.45087   67.01112 19905.28514
    178.83495   67.72519 19824.02285
    179.50132   67.67628 19831.21553
    182.75497   67.13004 19659.79825
    186.29853   68.64099 20030.64412
    186.32346   68.65879 19763.38137
    187.36173   68.15151 19942.66996
    187.20613   68.56189 19996.36607
    185.56358   67.97724 19729.32308
    183.27930   67.11318 19609.71345
    183.86498   67.82823 19500.44130
    184.07771   67.43429 19440.53790
    185.13370   67.15382 19390.60304
    189.15907   68.28252 19858.85853
    190.12209   68.55062 20044.29744
    193.65245   68.36878 19445.62469
    194.93514   68.34870 19158.93155
    180.36897   67.50058 18671.40780
    179.63278   67.51318 18657.59191
    180.75742   67.95530 18586.88913

    """

    if return_dtype is None:
        return_dtype = np.float64

    if not file_exists(filename):
        msg = "Could not find file = {0}".format(filename)
        raise IOError(msg)

    # check if pandas is available - much faster to read in the data
    # using pandas
    if pd is not None:
        df = pd.read_csv(filename, header=None,
                         engine="c",
                         dtype={"x": return_dtype,
                                "y": return_dtype,
                                "z": return_dtype},
                         delim_whitespace=True)
        x = np.asarray(df[0], dtype=return_dtype)
        y = np.asarray(df[1], dtype=return_dtype)
        z = np.asarray(df[2], dtype=return_dtype)

    else:
        x, y, z, _ = np.genfromtxt(filename, dtype=return_dtype, unpack=True)

    return x, y, z


def read_catalog(filebase=None, return_dtype=np.float64):
    """
    Reads a galaxy/randoms catalog and returns 3 XYZ arrays.

    Parameters
    -----------

    filebase: string (optional)
        The fully qualified path to the file. If omitted, reads the
        theory galaxy catalog under ../theory/tests/data/

    return_dtype: numpy dtype for returned arrays. Default ``numpy.float``
        Specifies the datatype for the returned arrays. Must be in
        {np.float64, np.float32}

    Returns
    --------

    ``x y z`` - Unpacked numpy arrays compatible with the installed
    version of ``Corrfunc``.


    .. note:: If the filename is omitted, then first the fast-food file
        is searched for, and then the ascii file. End-users should always
        supply the full filename.


    """

    if filebase is None:
        filename = pjoin(dirname(abspath(__file__)),
                         "../theory/tests/data/", "gals_Mr19")
        allowed_exts = {'.ff': read_fastfood_catalog,
                        '.txt': read_ascii_catalog,
                        '.dat': read_ascii_catalog,
                        '.csv': read_ascii_catalog
                        }

        for e in allowed_exts:
            if file_exists(filename + e):
                f = allowed_exts[e]
                x, y, z = f(filename + e, return_dtype)
                return x, y, z

        raise IOError("Could not locate {0} with any of these extensions \
        = {1}".format(filename, allowed_exts.keys()))
    else:
        # Likely an user-supplied value
        if file_exists(filebase):
            extension = splitext(filebase)[1]
            f = read_fastfood_catalog if '.ff' in extension else read_ascii_catalog

            # default return is double
            x, y, z = f(filebase, return_dtype)
            return x, y, z

        raise IOError("Could not locate file {0}".format(filebase))

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
