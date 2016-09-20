#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Routines to read galaxy catalogs from disk.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
from os.path import dirname, abspath, splitext, exists as file_exists,\
    join as pjoin
import numpy as np
try:
    from future.utils import bytes_to_native_str
except ImportError:
    print("\n\tPlease run python setup.py install before using "
          "the 'Corrfunc' package\n")
    raise


__all__ = ['read_catalog', 'read_fastfood_catalog', 'read_ascii_catalog', ]
if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]


def read_fastfood_catalog(filename, return_dtype=None, need_header=None):
    """
    Read a galaxy catalog from a fast-food binary file.

    Parameters
    -----------
    filename: string
        Filename containing the galaxy positions

    return_dtype: numpy dtype for returned arrays. Default ``numpy.float``
        Specifies the datatype for the returned arrays. Must be in
        {np.float, np.float32}

    need_header: boolean, default None.
        Returns the header found in the fast-food file in addition to the
        X/Y/Z arrays.

    Returns
    --------

    X, Y, Z: numpy arrays
        Returns the triplet of X/Y/Z positions as separate numpy arrays.

        If need_header is set, then the header is also returned

    Example
    --------
    >>> from Corrfunc.io import read_fastfood_catalog
    >>> filename = "galaxies.ff"
    >>> X, Y, Z = read_fastfood_catalog(filename, dtype=numpy.float32)


    """
    if return_dtype is None:
        return_dtype = np.float

    if return_dtype not in [np.float32, np.float]:
        msg = "Return data-type must be set and a valid numpy float"
        raise ValueError(msg)

    if not file_exists(filename):
        msg = "Could not find file = {0}".format(filename)
        raise IOError(msg)

    import struct
    with open(filename, "rb") as f:
        skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
        idat = struct.unpack(bytes_to_native_str(b'@iiiii'),
                             f.read(20))[0:5]
        skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
        assert skip1 == 20 and skip2 == 20,\
            "fast-food file seems to be incorrect (reading idat)"
        ngal = idat[1]

        if need_header is not None:
            # now read fdat
            skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            fdat = struct.unpack(bytes_to_native_str(b'@fffffffff'),
                                 f.read(36))[0:9]
            skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            assert skip1 == 36 and skip2 == 36,\
                "fast-food file seems to be incorrect (reading fdat )"

            skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            znow = struct.unpack(bytes_to_native_str(b'@f'), f.read(4))[0]
            skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            assert skip1 == 4 and skip2 == 4,\
                "fast-food file seems to be incorrect (reading redshift)"
        else:
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
        for field in 'xyz':
            skip1 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            assert skip1 == ngal * 4 or skip1 == ngal * 8, \
                "fast-food file seems to be corrupt (padding bytes a)"
            # the next division must be the integer division
            input_dtype = np.float32 if skip1 // ngal == 4 else np.float
            array = np.fromfile(f, input_dtype, ngal)
            skip2 = struct.unpack(bytes_to_native_str(b'@i'), f.read(4))[0]
            if return_dtype == input_dtype:
                pos[field] = array
            else:
                pos[field] = [return_dtype(a) for a in array]

    x = np.array(pos['x'])
    y = np.array(pos['y'])
    z = np.array(pos['z'])

    if need_header is not None:
        return idat, fdat, znow, x, y, z
    else:
        return x, y, z


def read_ascii_catalog(filename, return_dtype=None):
    """
    Read a galaxy catalog from an ascii file.

    Parameters
    -----------
    filename: string
        Filename containing the galaxy positions

    return_dtype: numpy dtype for returned arrays. Default ``numpy.float``
        Specifies the datatype for the returned arrays. Must be in
        {np.float, np.float32}

    Returns
    --------

    X, Y, Z: numpy arrays
        Returns the triplet of X/Y/Z positions as separate numpy arrays.

    Example
    --------
    >>> from Corrfunc.io import read_ascii_catalog
    >>> filename = "galaxies.dat"
    >>> X, Y, Z = read_ascii(filename, dtype=numpy.float)

    """

    if return_dtype is None:
        msg = 'Return data-type must be set and a valid numpy data-type'
        raise ValueError(msg)

    if not file_exists(filename):
        msg = "Could not find file = {0}".format(filename)
        raise IOError(msg)

    # check if pandas is available - much faster to read in the data
    # using pandas
    print("Reading in the data...")
    try:
        import pandas as pd
    except ImportError:
        pd = None

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
        x, y, z = np.genfromtxt(filename, dtype=return_dtype, unpack=True)

    return x, y, z


def read_catalog(filebase=None):
    """
    Reads a galaxy/randoms catalog and returns 3 XYZ arrays.

    Parameters
    -----------

    filebase: string (optional)
        The fully qualified path to the file. If omitted, reads the
        theory galaxy catalog under ../theory/tests/data/

    Returns
    --------

    ``x y z`` - Unpacked numpy arrays compatible with the installed
    version of ``Corrfunc``.

    **Note** If the filename is omitted, then first the fast-food file
    is searched for, and then the ascii file. End-users should always
    supply the full filename.
    """

    if filebase is None:
        filename = pjoin(dirname(abspath(__file__)),
                         "../theory/tests/data/", "gals_Mr19")
        dtype = np.float
        allowed_exts = {'.ff': read_fastfood_catalog,
                        '.txt': read_ascii_catalog,
                        '.dat': read_ascii_catalog,
                        '.csv': read_ascii_catalog
                        }

        for e in allowed_exts:
            if file_exists(filename + e):
                f = allowed_exts[e]
                x, y, z = f(filename + e, dtype)
                return x, y, z

        raise IOError("Could not locate {0} with any of these extensions \
        = {1}".format(filename, allowed_exts.keys()))
    else:
        # Likely an user-supplied value
        if file_exists(filebase):
            extension = splitext(filebase)[1]
            f = read_fastfood_catalog if '.ff' in extension else read_ascii_catalog

            # default return is double
            x, y, z = f(filebase, np.float)
            return x, y, z

        raise IOError("Could not locate file {0}".format(filebase))
