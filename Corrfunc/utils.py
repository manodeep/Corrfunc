#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.utils import bytes_to_native_str
import os

__all__ = ['read_text_file', 'read_catalog']


def read_text_file(filename, encoding="utf-8"):
    """
    Reads a file under python3 with encoding (default UTF-8).
    Also works under python2, without encoding.
    Uses the EAFP (https://docs.python.org/2/glossary.html#term-eafp)
    principle.
    """
    try:
        with open(filename, 'r', encoding) as f:
            r = f.read()
    except TypeError:
        with open(filename, 'r') as f:
            r = f.read()

    return r


def read_catalog(filebase=None):
    """
    Reads a galaxy/randoms catalog.

    :param filebase: (optional)
        The fully qualified path to the file. If omitted, reads the
        theory galaxy catalog under ../xi_theory/tests/data/

    Returns:
    * ``x y z`` - Unpacked numpy arrays compatible with the installed
    version of ``Corrfunc``.

    **Note** If the filename is omitted, then first the fast-food file
    is searched for, and then the ascii file. End-users should always
    supply the full filename.
    """

    import re
    import numpy as np

    def read_ascii(filename, return_dtype=None):
        if return_dtype is None:
            msg = 'Return data-type must be set and a valid numpy data-type'
            raise ValueError(msg)

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

    def read_fastfood(filename, return_dtype=None, need_header=None):
        if return_dtype is None or return_dtype not in [np.float32, np.float]:
            msg = "Return data-type must be set and a valid numpy float"
            raise ValueError(msg)

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
                pos[field] = array if return_dtype == input_dtype \
                             else [return_dtype(a) for a in array]

        x = pos['x']
        y = pos['y']
        z = pos['z']

        if need_header is not None:
            return idat, fdat, znow, x, y, z
        else:
            return x, y, z

    if filebase is None:
        filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "../xi_theory/tests/data/", "gals_Mr19")
        # Figure out the datatype, use the header file in the include directory
        # because that is most likely correct (common.mk might have been
        # modified but not recompiled)
        include_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "../include/", "countpairs.h")
        includes = read_text_file(include_file)
        vector_type = re.search(r'(\w+)\s*\*\s*rupp\s*\;', includes,
                                re.I).group(1)
        allowed_types = {"float": np.float32, "double": np.float}
        if vector_type not in list(allowed_types.keys()):
            msg = "Error: Unknown precision={0} found in header file {1}. \
            Allowed types are `{2}'".format(vector_type,
                                            include_file,
                                            allowed_types)
            raise AssertionError(msg)

        dtype = allowed_types[vector_type]
        allowed_exts = {'.ff': read_fastfood,
                        '.txt': read_ascii,
                        '.dat': read_ascii,
                        '.csv': read_ascii
                        }

        for e in allowed_exts:
            if os.path.exists(filename + e):
                f = allowed_exts[e]
                x, y, z = f(filename + e, dtype)
                return x, y, z

        raise IOError("Could not locate {0} with any of these extensions \
        = {1}".format(filename, allowed_exts.keys()))
    else:
        # Likely an user-supplied value
        if os.path.exists(filebase):
            extension = os.path.splitext(filebase)[1]
            f = read_fastfood if '.ff' in extension else read_ascii

            # default return is double
            x, y, z = f(filebase, np.float32)
            return x, y, z

        raise IOError("Could not locate file {0}".format(filebase))
