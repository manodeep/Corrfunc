#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.utils import bytes_to_native_str
import sys
import os.path as path

__all__ = ['rd', 'read_catalog']

if sys.version_info[0] >= 3:
    def rd(filename):
        """
        Reads a file under python3 assuming an UTF-8 encoding.
        """
        with open(filename, encoding="utf-8") as f:
            r = f.read()

        return r
else:
    def rd(filename):
        """
        Reads a file under python2.
        """
        with open(filename) as f:
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
        if return_dtype is None:
            msg = "Return data-type must be set and a valid numpy data-type"
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
        filename = path.join(path.dirname(path.abspath(__file__)),
                             "../xi_theory/tests/data/", "gals_Mr19")
        # Figure out the datatype, use the header file in the include directory
        # because that is most likely correct (common.mk might have been
        # modified but not recompiled)
        include_file = path.join(path.dirname(path.abspath(__file__)),
                                 "../include/", "countpairs.h")
        includes = rd(include_file)
        vector_type = re.search(r'(\w+)\s*\*\s*rupp\s*\;', includes,
                                re.I).group(1)
        allowed_types = {"float": np.float32, "double": np.float}
        if vector_type not in list(allowed_types.keys()):
            print("Error: Unknown precision={0} found in header file {1}. \
            Allowed types are `{2}'".format(vector_type,
                                            include_file,
                                            allowed_types))
            sys.exit()

        dtype = allowed_types[vector_type]
        allowed_exts = {'.ff': read_fastfood,
                        '.txt': read_ascii,
                        '.dat': read_ascii,
                        '.csv': read_ascii
                        }

        for e in allowed_exts:
            if path.exists(filename + e):
                f = allowed_exts[e]
                x, y, z = f(filename + e, dtype)
                return x, y, z

        raise IOError("Could not locate {0} with any of these extensions \
        = {1}".format(filename, allowed_exts.keys()))
    else:
        # Likely an user-supplied value
        if path.exists(filebase):
            extension = path.splitext(filebase)[1]
            f = read_fastfood if '.ff' in extension else read_ascii

            # default return is double
            x, y, z = f(filebase, np.float32)
            return x, y, z

        raise IOError("Could not locate file {0}".format(filebase))
