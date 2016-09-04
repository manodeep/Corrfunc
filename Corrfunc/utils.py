#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
try:
    from future.utils import bytes_to_native_str
except ImportError:
    print("\n\tPlease run python setup.py install before using "
          "the 'Corrfunc' package\n")
    raise

from os.path import dirname, abspath, splitext, exists as file_exists,\
    join as pjoin

__all__ = ['translate_isa_string_to_enum', 'read_text_file', 'read_catalog',
           'return_file_with_rbins', 'fix_ra_dec', 'fix_cz']
if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]


def return_file_with_rbins(rbins):
    """
    Helper function to ensure that the ``binfile`` required by the Corrfunc
    extensions is a actually a string (expecting that to be a filename).

    Checks if the input is a string and file; return if True. If not, and
    the input is an array, then a temporary file is created and the contents
    of rbins is written out.
    
    Returns a filename containing the bins and a flag stating if the file needs
    to be deleted after use.
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
        with tempfile.NamedTemporaryFile(delete=False) as f:
            for i in xrange(len(rbins)-1):
                f.write("{0} {1}\n".format(rbins[i], rbins[i+1]))

            tmpfilename = f.name

        delete_after_use = True
        return tmpfilename, delete_after_use

    msg = "Input `binfile` was not a valid array (>= 1 element)."\
          "Num elements = {0}".format(len(rbins))
    raise TypeError(msg)

    
def fix_cz(cz):
    """
    Multiplies the input array by speed of light, if the input values are
    too small. Essentially, converts redshift into `cz`, if the user passed
    redshifts instead of `cz`.
    """

    # if I find that max cz is smaller than this threshold,
    # then I will assume z has been supplied rather than cz
    max_cz_threshold = 10.0
    if max(cz) < max_cz_threshold:
        speed_of_light = 299800.0
        cz *= speed_of_light

    return cz


def fix_ra_dec(ra, dec):
    """
    Wraps input RA and DEC values into range expected by the extensions.
    
    RA: In range [0, 360.0]
    DEC: In range [-90.0, 90.0]
    """

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
    

def translate_isa_string_to_enum(isa):
    """
    Helper function to convert an user-supplied string to the
    underlying enum in the C-API. The extensions only have specific
    implementations for AVX, SSE42 and FALLBACK. Any other value,
    will raise a ValueError.
    
    The enum definition contains all the valid instruction sets I know
    of. Here to facilitate easy extension when a new instruction set has
    been added in.

    Parameters:
    ----------
    isa: string
       A string containing the desired instruction set. Valid values are
       ['AVX', 'SSE42', 'FALLBACK', 'FASTEST']
    
    Returns:
    -------
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

    def read_fastfood(filename, return_dtype=np.float, need_header=None):
        if return_dtype not in [np.float32, np.float]:
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

    if filebase is None:
        filename = pjoin(dirname(abspath(__file__)),
                         "../xi_theory/tests/data/", "gals_Mr19")
        dtype = np.float
        allowed_exts = {'.ff': read_fastfood,
                        '.txt': read_ascii,
                        '.dat': read_ascii,
                        '.csv': read_ascii
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
            f = read_fastfood if '.ff' in extension else read_ascii

            # default return is double
            x, y, z = f(filebase, np.float)
            return x, y, z

        raise IOError("Could not locate file {0}".format(filebase))
