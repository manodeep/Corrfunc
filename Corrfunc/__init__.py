#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import sys

__all__ = ["_countpairs", "_countpairs_mocks", "utils"]

# from Corrfunc import * throws: TypeError: Item in ``from list'' not a string
# following the accepted answer in:
# http://stackoverflow.com/questions/19913653/no-unicode-in-all-
# for-a-packages-init
if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]

__version__ = "1.2.0"


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


def write_text_file(filename, contents, encoding="utf-8"):
    """
    Writes a file under python3 with encoding (default UTF-8).
    Also works under python2, without encoding.
    Uses the EAFP (https://docs.python.org/2/glossary.html#term-eafp)
    principle.
    """
    try:
        with open(filename, 'w', encoding) as f:
            f.write(contents)
            
    except TypeError:
        with open(filename, 'w') as f:
            f.write(contents)


