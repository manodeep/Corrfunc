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

__version__ = "1.1.0"

if sys.version_info[0] >= 3:
    def rd(filename):
        with open(filename, encoding="utf-8") as f:
            r = f.read()

        return r
else:
    def rd(filename):
        with open(filename) as f:
            r = f.read()

        return r
