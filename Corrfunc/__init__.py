#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                                                unicode_literals)
import sys

__all__ = ['_countpairs','_countpairs_mocks','utils']
__version__ = "0.0.1"

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
