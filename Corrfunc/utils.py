#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

__all__ = ['rd']

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
