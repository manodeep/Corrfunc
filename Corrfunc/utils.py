#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

__all__ = ['rd']

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
