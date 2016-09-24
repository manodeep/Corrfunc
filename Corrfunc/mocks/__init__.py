#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wrapper for all clustering statistic calculations on galaxies
in a mock catalog.
"""
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys

__author__ = ('Manodeep Sinha')
__all__ = ("DDrppi_mocks", "DDtheta_mocks", "vpf_mocks", )

if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]
