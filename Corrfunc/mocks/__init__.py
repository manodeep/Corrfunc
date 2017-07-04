#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrapper for all clustering statistic calculations on galaxies
in a mock catalog.
"""
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ("DDrppi_mocks", "DDtheta_mocks", "vpf_mocks", "DDsmu_mocks" )

import sys
from .DDrppi_mocks import DDrppi_mocks
from .DDtheta_mocks import DDtheta_mocks
from .vpf_mocks import vpf_mocks
from .DDsmu_mocks import DDsmu_mocks

if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]
