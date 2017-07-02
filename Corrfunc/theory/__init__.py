#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrapper for all clustering statistic calculations on galaxies
in a simulation volume.
"""
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = ('Manodeep Sinha')
__all__ = ('DD', 'DDrppi', 'wp', 'xi', 'vpf', 'DDsmu',)

import sys

from .DD import DD
from .DDrppi import DDrppi
from .wp import wp
from .xi import xi
from .vpf import vpf
from .DDsmu import DDsmu

if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]
