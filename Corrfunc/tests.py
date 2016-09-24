#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys


__all__ = ['tests', ]
if sys.version_info[0] < 3:
    __all__ = [n.encode('ascii') for n in __all__]


def tests():
    """
    Wrapper to run the two scripts that should have been installed
    with the Corrfunc package.

    If the two scripts (one for theory extensions, one for mocks extensions)
    run successfully, then the package is working correctly.
    """

    # Import the script for calling the theory extensions
    from Corrfunc import call_correlation_functions as ct
    # Import the script for calling the mocks extensions
    from Corrfunc import call_correlation_functions_mocks as cm

    # Run the theory script
    ct.main()

    # Run the mocks script
    cm.main()


if __name__ == '__main__':
    tests()
