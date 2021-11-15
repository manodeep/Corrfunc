#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import numpy as np

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


def allclose(a, *others):
    toret = True
    for b in others:
        for name in ['npairs','weightavg']:
            if not np.all(np.allclose(a[name], b[name])):
                print("Mis-match for {0}: (a, b) = {1} ".format(name, list(zip(a[name], b[name]))))
            toret &= np.allclose(a[name], b[name])

    return toret


def test_linear_binning_mocks(isa='fallback'):
    """Here we test that the special treatment for linear binning returns the correct result."""
    from os.path import dirname, abspath, join as pjoin
    import time
    import Corrfunc
    from Corrfunc.io import read_catalog
    from Corrfunc.mocks import DDrppi_mocks, DDsmu_mocks, DDtheta_mocks
    tstart = time.time()
    filename = pjoin(dirname(abspath(__file__)),
                     "../mocks/tests/data/", "Mr19_mock_northonly.rdcz.ff")

    t0 = time.time()
    ra, dec, cz = np.array(read_catalog(filename))[:,:1000]
    rng = np.random.RandomState(seed=42)
    w = rng.uniform(0.5, 1., (1, len(ra)))
    t1 = time.time()
    print("RA min  = {0} max = {1}".format(np.min(ra), np.max(ra)))
    print("DEC min = {0} max = {1}".format(np.min(dec), np.max(dec)))
    print("cz min  = {0} max = {1}".format(np.min(cz), np.max(cz)))
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")

    nthreads = 4
    pimax = 40.0
    binfile = np.linspace(5.1, 50.1, 46)
    autocorr = 1
    numbins_to_print = 5
    cosmology = 1
    bin_types = ['custom', 'lin']

    def test_bin_type(pair_counter, *args, **kwargs):
        res = (pair_counter(*args, bin_type=bin_type,
               weight_type="pair_product", **kwargs) for bin_type in bin_types)
        assert allclose(*res)

    print("\nRunning 2-D correlation function xi(rp,pi)")
    test_bin_type(DDrppi_mocks, autocorr, cosmology, nthreads,
                  pimax, binfile,
                  ra, dec, cz, weights1=w,
                  output_rpavg=True, verbose=True,
                  isa=isa)

    nmu_bins = 10
    mu_max = 1.0
    print("\nRunning 2-D correlation function xi(s,mu)")
    test_bin_type(DDsmu_mocks, autocorr, cosmology, nthreads,
                  mu_max, nmu_bins, binfile,
                  ra, dec, cz, weights1=w,
                  output_savg=True, verbose=True,
                  isa=isa)

    binfile = np.linspace(3.1, 13.1, 11)
    print("\nRunning angular correlation function DD(theta)")
    test_bin_type(DDtheta_mocks, autocorr, nthreads, binfile,
                  ra, dec, RA2=ra, DEC2=dec,
                  weights1=w,
                  weights2=w,
                  output_thetaavg=True, fast_acos=False,
                  verbose=True, isa=isa)


def test_linear_binning_theory(isa='fallback'):
    """
    Here we test that the special treatment for linear binning
    returns the correct result.
    """
    from os.path import dirname, abspath, join as pjoin
    import time
    import numpy as np
    import Corrfunc
    from Corrfunc.io import read_catalog
    from Corrfunc.theory import DD, DDrppi, DDsmu, wp, xi

    tstart = time.time()
    t0 = tstart
    x, y, z = np.array(read_catalog())[:,:1000]
    rng = np.random.RandomState(seed=42)
    w = rng.uniform(0.5, 1., (1, len(x)))
    boxsize = 420.0
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))

    numbins_to_print = 5

    print("Beginning Theory Correlation functions calculations")
    nthreads = 4
    pimax = 40.0
    binfile = np.linspace(5.1, 30.1, 26)
    autocorr = 1
    periodic = 1
    bin_types = ['custom', 'lin']

    def test_bin_type(pair_counter, *args, **kwargs):
        res = (pair_counter(*args, bin_type=bin_type,
               weight_type="pair_product", **kwargs) for bin_type in bin_types)
        assert allclose(*res)


    print("Running 3-D correlation function DD(r) (isa={})".format(isa))
    test_bin_type(DD, autocorr, nthreads, binfile, x, y, z, weights1=w,
                  verbose=True, periodic=periodic, boxsize=boxsize, isa=isa)
    print("Running 3-D correlation function DD(r) (isa={})...PASSED".format(isa))


    print("\nRunning 2-D correlation function DD(rp,pi) (isa={})".format(isa))
    test_bin_type(DDrppi, autocorr, nthreads, pimax,
                  binfile, x, y, z, weights1=w,
                  verbose=True, periodic=periodic,
                  boxsize=boxsize, isa=isa)
    print("\nRunning 2-D correlation function DD(rp,pi) (isa={})...PASSED".format(isa))

    mu_max = 0.5
    nmu_bins = 10

    print("\nRunning 2-D correlation function DD(s,mu)")
    test_bin_type(DDsmu, autocorr, nthreads, binfile,
                  mu_max, nmu_bins,
                  x, y, z, weights1=w,
                  verbose=True, periodic=periodic,
                  boxsize=boxsize, output_savg=True, isa=isa)

    print("\nRunning 2-D projected correlation function wp(rp)")
    test_bin_type(wp, boxsize, pimax, nthreads,
                  binfile, x, y, z, weights=w,
                  verbose=True, isa=isa)

    print("\nRunning 3-D auto-correlation function xi(r)")
    test_bin_type(xi, boxsize, nthreads, binfile,
                  x, y, z, weights=w,
                  verbose=True, isa=isa)


if __name__ == '__main__':

    tests()
    for isa in ['fallback','sse42','avx','avx512f']:
        test_linear_binning_theory(isa=isa)
        test_linear_binning_mocks(isa=isa)
