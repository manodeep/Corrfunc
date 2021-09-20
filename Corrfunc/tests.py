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


def test_linear_binning_mocks(isa='fallback'):
    """Here we test that the special treatment for linear binning returns the correct result."""
    from os.path import dirname, abspath, join as pjoin
    import numpy as np
    import time
    import Corrfunc
    from Corrfunc.io import read_catalog
    from Corrfunc.mocks import DDrppi_mocks, DDsmu_mocks, DDtheta_mocks
    tstart = time.time()
    filename = pjoin(dirname(abspath(Corrfunc.__file__)),
                     "../mocks/tests/data/", "Mr19_mock_northonly.rdcz.ff")

    t0 = time.time()
    ra, dec, cz = read_catalog(filename)
    w = np.ones((1, len(ra)), dtype=ra.dtype)
    t1 = time.time()
    print("RA min  = {0} max = {1}".format(np.min(ra), np.max(ra)))
    print("DEC min = {0} max = {1}".format(np.min(dec), np.max(dec)))
    print("cz min  = {0} max = {1}".format(np.min(cz), np.max(cz)))
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")

    nthreads = 4
    pimax = 40.0
    binfile = pjoin(dirname(abspath(__file__)),
                    "../mocks/tests/", "bins_lin")
    rbins = np.linspace(0.1, 45.1, 46)
    with open(binfile,'w') as file:
        for low,hi in zip(rbins[:-1], rbins[1:]):
            file.write("{0} {1}\n".format(low, hi))
    autocorr = 1
    numbins_to_print = 5
    cosmology = 1

    def allclose(a, b):
        return all(np.allclose(a[name], b[name]) for name in a.dtype.names)

    print("\nRunning 2-D correlation function xi(rp,pi)")
    results_DDrppi_ref = DDrppi_mocks(autocorr, cosmology, nthreads,
                                      pimax, binfile,
                                      ra, dec, cz,
                                      weights1=w, weight_type='pair_product',
                                      output_rpavg=True, verbose=True, isa=isa, bin_type='custom')
    results_DDrppi = DDrppi_mocks(autocorr, cosmology, nthreads,
                                  pimax, binfile,
                                  ra, dec, cz,
                                  weights1=w, weight_type='pair_product',
                                  output_rpavg=True, verbose=False, isa=isa, bin_type='lin')
    assert allclose(results_DDrppi, results_DDrppi_ref)

    nmu_bins = 10
    mu_max = 1.0

    print("\nRunning 2-D correlation function xi(s,mu)")
    results_DDsmu_ref = DDsmu_mocks(autocorr, cosmology, nthreads,
                                    mu_max, nmu_bins, binfile,
                                    ra, dec, cz, weights1=w,
                                    output_savg=True, verbose=True,
                                    weight_type='pair_product', isa=isa, bin_type='custom')
    results_DDsmu = DDsmu_mocks(autocorr, cosmology, nthreads,
                                mu_max, nmu_bins, binfile,
                                ra, dec, cz, weights1=w,
                                output_savg=True, verbose=False,
                                weight_type='pair_product', isa=isa, bin_type='lin')
    assert allclose(results_DDsmu, results_DDsmu_ref)


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
    x, y, z = read_catalog()
    w = np.ones((1, len(x)), dtype=x.dtype)
    boxsize = 420.0
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))

    numbins_to_print = 5

    print("Beginning Theory Correlation functions calculations")
    nthreads = 4
    pimax = 40.0
    binfile = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "bins_lin")
    rbins = np.linspace(0.1, 25.1, 26)
    with open(binfile,'w') as file:
        for low, hi in zip(rbins[:-1], rbins[1:]):
            file.write("{0} {1}\n".format(low, hi))
    autocorr = 1
    periodic = 1

    def allclose(a, b):
        return all(np.allclose(a[name], b[name]) for name in a.dtype.names)

    print("Running 3-D correlation function DD(r)")
    results_DD_ref = DD(autocorr, nthreads, binfile, x, y, z,
                        weights1=w, weight_type='pair_product',
                        verbose=True, periodic=periodic, boxsize=boxsize,
                        isa=isa, bin_type='custom')
    results_DD = DD(autocorr, nthreads, binfile, x, y, z,
                    weights1=w, weight_type='pair_product',
                    verbose=False, periodic=periodic, boxsize=boxsize,
                    isa=isa, bin_type='auto')
    allclose(results_DD, results_DD_ref)

    print("\nRunning 2-D correlation function DD(rp,pi)")
    results_DDrppi_ref = DDrppi(autocorr, nthreads, pimax,
                                binfile, x, y, z,
                                weights1=w, weight_type='pair_product',
                                verbose=False, periodic=periodic,
                                boxsize=boxsize, isa=isa, bin_type='custom')
    results_DDrppi = DDrppi(autocorr, nthreads, pimax,
                            binfile, x, y, z,
                            weights1=w, weight_type='pair_product',
                            verbose=False, periodic=periodic,
                            boxsize=boxsize, isa=isa, bin_type='lin')
    allclose(results_DDrppi, results_DDrppi_ref)

    mu_max = 0.5
    nmu_bins = 10

    print("\nRunning 2-D correlation function DD(s,mu)")
    results_DDsmu_ref = DDsmu(autocorr, nthreads, binfile,
                              mu_max, nmu_bins,
                              x, y, z,
                              weights1=w, weight_type='pair_product',
                              verbose=True, periodic=periodic,
                              boxsize=boxsize, output_savg=True,
                              isa=isa, bin_type='custom')
    results_DDsmu = DDsmu(autocorr, nthreads, binfile,
                          mu_max, nmu_bins,
                          x, y, z,
                          weights1=w, weight_type='pair_product',
                          verbose=False, periodic=periodic,
                          boxsize=boxsize, output_savg=True,
                          isa=isa, bin_type='lin')
    allclose(results_DDsmu, results_DDsmu_ref)

    print("\nRunning 2-D projected correlation function wp(rp)")
    results_wp_ref = wp(boxsize, pimax, nthreads,
                        binfile, x, y, z,
                        weights=w, weight_type='pair_product',
                        verbose=True, isa=isa, bin_type='custom')
    results_wp = wp(boxsize, pimax, nthreads,
                    binfile, x, y, z,
                    weights=w, weight_type='pair_product',
                    verbose=False, isa=isa, bin_type='custom')
    allclose(results_wp, results_wp_ref)

    print("\nRunning 3-D auto-correlation function xi(r)")
    results_xi_ref = xi(boxsize, nthreads, binfile,
                        x, y, z,
                        weights=w, weight_type='pair_product',
                        verbose=True, isa=isa, bin_type='custom')
    results_xi = xi(boxsize, nthreads, binfile,
                    x, y, z,
                    weights=w, weight_type='pair_product',
                    verbose=False, isa=isa, bin_type='lin')
    allclose(results_xi, results_xi_ref)


if __name__ == '__main__':

    tests()
    for isa in ['fallback','sse42','avx','avx512f']:
        test_linear_binning_mocks(isa=isa)
        test_linear_binning_theory(isa=isa)
