#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example python code to call the theory clustering functions
from python. This script calls the python extensions directly;
however the recommended use is via the wrappers provided
in :py:mod:`Corrfunc.theory`.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from os.path import dirname, abspath, join as pjoin
import time
import numpy as np

import Corrfunc
from Corrfunc.io import read_catalog
from Corrfunc._countpairs import\
    countpairs as DD_extn,\
    countpairs_rp_pi as DDrppi_extn,\
    countpairs_wp as wp_extn,\
    countpairs_xi as xi_extn,\
    countspheres_vpf as vpf_extn,\
    countpairs_s_mu as DDsmu_extn


def main():
    tstart = time.time()
    t0 = tstart
    x, y, z = read_catalog()
    w = np.ones((1,len(x)), dtype=x.dtype)
    boxsize = 420.0
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))

    numbins_to_print = 5

    print("Beginning Theory Correlation functions calculations")
    nthreads = 4
    pimax = 40.0
    binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                    "../theory/tests/", "bins")
    autocorr = 1
    periodic = 1

    print("Running 3-D correlation function DD(r)")
    results_DD, _ = DD_extn(autocorr, nthreads, binfile, x, y, z,
                            weights1=w, weight_type='pair_product',
                            verbose=True, periodic=periodic, boxsize=boxsize)
    print("\n#      **** DD(r): first {0} bins  *******       "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg       npairs    weightavg")
    print("#############################################################")
    for ibin in range(numbins_to_print):
        items = results_DD[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d} {4:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-------------------------------------------------------------")

    print("\nRunning 2-D correlation function DD(rp,pi)")
    results_DDrppi, _ = DDrppi_extn(autocorr, nthreads, pimax,
                                    binfile, x, y, z,
                                    weights1=w, weight_type='pair_product',
                                    verbose=True, periodic=periodic,
                                    boxsize=boxsize)
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs    weightavg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")

    mu_max = 0.5
    nmu_bins = 10

    print("\nRunning 2-D correlation function DD(s,mu)")
    results_DDsmu, _ = DDsmu_extn(autocorr, nthreads, binfile,
                                    mu_max, nmu_bins,
                                    x, y, z,
                                    weights1=w, weight_type='pair_product',
                                    verbose=True, periodic=periodic,
                                    boxsize=boxsize, output_savg=True)
    print("\n#            ****** DD(s,mu): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      smin        smax       savg     mu_max     npairs    weightavg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDsmu[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")

    print("\nRunning 2-D projected correlation function wp(rp)")
    results_wp, _, _ = wp_extn(boxsize, pimax, nthreads,
                            binfile, x, y, z,
                            weights=w, weight_type='pair_product',
                            verbose=True)
    print("\n#            ******    wp: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        wp       npairs    weightavg")
    print("#######################################################################")
    for ibin in range(numbins_to_print):
        items = results_wp[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("-----------------------------------------------------------------------")

    print("\nRunning 3-D auto-correlation function xi(r)")
    results_xi, _ = xi_extn(boxsize, nthreads, binfile,
                            x, y, z,
                            weights=w, weight_type='pair_product',
                            verbose=True)

    print("\n#            ******    xi: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        xi       npairs    weightavg")
    print("#######################################################################")
    for ibin in range(numbins_to_print):
        items = results_xi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("-----------------------------------------------------------------------")
    print("Done with all four correlation calculations.")

    print("\nRunning VPF pN(r)")
    rmax = 10.0
    nbin = 10
    nspheres = 10000
    num_pN = 3
    seed = -1
    results_vpf, _ = vpf_extn(rmax, nbin, nspheres, num_pN,
                              seed, x, y, z, verbose=True, periodic=periodic,
                              boxsize=boxsize)

    print("\n#            ******    pN: first {0} bins  *******         "
          .format(numbins_to_print))
    print('#       r    ', end="")

    for ipn in range(num_pN):
        print('        p{0:0d}      '.format(ipn), end="")

    print("")

    print("###########", end="")
    for ipn in range(num_pN):
        print('################', end="")
    print("")

    for ibin in range(numbins_to_print):
        items = results_vpf[ibin]
        print('{0:10.2f} '.format(items[0]), end="")
        for ipn in range(num_pN):
            print(' {0:15.4e}'.format(items[ipn + 1]), end="")
        print("")

    print("-----------------------------------------------------------")

    tend = time.time()
    print("Done with all functions. Total time taken = {0:10.1f} seconds. \
    Read-in time = {1:10.1f} seconds.".format(tend - tstart, t1 - t0))

if __name__ == "__main__":
    main()
