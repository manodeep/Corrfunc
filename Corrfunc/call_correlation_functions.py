#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example python code to call the 3 correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os.path as path
import time
from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog


def main():
    tstart = time.time()
    t0 = tstart
    x, y, z = read_catalog()
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")
    boxsize = 420.0
    nthreads = 4
    pimax = 40.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../xi_theory/tests/", "bins")
    autocorr = 1
    numbins_to_print = 5

    print("Running 3-D correlation function DD(r)")
    results_DD = _countpairs.countpairs(autocorr, nthreads, binfile, x, y, z,
                                        x, y, z)
    print("\n#      **** DD(r): first {0} bins  *******       "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg       npairs")
    print("################################################")
    for ibin in range(numbins_to_print):
        items = results_DD[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d}"
              .format(items[0], items[1], items[2], items[3]))
    print("------------------------------------------------")

    print("\nRunning 2-D correlation function DD(rp,pi)")
    results_DDrppi = _countpairs.countpairs_rp_pi(autocorr, nthreads, pimax,
                                                  binfile, x, y, z, x, y, z)
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs")
    print("###########################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-----------------------------------------------------------")

    print("\nRunning 2-D projected correlation function wp(rp)")
    results_wp = _countpairs.countpairs_wp(boxsize, pimax, nthreads,
                                           binfile, x, y, z)
    print("\n#            ******    wp: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        wp       npairs")
    print("##########################################################")
    for ibin in range(numbins_to_print):
        items = results_wp[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-----------------------------------------------------------")

    print("\nRunning 3-D auto-correlation function xi(r)")
    results_xi = _countpairs.countpairs_xi(boxsize, nthreads, binfile, x, y, z)
    print("\n#            ******    xi: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        xi       npairs")
    print("##########################################################")
    for ibin in range(numbins_to_print):
        items = results_xi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-----------------------------------------------------------")
    print("Done with all four correlation calculations.")

    print("\nRunning VPF pN(r)")
    rmax = 10.0
    nbin = 10
    nspheres = 10000
    num_pN = 3
    seed = -1
    results_vpf = _countpairs.countspheres_vpf(rmax, nbin, nspheres, num_pN,
                                               seed, x, y, z)

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
