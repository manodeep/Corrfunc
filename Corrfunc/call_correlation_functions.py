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


def main():
    from os.path import dirname, abspath, join as pjoin
    import time
    import numpy as np
    import Corrfunc
    from Corrfunc.io import read_fastfood_catalog
    from Corrfunc._countpairs import\
        countpairs as DD_extn,\
        countpairs_rp_pi as DDrppi_extn,\
        countpairs_wp as wp_extn,\
        countpairs_xi as xi_extn,\
        countspheres_vpf as vpf_extn,\
        countpairs_s_mu as DDsmu_extn

    def check_against_reference(results, filename, results_cols=(-2, -1, 2), ref_cols=(0, 4, 1)):
        # results is output of Python function
        # filename the reference counts
        # ref_cols in order of npairs, weightavg, rpavg
        refs = np.loadtxt(filename, unpack=True, usecols=ref_cols)
        for ii, items in enumerate(results):
            for icol, ref in zip(results_cols, refs):
                assert np.allclose(items[icol], ref[ii])

    tstart = time.time()
    t0 = tstart
    #x, y, z = read_catalog()
    #w = np.ones((1,len(x)), dtype=x.dtype)
    filename = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/data", "gals_Mr19.ff")
    x, y, z, w = read_fastfood_catalog(filename, need_weights=True)
    w = w[None,:] # shape (n_weights_per_particle,n_particles)
    boxsize = 420.0
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))

    numbins_to_print = 5

    print("Beginning Theory Correlation functions calculations")
    nthreads = 4
    pimax = 40.0
    binfile = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "bins")
    autocorr = 1
    periodic = 1

    print("Running 3-D correlation function DD(r)")
    results_DD, _ = DD_extn(autocorr, nthreads, binfile, x, y, z,
                            weights1=w, weight_type='pair_product',
                            periodic=periodic, boxsize=boxsize,
                            output_ravg=True, verbose=True)
    print("\n#      **** DD(r): first {0} bins  *******       "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg       npairs    weightavg")
    print("#############################################################")
    for ibin in range(numbins_to_print):
        items = results_DD[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d} {4:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "Mr19_DD_periodic")
    check_against_reference(results_DD, file_ref, results_cols=(3, 4, 2), ref_cols=(0, 4, 1))

    print("\nRunning 2-D correlation function DD(rp,pi)")
    results_DDrppi, _ = DDrppi_extn(autocorr, nthreads, pimax,
                                    binfile, x, y, z,
                                    weights1=w, weight_type='pair_product',
                                    periodic=periodic, boxsize=boxsize,
                                    output_rpavg=True, verbose=True)
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs    weightavg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "Mr19_DDrppi_periodic")
    check_against_reference(results_DDrppi, file_ref, results_cols=(4, 5, 2), ref_cols=(0, 4, 1))

    mu_max = 0.5
    nmu_bins = 10

    print("\nRunning 2-D correlation function DD(s,mu)")
    results_DDsmu, _ = DDsmu_extn(autocorr, nthreads, binfile,
                                    mu_max, nmu_bins,
                                    x, y, z,
                                    weights1=w, weight_type='pair_product',
                                    periodic=periodic, boxsize=boxsize,
                                    output_savg=True, verbose=True)
    print("\n#            ****** DD(s,mu): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      smin        smax       savg     mu_max     npairs    weightavg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDsmu[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "Mr19_DDsmu_periodic")
    check_against_reference(results_DDsmu, file_ref, results_cols=(4, 5, 2), ref_cols=(0, 4, 1))

    print("\nRunning 2-D projected correlation function wp(rp)")
    results_wp, _, _ = wp_extn(boxsize, pimax, nthreads,
                            binfile, x, y, z,
                            weights=w, weight_type='pair_product',
                            output_rpavg=True, verbose=True)
    print("\n#            ******    wp: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        wp       npairs    weightavg")
    print("#######################################################################")
    for ibin in range(numbins_to_print):
        items = results_wp[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("-----------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "Mr19_wp")
    check_against_reference(results_wp, file_ref, results_cols=(4, 5, 2, 3), ref_cols=(4, 5, 1, 0))

    print("\nRunning 3-D auto-correlation function xi(r)")
    results_xi, _ = xi_extn(boxsize, nthreads, binfile,
                            x, y, z,
                            weights=w, weight_type='pair_product',
                            output_ravg=True, verbose=True)

    print("\n#            ******    xi: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        xi       npairs    weightavg")
    print("#######################################################################")
    for ibin in range(numbins_to_print):
        items = results_xi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("-----------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "Mr19_xi")
    check_against_reference(results_xi, file_ref, results_cols=(4, 5, 2, 3), ref_cols=(4, 5, 1, 0))
    print("Done with all four correlation calculations.")

    print("\nRunning VPF pN(r)")
    rmax = 10.0
    nbin = 10
    nspheres = 10000
    num_pN = 6
    seed = -1234
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
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../theory/tests/", "Mr19_vpf_periodic")
    check_against_reference(results_vpf, file_ref, ref_cols=range(7), results_cols=range(7))

    tend = time.time()
    print("Done with all functions. Total time taken = {0:10.1f} seconds. \
    Read-in time = {1:10.1f} seconds.".format(tend - tstart, t1 - t0))


if __name__ == "__main__":
    main()
