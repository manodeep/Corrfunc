#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example python code to call the mocks clustering functions
from python. This script calls the python extensions
directly; however the recommended use is via the wrappers provided
in :py:mod:`Corrfunc.mocks`.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


def main():
    from os.path import dirname, abspath, join as pjoin
    import numpy as np
    import time
    import Corrfunc
    from Corrfunc.io import read_fastfood_catalog
    from Corrfunc._countpairs_mocks import\
        countpairs_rp_pi_mocks as rp_pi_mocks_extn,\
        countpairs_s_mu_mocks as s_mu_mocks_extn,\
        countpairs_theta_mocks as theta_mocks_extn,\
        countspheres_vpf_mocks as vpf_mocks_extn

    def check_against_reference(results, filename, results_cols=(-2, -1, 2), ref_cols=(0, 4, 1)):
        # results is output of Python function
        # filename the reference counts
        # ref_cols in order of npairs, weightavg, rpavg
        refs = np.loadtxt(filename, unpack=True, usecols=ref_cols)
        for ii, items in enumerate(results):
            for icol, ref in zip(results_cols, refs):
                assert np.allclose(items[icol], ref[ii])

    tstart = time.time()
    t0 = time.time()
    filename = pjoin(dirname(abspath(Corrfunc.__file__)),
                     "../mocks/tests/data/", "Mr19_mock_northonly.rdcz.ff")
    ra, dec, cz, w = read_fastfood_catalog(filename, need_weights=True)
    w = w[None,:]

    filename = pjoin(dirname(abspath(Corrfunc.__file__)),
                     "../mocks/tests/data/", "Mr19_randoms_northonly.rdcz.ff")
    ra_ran, dec_ran, cz_ran, w_ran = read_fastfood_catalog(filename, need_weights=True)
    w_ran = w_ran[None,:]

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
                    "../mocks/tests/", "bins")
    autocorr = 1
    numbins_to_print = 5
    cosmology = 1

    print("\nRunning 2-D correlation function xi(rp,pi)")
    results_DDrppi, _ = rp_pi_mocks_extn(autocorr, cosmology, nthreads,
                                         pimax, binfile,
                                         ra, dec, cz,
                                         weights1=w, weight_type='pair_product',
                                         output_rpavg=True, verbose=True)
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs   weightavg")
    print("########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:10.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("------------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../mocks/tests/", "Mr19_mock.DD")
    check_against_reference(results_DDrppi, file_ref, results_cols=(4, 5, 2), ref_cols=(0, 4, 1))

    nmu_bins = 10
    mu_max = 1.0

    print("\nRunning 2-D correlation function xi(s,mu)")
    results_DDsmu, _ = s_mu_mocks_extn(0, cosmology, nthreads,
                                       mu_max, nmu_bins, binfile,
                                       ra, dec, cz, weights1=w,
                                       RA2=ra_ran, DEC2=dec_ran, CZ2=cz_ran, weights2=w_ran,
                                       output_savg=True, verbose=True,
                                       weight_type='pair_product')
    print("\n#            ****** DD(s,mu): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      smin        smax       savg     mu_upper       npairs    weight_avg")
    print("###########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDsmu[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:12.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))
    print("--------------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../mocks/tests/", "Mr19_mock_DDsmu.DR")
    check_against_reference(results_DDsmu, file_ref, results_cols=(4, 5, 2), ref_cols=(0, 4, 1))

    binfile = pjoin(dirname(abspath(__file__)),
                    "../mocks/tests/", "angular_bins")
    print("\nRunning angular correlation function DD(theta)")
    results_wtheta, _ = theta_mocks_extn(autocorr, nthreads, binfile,
                                         ra, dec, RA2=ra, DEC2=dec,
                                         weights1=w,
                                         weights2=w,
                                         weight_type='pair_product',
                                         output_thetaavg=True, fast_acos=True,
                                         verbose=1)
    print("\n#         ******  DD(theta): first {0} bins  *******        "
          .format(numbins_to_print))
    print("#      thetamin        thetamax       thetaavg        npairs      weightavg")
    print("############################################################################")
    for ibin in range(numbins_to_print):
        items = results_wtheta[ibin]
        print("{0:14.4f} {1:14.4f} {2:14.4f} {3:14d} {4:14.4f}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-----------------------------------------------------------------------")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../mocks/tests/", "Mr19_mock_wtheta.DD")
    check_against_reference(results_wtheta, file_ref, results_cols=(3, 4, 2), ref_cols=(0, 4, 1))

    print("Beginning the VPF")
    # Max. sphere radius of 10 Mpc
    rmax = 10.0
    # 10 bins..so counts in spheres of radius 1, 2, 3, 4...10 Mpc spheres
    nbin = 10
    num_spheres = 10000
    num_pN = 6
    threshold_neighbors = 1  # does not matter since we have the centers
    centers_file = pjoin(dirname(abspath(__file__)),
                         "../mocks/tests/data/",
                         "Mr19_centers_xyz_forVPF_rmax_10Mpc.txt")
    results_vpf, _ = vpf_mocks_extn(rmax, nbin, num_spheres, num_pN,
                                    threshold_neighbors, centers_file,
                                    cosmology,
                                    ra, dec, cz, ra, dec, cz, verbose=True)

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
    print("Done with the VPF.")
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../mocks/tests/", "Mr19_mock_vpf")
    check_against_reference(results_vpf, file_ref, ref_cols=range(7), results_cols=range(7))

    tend = time.time()
    print("Done with all the MOCK clustering calculations. Total time \
    taken = {0:0.2f} seconds.".format(tend - tstart))


if __name__ == "__main__":
    main()
