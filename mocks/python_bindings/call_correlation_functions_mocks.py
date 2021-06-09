"""
Example python code to call the 2 mocks correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import print_function
from os.path import dirname, abspath, join as pjoin
import time
import numpy as np

from _countpairs_mocks import \
    countpairs_rp_pi_mocks as rp_pi_mocks,\
    countpairs_theta_mocks as theta_mocks,\
    countspheres_vpf_mocks as vpf_mocks, \
    countpairs_s_mu_mocks as s_mu_mocks


try:
    import pandas as pd
except ImportError:
    pd = None


def read_text_file(filename, encoding="utf-8"):
    """
    Reads a file under python3 with encoding (default UTF-8).
    Also works under python2, without encoding.
    Uses the EAFP (https://docs.python.org/2/glossary.html#term-eafp)
    principle.
    """
    try:
        with open(filename, 'r', encoding) as f:
            r = f.read()
    except TypeError:
        with open(filename, 'r') as f:
            r = f.read()
    return r


def main():
    tstart = time.time()
    filename = pjoin(dirname(abspath(__file__)),
                     "../tests/data/", "Mr19_mock_northonly.rdcz.dat")
    # Double-precision calculations
    # (if you want single-prec, just change the following line
    # to dtype = np.float32)
    dtype = np.float64

    # Check if pandas is available - much faster to read in the
    # data through pandas
    t0 = time.time()
    print("Reading in the data...")
    if pd is not None:
        df = pd.read_csv(filename, header=None, engine="c",
                         dtype={"x": dtype, "y": dtype, "z": dtype},
                         delim_whitespace=True)
        ra = np.asarray(df[0], dtype=dtype)
        dec = np.asarray(df[1], dtype=dtype)
        cz = np.asarray(df[2], dtype=dtype)
        weights = np.asarray(df[3], dtype=dtype)
    else:
        ra, dec, cz, weights = np.genfromtxt(filename, dtype=dtype,
                                             unpack=True)

    weights = weights.reshape(1,-1)

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
                    "../tests/", "bins")
    autocorr = 1
    numbins_to_print = 5
    cosmology = 1

    print("\nRunning 2-D correlation function xi(rp,pi)")
    results_DDrppi, _ = rp_pi_mocks(autocorr, cosmology, nthreads,
                                    pimax, binfile,
                                    ra, dec, cz, weights1=weights,
                                    output_rpavg=True, verbose=True,
                                    weight_type='pair_product')
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs     weight_avg")
    print("##########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:12.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))

    print("--------------------------------------------------------------------------")

    print("\nRunning 2-D correlation function xi(rp,pi) with different bin refinement")
    results_DDrppi, _ = rp_pi_mocks(autocorr, cosmology, nthreads,
                                    pimax, binfile,
                                    ra, dec, cz,
                                    output_rpavg=True,
                                    xbin_refine_factor=3,
                                    ybin_refine_factor=3,
                                    zbin_refine_factor=2,
                                    verbose=True)
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs")
    print("###########################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}"
              .format(items[0], items[1], items[2], items[3], items[4]))

    print("-----------------------------------------------------------")

    nmu_bins = 10
    mu_max = 1.0

    print("\nRunning 2-D correlation function xi(s,mu)")
    results_DDsmu, _ = s_mu_mocks(autocorr, cosmology, nthreads,
                                  mu_max, nmu_bins, binfile,
                                  ra, dec, cz, weights1=weights,
                                  output_savg=True, verbose=True,
                                  weight_type='pair_product')
    print("\n#            ****** DD(s,mu): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      smin        smax       savg     mu_upper    npairs     weight_avg")
    print("##########################################################################")
    for ibin in range(numbins_to_print):
        items = results_DDsmu[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d} {5:12.4f}"
              .format(items[0], items[1], items[2], items[3], items[4], items[5]))

    print("--------------------------------------------------------------------------")

    binfile = pjoin(dirname(abspath(__file__)),
                    "../tests/", "angular_bins")
    print("\nRunning angular correlation function w(theta)")
    results_wtheta, _ = theta_mocks(autocorr, nthreads, binfile,
                                    ra, dec, weights1=weights,
                                    RA2=ra, DEC2=dec, weights2=weights,
                                    output_thetaavg=True, fast_acos=True,
                                    verbose=1, weight_type='pair_product')
    print("\n#         ******  wtheta: first {0} bins  *******        "
          .format(numbins_to_print))
    print("#      thetamin        thetamax       thetaavg      npairs    weightavg")
    print("#######################################################################")
    for ibin in range(numbins_to_print):
        items = results_wtheta[ibin]
        print("{0:14.4f} {1:14.4f} {2:14.4f} {3:14d} {4:14.4f}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-----------------------------------------------------------------------")

    print("Beginning the VPF")
    # Max. sphere radius of 10 Mpc
    rmax = 10.0
    # 10 bins..so counts in spheres of radius 1, 2, 3, 4...10 Mpc spheres
    nbin = 10
    num_spheres = 10000
    num_pN = 6
    threshold_neighbors = 1  # does not matter since we have the centers
    centers_file = pjoin(dirname(abspath(__file__)),
                         "../tests/data/",
                         "Mr19_centers_xyz_forVPF_rmax_10Mpc.txt")
    results_vpf, _ = vpf_mocks(rmax, nbin, num_spheres, num_pN,
                               threshold_neighbors, centers_file, cosmology,
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
    tend = time.time()
    print("Done with all the MOCK clustering calculations. Total time \
    taken = {0:0.2f} seconds.".format(tend - tstart))

if __name__ == "__main__":
    main()
