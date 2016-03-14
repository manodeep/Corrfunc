"""
Example python code to call the 2 mocks correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import print_function
import os.path as path
import sys
import re
import numpy as np
import time


# Import from current directory first,
# and then from the package.
from _countpairs_mocks import countpairs_rp_pi_mocks as rp_pi_mocks
from _countpairs_mocks import countpairs_theta_mocks as theta_mocks
from _countpairs_mocks import countspheres_vpf_mocks as vpf_mocks

if sys.version_info[0] >= 3:
    def rd(filename):
        with open(filename, encoding="utf-8") as f:
            r = f.read()

        return r
else:
    def rd(filename):
        with open(filename) as f:
            r = f.read()

        return r

try:
    import pandas as pd
except ImportError:
    pd = None


def main():
    tstart = time.time()
    file = path.join(path.dirname(path.abspath(__file__)),
                     "../tests/data/", "Mr19_mock_northonly.rdcz.dat")
    # Figure out the datatype, use the header file in the include directory
    # because that is most likely correct (common.mk might have been modified
    # but not recompiled)
    include_file = path.join(path.dirname(path.abspath(__file__)),
                             "../../include/", "countpairs_rp_pi_mocks.h")
    includes = rd(include_file)
    vector_type = re.search(r'(\w+)\s*\*\s*rupp\s*\;', includes, re.I).group(1)
    allowed_types = {"float": np.float32, "double": np.float}
    if vector_type not in list(allowed_types.keys()):
        print("Error: Unknown precision={0} found in header file {1}. \
        Allowed types are `{2}'"
              .format(vector_type, include_file, allowed_types))
        sys.exit()

    dtype = allowed_types[vector_type]

    # Check if pandas is available - much faster to read in the
    # data through pandas
    t0 = time.time()
    print("Reading in the data...")
    if pd is not None:
        df = pd.read_csv(file, header=None, engine="c",
                         dtype={"x": dtype, "y": dtype, "z": dtype},
                         delim_whitespace=True)
        ra = np.asarray(df[0], dtype=dtype)
        dec = np.asarray(df[1], dtype=dtype)
        cz = np.asarray(df[2], dtype=dtype)
    else:
        ra, dec, cz = np.genfromtxt(file, dtype=dtype, unpack=True)

    t1 = time.time()
    print("RA min  = {0} max = {1}".format(np.min(ra), np.max(ra)))
    print("DEC min = {0} max = {1}".format(np.min(dec), np.max(dec)))
    print("cz min  = {0} max = {1}".format(np.min(cz), np.max(cz)))
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")

    nthreads = 4
    pimax = 40.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../tests/", "bins")
    autocorr = 1
    numbins_to_print = 5
    cosmology = 1

    print("\nRunning 2-D correlation function xi(rp,pi)")
    results_DDrppi = rp_pi_mocks(autocorr, cosmology, nthreads,
                                 pimax, binfile,
                                 ra, dec, cz, ra, dec, cz)
    print("\n#            ****** DD(rp,pi): first {0} bins  *******      "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg     pi_upper     npairs")
    print("###########################################################")
    for ibin in range(numbins_to_print):
        items = results_DDrppi[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}"
              .format(items[0], items[1], items[2], items[3], items[4]))

    print("-----------------------------------------------------------")

    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../tests/", "angular_bins")
    print("\nRunning angular correlation function w(theta)")
    results_wtheta = theta_mocks(autocorr, cosmology, nthreads, binfile,
                                 ra, dec, ra, dec)
    print("\n#         ******  wtheta: first {0} bins  *******        "
          .format(numbins_to_print))
    print("#      thetamin        thetamax       thetaavg      npairs")
    print("##########################################################")
    for ibin in range(numbins_to_print):
        items = results_wtheta[ibin]
        print("{0:14.4f} {1:14.4f} {2:14.4f} {3:14d}"
              .format(items[0], items[1], items[2], items[3]))
    print("-----------------------------------------------------------")

    print("Beginning the VPF")
    # Max. sphere radius of 10 Mpc
    rmax = 10.0
    # 10 bins..so counts in spheres of radius 1, 2, 3, 4...10 Mpc spheres
    nbin = 10
    num_spheres = 10000
    num_pN = 6
    threshold_neighbors = 1  # does not matter since we have the centers
    centers_file = path.join(path.dirname(path.abspath(__file__)),
                             "../tests/data/",
                             "Mr19_centers_xyz_forVPF_rmax_10Mpc.txt")
    results_vpf = vpf_mocks(rmax, nbin, num_spheres, num_pN,
                            threshold_neighbors, centers_file, cosmology,
                            ra, dec, cz, ra, dec, cz)
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
