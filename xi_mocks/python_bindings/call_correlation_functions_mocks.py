"""
Example python code to call the 2 mocks correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import print_function
import os
import sys
import numpy as np
import _countpairs_mocks

file="../tests/data/Mr19_mock_northonly.rdcz.dat"
### Make sure the precision agrees with the definition in ../common.mk.
### Otherwise, you will get a runtime error -- 
### TypeError TypeError: array cannot be safely cast to required type
dtype=np.float

### Check if pandas is available - much faster to read in the data through pandas
try:
    import pandas as pd
    df = pd.read_csv(file,header=None,engine="c",dtype={"x":dtype,"y":dtype,"z":dtype},delim_whitespace=True)
    ra  = np.asarray(df[0],dtype=dtype)
    dec = np.asarray(df[1],dtype=dtype)
    cz  = np.asarray(df[2],dtype=dtype)
except:
    ra,dec,cz = np.genfromtxt(file,dtype=dtype,unpack=True)

nthreads=4
pimax=40.0
binfile="../tests/bins"
autocorr=1
numbins_to_print=5
cosmology=1

print("RA min = {} max = {}".format(np.min(ra),np.max(ra)))
print("DEC min = {} max = {}".format(np.min(dec),np.max(dec)))
print("cz min = {} max = {}".format(np.min(cz),np.max(cz)))

print("\nRunning 2-D correlation function xi(rp,pi)")
results_DDrppi = _countpairs_mocks.countpairs_rp_pi_mocks(autocorr, cosmology,nthreads,pimax,binfile,ra,dec,cz,ra,dec,cz)
print("\n#            ****** DD(rp,pi): first {} bins  *******      ".format(numbins_to_print))
print("#      rmin        rmax       rpavg     pi_upper     npairs")
print("###########################################################")
for ibin in xrange(numbins_to_print):
    items = results_DDrppi[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")


print("\nRunning angular correlation function w(theta)")
results_wp = _countpairs_mocks.countpairs_theta_mocks(autocorr, cosmology, nthreads, binfile, ra,dec, ra,dec)
print("\n#         ******  wtheta: first {} bins  *******        ".format(numbins_to_print))
print("#      thetamin        thetamax       thetaavg      npairs")
print("##########################################################")
for ibin in xrange(numbins_to_print):
    items = results_wp[ibin]
    print("{0:14.4f} {1:14.4f} {2:14.4f} {3:14d}".format(items[0],items[1],items[2],items[3]))
print("-----------------------------------------------------------")

print("Done with all MOCK correlation calculations.")





