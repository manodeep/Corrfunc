"""
Example python code to call the 3 correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

Requires: numpy

"""
from __future__ import print_function
import os
import sys
import numpy as np
import _countpairs

# try:
#     import pandas as pd
# except ImportError:
#     pd = None

file="../tests/data/gals_Mr19.txt"
### Make sure the precision agrees with the definition in ../common.mk.
### Otherwise, you will get a runtime error -- 
### TypeError TypeError: array cannot be safely cast to required type
dtype=np.float32

### Check if pandas is available - much faster to read in the data through pandas
try:
    import pandas as pd
    df = pd.read_csv(file,header=None,engine="c",dtype={"x":dtype,"y":dtype,"z":dtype},delim_whitespace=True)
    x = np.asarray(df[0],dtype=dtype)
    y = np.asarray(df[1],dtype=dtype)
    z = np.asarray(df[2],dtype=dtype)
except:
    x,y,z = np.genfromtxt(file,dtype=dtype,unpack=True)

boxsize=420.0
nthreads=4
pimax=40.0
binfile="../tests/bins"
autocorr=1
numbins_to_print=5

print("Running 3-D correlation function xi(r)")
results_DD = _countpairs.countpairs(autocorr,nthreads,binfile,x,y,z,x,y,z)
print("\n#      **** xi(r): first {} bins  *******       ".format(numbins_to_print))
print("#      rmin        rmax       rpavg       npairs")
print("################################################")
for ibin in xrange(numbins_to_print):
    items = results_DD[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d}".format(items[0],items[1],items[2],items[3]))
print("------------------------------------------------")


print("\nRunning 2-D correlation function xi(rp,pi)")
results_DDrppi = _countpairs.countpairs_rp_pi(autocorr,nthreads,pimax,binfile,x,y,z,x,y,z)
print("\n#            ****** xi(rp,pi): first {} bins  *******      ".format(numbins_to_print))
print("#      rmin        rmax       rpavg     pi_upper     npairs")
print("###########################################################")
for ibin in xrange(numbins_to_print):
    items = results_DDrppi[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")


print("\nRunning 2-D projected correlation function wp(rp)")
results_wp = _countpairs.countpairs_wp(boxsize,pimax,nthreads,binfile,x,y,z)
print("\n#            ******    wp: first {} bins  *******         ".format(numbins_to_print))
print("#      rmin        rmax       rpavg        wp       npairs")
print("##########################################################")
for ibin in xrange(numbins_to_print):
    items = results_wp[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")

print("Done with all three correlation calculations.")





