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
import time

tstart=time.time()
file="../tests/data/gals_Mr19.txt"
### Make sure the precision agrees with the definition in ../common.mk.
### Otherwise, you will get a runtime error -- 
### TypeError TypeError: array cannot be safely cast to required type
dtype=np.float32

### Check if pandas is available - much faster to read in the data through pandas
t0=time.time()
print("Reading in the data...")
try:
    import pandas as pd
    df = pd.read_csv(file,header=None,engine="c",dtype={"x":dtype,"y":dtype,"z":dtype},delim_whitespace=True)
    x = np.asarray(df[0],dtype=dtype)
    y = np.asarray(df[1],dtype=dtype)
    z = np.asarray(df[2],dtype=dtype)
except:
    print("Warning: Could not read in data with pandas -- due to error : {}. Falling back to slower numpy.".format(sys.exc_info()[0]))
    x,y,z = np.genfromtxt(file,dtype=dtype,unpack=True)

t1=time.time()    
print("Done reading the data - time taken = {0:10.1f} seconds.\nBeginning Correlation functions calculations".format(t1-t0))

boxsize=420.0
nthreads=4
pimax=40.0
binfile="../tests/bins"
autocorr=1
numbins_to_print=5

print("Running 3-D correlation function DD(r)")
results_DD = _countpairs.countpairs(autocorr,nthreads,binfile,x,y,z,x,y,z)
print("\n#      **** DD(r): first {} bins  *******       ".format(numbins_to_print))
print("#      rmin        rmax       rpavg       npairs")
print("################################################")
for ibin in range(numbins_to_print):
    items = results_DD[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10d}".format(items[0],items[1],items[2],items[3]))
print("------------------------------------------------")


print("\nRunning 2-D correlation function DD(rp,pi)")
results_DDrppi = _countpairs.countpairs_rp_pi(autocorr,nthreads,pimax,binfile,x,y,z,x,y,z)
print("\n#            ****** DD(rp,pi): first {} bins  *******      ".format(numbins_to_print))
print("#      rmin        rmax       rpavg     pi_upper     npairs")
print("###########################################################")
for ibin in range(numbins_to_print):
    items = results_DDrppi[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")


print("\nRunning 2-D projected correlation function wp(rp)")
results_wp = _countpairs.countpairs_wp(boxsize,pimax,nthreads,binfile,x,y,z)
print("\n#            ******    wp: first {} bins  *******         ".format(numbins_to_print))
print("#      rmin        rmax       rpavg        wp       npairs")
print("##########################################################")
for ibin in range(numbins_to_print):
    items = results_wp[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")

print("\nRunning 3-D auto-correlation function xi(r)")
results_xi = _countpairs.countpairs_xi(boxsize,nthreads,binfile,x,y,z)
print("\n#            ******    xi: first {} bins  *******         ".format(numbins_to_print))
print("#      rmin        rmax       rpavg        xi       npairs")
print("##########################################################")
for ibin in range(numbins_to_print):
    items = results_xi[ibin]
    print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}".format(items[0],items[1],items[2],items[3],items[4]))
print("-----------------------------------------------------------")

print("Done with all four correlation calculations.")

print("\nRunning VPF pN(r)")
rmax=10.0
nbin=10
nspheres=10000
num_pN=3
seed=-1
results_vpf = _countpairs.countspheres_vpf(rmax,nbin,nspheres,num_pN,seed,x,y,z)

print("\n#            ******    pN: first {} bins  *******         ".format(numbins_to_print))
print('#       r    ',end="")

for ipn in range(num_pN):
    print('        p{:0d}      '.format(ipn),end="")

print("")

print("###########",end="")
for ipn in range(num_pN):
    print('################',end="")
print("")


for ibin in range(numbins_to_print):
    items = results_vpf[ibin]
    print('{0:10.2f} '.format(items[0]),end="")
    for ipn in range(num_pN):
        print(' {0:15.4e}'.format(items[ipn+1]),end="")
    print("")

print("-----------------------------------------------------------")

tend=time.time()
print("Done with all functions. Total time taken = {0:10.1f} seconds. Read-in time = {1:10.1f} seconds.".format(tend-tstart,t1-t0))


