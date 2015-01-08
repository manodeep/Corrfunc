"""
Example python code to call the 3 correlation function
routines from python. (The codes are written in C)

Author: Manodeep Sinha <manodeep@gmail.com>

"""
from __future__ import print_function
import os
import sys
import numpy as np
import _countpairs


file="../tests/data/gals_Mr19.txt"
### Make sure the precision agrees with the definition in ../common.mk.
### Otherwise, you will get TypeError TypeError: array cannot be safely cast to required type
dtype='float32' 
x,y,z = np.loadtxt(file,dtype=dtype,unpack=True)
array = np.hstack((x,y,z))
boxsize=np.max(array)
nthreads=4
pimax=40.0
binfile="../tests/bins"
autocorr=1
results_wp = _countpairs.countpairs_wp(boxsize,pimax,nthreads,binfile,x,y,z)
print("{}".format(results_wp))
results_DDrppi = _countpairs.countpairs_rp_pi(autocorr,nthreads,pimax,binfile,x,y,z,x,y,z)
print("{}".format(results_DDrppi))
results_DD = _countpairs.countpairs(autocorr,nthreads,binfile,x,y,z,x,y,z)
print("{}".format(results_DD))

        




