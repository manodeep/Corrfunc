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
binfile="../tests/bins"
dtype='float' ### Make sure the precision agrees with the definition in ../common.mk

x,y,z = np.loadtxt(file,dtype=dtype,unpack=True)

array = np.hstack((x,y,z))
print("max = {} min = {} ".format(np.max(array),np.min(array)))
boxsize=np.max(array)
nthreads=4
pimax=40.0
results = _countpairs.countpairs_wp(boxsize,pimax,nthreads,binfile,x,y,z)





