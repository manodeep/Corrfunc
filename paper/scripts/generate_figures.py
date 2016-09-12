from __future__ import print_function
import numpy as np

import Corrfunc
from Corrfunc._countpairs import countpairs as DD, \
    countpairs_rp_pi as DDrppi, \
    countpairs_wp as wp, \
    countpairs_xi as xi

from Corrfunc.io import read_catalog
from os.path import join as pjoin, abspath, dirname
import time


def benchmark_theory_threads_all(min_threads=1, max_threads=16,
                                 keys = ['DDrppi', 'DD', 'wp', 'xi']):

    allkeys =  ['DDrppi', 'DD', 'wp', 'xi']
    if keys is None:
        keys = allkeys
    else:
        for k in keys:
            if k not in allkeys:
                msg = "Valid routines to benchmark are: {0}\nFound routine"\
                    " = {1}".format(allkeys, k)
                raise ValueError(msg)

    print("Benchmarking routines = {0}".format(keys))
    x, y, z = read_catalog()
    binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                    "../xi_theory/tests/", "bins")
    autocorr = 1
    pimax = 40.0
    boxsize = 420.0
    dtype = np.dtype([('name', 'S16'),
                      ('nthreads', np.int),
                      ('runtime', np.float)])

    totN = (max_threads - min_threads + 1) * len(keys)
    all_runtimes = np.empty(totN, dtype=dtype)
    index = 0L
    for nthreads in xrange(min_threads, max_threads + 1):

        if 'DDrppi' in keys:
            t0 = time.time()
            _ = DDrppi(autocorr, nthreads, pimax, binfile, x, y, z, x, y, z)
            t1 = time.time()
            all_runtimes['name'][index] = 'DDrppi'
            all_runtimes['nthreads'][index] = nthreads
            all_runtimes['runtime'][index] = t1 - t0
            index += 1

        if 'DD' in keys:
            t0 = time.time()
            _ = DD(autocorr, nthreads, binfile, x, y, z, x, y, z)
            t1 = time.time()
            all_runtimes['name'][index] = 'DD'
            all_runtimes['nthreads'][index] = nthreads
            all_runtimes['runtime'][index] = t1 - t0
            index += 1

        if 'wp' in keys:
            t0 = time.time()
            _ = wp(boxsize, pimax, nthreads, binfile, x, y, z)
            t1 = time.time()
            all_runtimes['name'][index] = 'wp'
            all_runtimes['nthreads'][index] = nthreads
            all_runtimes['runtime'][index] = t1 - t0
            index += 1

        if 'xi' in keys:
            t0 = time.time()
            _ = xi(boxsize, nthreads, binfile, x, y, z)
            t1 = time.time()
            all_runtimes['name'][index] = 'xi'
            all_runtimes['nthreads'][index] = nthreads
            all_runtimes['runtime'][index] = t1 - t0
            index += 1

    print("index = {0} totN = {1}".format(index, totN))
    return keys, all_runtimes
        

keys, all_runtimes = benchmark_theory_threads_all(keys=['wp'])
print("all_runtimes = {0}".format(all_runtimes))
nthreads = set(nthreads for nthreads in all_runtimes['nthreads'])
print("keys = {0}".format(keys))
print("nthreads = {0}".format(nthreads))
if min(nthreads) > 1:
    msg = "Can not scale to equivalent serial run. Min. nthreads "\
          "must be set to 1"
    raise ValueError(msg)

serial_run = all_runtimes[all_runtimes['nthreads'] == 1]
with open('nthreads_scaling_wp.txt', 'w') as f:
    print("##############################################################",
          file=f)
    print("#   Nthreads      ", end="", file=f)
    for k in keys:
        print("{0:12s} ".format(k), end="", file=f)
    
    print("", file=f)
    print("##############################################################",
          file=f)
    for nthread in nthreads:
        print("   {0:5d} ".format(nthread), end="", file=f)
        ind = all_runtimes['nthreads'] == nthread
        this_runtimes = all_runtimes[ind]

        efficiency = serial_run['runtime']/(nthread*this_runtimes['runtime'])
        efficiency *= 100.0
        for e in efficiency:
            print("{0:12.1f} ".format(e), end="", file=f)
        print("", file=f)
    

