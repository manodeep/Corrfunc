
from __future__ import print_function
import numpy as np

import Corrfunc
from Corrfunc.theory import DD, DDrppi, wp, xi

from Corrfunc.io import read_catalog
from os.path import join as pjoin, abspath, dirname
import time
import sys
import os
    

# Taken directly from answer by  J.F. Sebastian 
# http://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python/17954769#17954769
from contextlib import contextmanager
@contextmanager
def stderr_redirected(to=os.devnull):
    """
    import os
    
    with stderr_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    """
    fd = sys.stderr.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stderr")) == fd == 1
    
    def _redirect_stderr(to):
        sys.stderr.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stderr = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stderr:
        with open(to, 'w') as file:
            _redirect_stderr(to=file)
        try:
            yield # allow code to be run with the redirected stderr
        finally:
            _redirect_stderr(to=old_stderr) # restore stderr.
            # buffering and flags such as
            # CLOEXEC may be different
            
    
def benchmark_theory_threads_all(min_threads=1, max_threads=16,
                                 nrepeats=1,
                                 keys=None,
                                 isa=None):

    allkeys = ['DDrppi', 'DD', 'wp', 'xi']
    allisa = ['avx', 'sse42', 'fallback']
    if keys is None:
        keys = allkeys
    else:
        for k in keys:
            if k not in allkeys:
                msg = "Valid routines to benchmark are: {0}\nFound routine"\
                    " = {1}".format(allkeys, k)
                raise ValueError(msg)

    if isa is None:
        isa = allisa
    else:
        for i in isa:
            if i not in allisa:
                msg = "Valid instructions sets benchmark are: {0}\n"\
                      "Found routine = {1}".format(allisa, i)
                raise ValueError(msg)

            
    def _get_time_from_stderr(filename='stderr.txt'):
        serial_time = 0.0
        with open(filename,'r') as f:
            for l in f:
                if 'gridlink' in l:
                    splits = l.split()
                    serial_time = np.float(splits[-2])
                    break

        return serial_time

            
    print("Benchmarking routines = {0} with isa = {1}".format(keys, isa))
    x, y, z = read_catalog()
    binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
                    "../theory/tests/", "bins")
    autocorr = 1
    pimax = 40.0
    boxsize = 420.0
    dtype = np.dtype([('repeat', np.int),
                      ('name', 'S16'),
                      ('isa', 'S16'),
                      ('nthreads', np.int),
                      ('runtime', np.float),
                      ('serial_time', np.float),
                      ('api_time', np.float)])

    totN = (max_threads - min_threads + 1) * len(keys) * len(isa) * nrepeats
    all_runtimes = np.empty(totN, dtype=dtype)
    index = 0L

    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for nthreads in range(min_threads, max_threads + 1):
            print("Working on nthreads = {0}".format(nthreads), file=sys.stderr)
            start_thread_index = index
            if 'DD' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        t0 = time.time()
                        _, api_time = DD(autocorr, nthreads, binfile, x, y, z,
                                         verbose=True, c_api_timer=True, isa=run_isa)
                        t1 = time.time()
                        serial_time = _get_time_from_stderr(stderr_filename)
                        all_runtimes['name'][index] = 'DD'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            
            if 'DDrppi' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):                
                        t0 = time.time()
                        _, api_time = DDrppi(autocorr, nthreads, pimax, binfile, x, y, z,
                                             verbose=True, c_api_timer=True, isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'DDrppi'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            if 'wp' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):                                
                        t0 = time.time()
                        _, api_time = wp(boxsize, pimax, nthreads, binfile, x, y, z,
                                         verbose=True, c_api_timer=True, isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'wp'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)                
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            if 'xi' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):                
                        t0 = time.time()
                        _, api_time = xi(boxsize, nthreads, binfile, x, y, z,
                                         verbose=True, c_api_timer=True, isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'xi'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)                
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            print("{0}".format(all_runtimes[start_thread_index:index]))
            
    print("index = {0} totN = {1}".format(index, totN))
    return keys, isa, all_runtimes


keys, isa, all_runtimes = benchmark_theory_threads_all(nrepeats=10)
np.savez('runtimes.npz',keys=keys,isa=isa,all_runtimes=all_runtimes)
print("all_runtimes = {0}".format(all_runtimes))
nthreads = set(nthreads for nthreads in all_runtimes['nthreads'])
print("keys = {0}".format(keys))
print("isa = {0}".format(isa))
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


