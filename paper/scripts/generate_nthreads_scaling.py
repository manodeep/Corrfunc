
from __future__ import print_function
import numpy as np

import Corrfunc

from Corrfunc.io import read_catalog
from os.path import join as pjoin, abspath, dirname
import time
import sys
import os
# Taken directly from answer by  J.F. Sebastian
# http://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python/17954769#17954769
from contextlib import contextmanager

import multiprocessing
max_threads = multiprocessing.cpu_count()


@contextmanager
def stderr_redirected(to=os.devnull):
    """
    import os
    
    with stderr_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    """
    fd = sys.stderr.fileno()

    # assert that Python and C stdio write using the same file descriptor
    # assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stderr")) == fd == 1
    
    def _redirect_stderr(to):
        sys.stderr.close()  # + implicit flush()
        os.dup2(to.fileno(), fd)  # fd writes to 'to' file
        sys.stderr = os.fdopen(fd, 'w')  # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stderr:
        with open(to, 'w') as file:
            _redirect_stderr(to=file)
        try:
            yield  # allow code to be run with the redirected stderr
        finally:
            _redirect_stderr(to=old_stderr)    # restore stderr
            # buffering and flags such as
            # CLOEXEC may be different
            
    
def benchmark_theory_threads_all(min_threads=1, max_threads=max_threads,
                                 nrepeats=1,
                                 keys=None,
                                 isa=None):

    from Corrfunc.theory import DD, DDrppi, wp, xi
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
        with open(filename, 'r') as f:
            for l in f:
                if 'gridlink' in l:
                    splits = l.split()
                    serial_time = np.float(splits[-2])
                    break

        return serial_time

    print("Benchmarking theory routines = {0} with isa = {1}".format(keys,
                                                                     isa))
    x, y, z = read_catalog()
    rmax = 10.0
    rmin = 0.1
    nbins = 20
    bins = np.logspace(np.log10(rmin),
                       np.log10(rmax),
                       nbins)
    autocorr = 1
    pimax = rmax  # Set to rmax for comparisons between wp and xi
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
    index = 0
    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for nthreads in range(min_threads, max_threads + 1):
            print("Working on nthreads = {0}".format(nthreads),
                  file=sys.stderr)
            start_thread_index = index
            if 'DD' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        t0 = time.time()
                        _, api_time = DD(autocorr, nthreads, bins, x, y, z,
                                         verbose=True, c_api_timer=True,
                                         isa=run_isa)
                        t1 = time.time()
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
                        _, api_time = DDrppi(autocorr, nthreads, pimax,
                                             bins, x, y, z,
                                             verbose=True, c_api_timer=True,
                                             isa=run_isa)
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
                        _, api_time = wp(boxsize, pimax, nthreads,
                                         bins, x, y, z,
                                         verbose=True, c_api_timer=True,
                                         isa=run_isa)
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
                        _, api_time = xi(boxsize, nthreads, bins, x, y, z,
                                         verbose=True, c_api_timer=True,
                                         isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'xi'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)                
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            print("{0}".format(all_runtimes[start_thread_index:index]))
            sys.stdout.flush()
            
    print("index = {0} totN = {1}".format(index, totN))
    return keys, isa, all_runtimes


def benchmark_mocks_threads_all(min_threads=1, max_threads=max_threads,
                                nrepeats=1,
                                keys=None,
                                isa=None):
    from Corrfunc.mocks import DDrppi_mocks, DDtheta_mocks
    allkeys = ['DDrppi (DD)',
               'DDrppi (DR)',
               'DDtheta (DD)',
               'DDtheta (DR)']
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
        with open(filename, 'r') as f:
            for l in f:
                if 'gridlink' in l:
                    splits = l.split()
                    serial_time = np.float(splits[-2])
                    break

        return serial_time

    print("Benchmarking mocks routines = {0} with isa = {1}".format(keys, isa))
    mocks_file = pjoin(dirname(abspath(Corrfunc.__file__)),
                       "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
    ra, dec, cz = read_catalog(mocks_file)

    rand_file = pjoin(dirname(abspath(Corrfunc.__file__)),
                      "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
    rand_ra, rand_dec, rand_cz = read_catalog(rand_file)
    cosmology = 1
    nbins = 20
    rmin = 0.1
    rmax = 10.0
    bins = np.logspace(np.log10(rmin), np.log10(rmax), nbins)
    pimax = rmax    # set to rmax for easier handling of
                    # scaling with number of  particles
    dtype = np.dtype([('repeat', np.int),
                      ('name', 'S16'),
                      ('isa', 'S16'),
                      ('nthreads', np.int),
                      ('runtime', np.float),
                      ('serial_time', np.float),
                      ('api_time', np.float)])

    totN = (max_threads - min_threads + 1) * len(keys) * len(isa) * nrepeats
    all_runtimes = np.empty(totN, dtype=dtype)
    index = 0
    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for nthreads in range(min_threads, max_threads + 1):
            print("Working on nthreads = {0}".format(nthreads),
                  file=sys.stderr)
            start_thread_index = index
            if 'DDtheta (DD)' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 1
                        t0 = time.time()
                        _, api_time = DDtheta_mocks(autocorr, nthreads, bins,
                                                    ra, dec,
                                                    verbose=True,
                                                    c_api_timer=True,
                                                    isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'DDtheta (DD)'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            if 'DDtheta (DR)' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 0
                        t0 = time.time()
                        _, api_time = DDtheta_mocks(autocorr, nthreads, bins,
                                                    ra, dec,
                                                    RA2=rand_ra,
                                                    DEC2=rand_dec,
                                                    verbose=True,
                                                    c_api_timer=True,
                                                    isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'DDtheta (DR)'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)
                        all_runtimes['api_time'][index] = api_time
                        index += 1
                        
            if 'DDrppi (DD)' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 1
                        t0 = time.time()
                        _, api_time = DDrppi_mocks(autocorr, cosmology,
                                                   nthreads, pimax,
                                                   bins, ra, dec, cz,
                                                   verbose=True,
                                                   c_api_timer=True,
                                                   isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'DDrppi (DD)'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            if 'DDrppi (DR)' in keys:
                for repeat in range(nrepeats):
                    all_runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 0
                        t0 = time.time()
                        _, api_time = DDrppi_mocks(autocorr, cosmology,
                                                   nthreads, pimax,
                                                   bins, ra, dec, cz,
                                                   RA2=rand_ra,
                                                   DEC2=rand_dec,
                                                   CZ2=rand_cz,
                                                   verbose=True,
                                                   c_api_timer=True,
                                                   isa=run_isa)
                        t1 = time.time()
                        all_runtimes['name'][index] = 'DDrppi (DR)'
                        all_runtimes['isa'][index] = run_isa
                        all_runtimes['nthreads'][index] = nthreads
                        all_runtimes['runtime'][index] = t1 - t0
                        all_runtimes['serial_time'][index] = _get_time_from_stderr(stderr_filename)
                        all_runtimes['api_time'][index] = api_time
                        index += 1

            print("{0}".format(all_runtimes[start_thread_index:index]))
            sys.stdout.flush()
            
    print("index = {0} totN = {1}".format(index, totN))
    return keys, isa, all_runtimes

if len(sys.argv) == 1:
    # print("Running theory benchmarks")
    # keys, isa, all_runtimes = benchmark_theory_threads_all(nrepeats=10)
    # np.savez('theory_runtimes.npz', keys=keys, isa=isa,
    #          all_runtimes=all_runtimes)
    # print("Theory: all_runtimes = {0}".format(all_runtimes))

    print("Running mocks benchmarks")
    keys, isa, all_runtimes = benchmark_mocks_threads_all(nrepeats=10)
    np.savez('mocks_runtimes.npz', keys=keys, isa=isa,
             all_runtimes=all_runtimes)
    print("Mocks: all_runtimes = {0}".format(all_runtimes))
    
else:
    timings_file = sys.argv[1]
    print("Loading benchmarks from file = {0}".format(timings_file))
    xx = np.load(timings_file)
    try:
        keys = xx['keys']
        isa = xx['isa']
        all_runtimes = xx['all_runtimes']
    except KeyError:
        print("Error: Invalid timings file = `{0}' passed in the "
              "command-line ".format(timings_file))
        raise
    
    nthreads = set(nthreads for nthreads in all_runtimes['nthreads'])
    if min(nthreads) > 1:
        msg = "Can not scale to equivalent serial run. Min. nthreads "\
              "must be set to 1"
        raise ValueError(msg)

    print("# Nthreads ", end='')
    for mod in keys:
        print("{0:16s}{1:9s}{0:16s}".format("", mod), end='')
    print("")
    print("# {0:11s} ".format(""), end='')
    for mod in keys:
        for run_isa in isa:
            print("{0:12s}".format(run_isa), end='')
        print("     ", end='')
        
    print("")
    
    for it in nthreads:
        print(" {0:5d} ".format(it), end='')
        for mod in keys:
            for run_isa in isa:
                ind = (all_runtimes['nthreads'] == it) & \
                      (all_runtimes['name'] == mod) & \
                      (all_runtimes['isa'] == run_isa)
                serial_ind = (all_runtimes['nthreads'] == 1) & \
                             (all_runtimes['name'] == mod) & \
                             (all_runtimes['isa'] == run_isa)
                serial_avg = np.mean(all_runtimes['runtime'][serial_ind])
                parallel_avg = np.mean(all_runtimes['runtime'][ind])

                # If you want the Amdahl's law limited max. effiency, uncomment
                # the following lines.
                # serial_runtime_in_parallel_runs = np.mean(all_runtimes['serial_time'][ind])
                # serial_runtime_in_serial_runs = np.mean(all_runtimes['serial_time'][serial_ind])
                # theoretical_best_time = serial_runtime_in_parallel_runs + \
                #                         (serial_avg-serial_runtime_in_serial_runs)/it
                # print("{0:7.1f}({1:3.1f})".format((serial_avg/it)/parallel_avg*100,
                #                               (serial_avg/it)/theoretical_best_time*100.0),
                #       end='')
                print("{0:12.1f}".format((serial_avg/it)/parallel_avg*100),
                      end='')
                
            print("    |", end='')
            
        print("")
        
    # with open('nthreads_scaling_wp.txt', 'w') as f:
    #     print("##############################################################",
    #           file=f)
    #     print("#   Nthreads      ", end="", file=f)
    #     for k in keys:
    #         print("{0:12s} ".format(k), end="", file=f)

    #     print("", file=f)
    #     print("##############################################################",
    #           file=f)
    #     for nthread in nthreads:
    #         print("   {0:5d} ".format(nthread), end="", file=f)
    #         ind = all_runtimes['nthreads'] == nthread
    #         this_runtimes = all_runtimes[ind]

    #         efficiency = serial_run['runtime']/(nthread*this_runtimes['runtime'])
    #         efficiency *= 100.0
    #         for e in efficiency:
    #             print("{0:12.1f} ".format(e), end="", file=f)
    #         print("", file=f)


