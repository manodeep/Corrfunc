
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

            
def _get_times(filename='stderr.txt'):
    serial_time = 0.0
    pair_time = 0.0
    bruteforce = False

    # Can be at most 2, for cross-correlations
    # Once for dataset1, and once for dataset2
    with open(filename, 'r') as f:
        for l in f:
            l = l.lower()
            if 'gridlink' in l and ('better' not in l and
                                    'error' not in l):
                splits = l.split()
                try:
                    serial_time += np.float(splits[-2])
                except:
                    print("Could not convert string = '{0}' "
                          "into time".format(l))
                    
                    raise
            if 'brute' in l:
                bruteforce = True

            if l.startswith('0%'):
                splits = l.split()
                try:
                    pair_time = np.float(splits[-2])
                except:
                    print("Could not convert string = {0} "
                          "into time".format(l))
                    raise

            if bruteforce:
                serial_time = pair_time

    return (serial_time, pair_time)


def benchmark_theory_threads_all(rmax_array=[5.0, 10.0, 15.0, 20.0, 25.0,
                                             40.0, 50.0, 60.0, 70.0, 80.0],
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
    rmax_array = np.array(rmax_array)
    print("Benchmarking theory routines {0} for isa = {1}".format(keys, isa))
    x, y, z = read_catalog()

    rmin = 0.1
    nbins = 20
    autocorr = 1
    boxsize = 420.0
    nthreads = max_threads
    
    dtype = np.dtype([('repeat', np.int),
                      ('name', 'S16'),
                      ('isa', 'S16'),
                      ('rmax', np.float),
                      ('nthreads', np.int),
                      ('runtime', np.float),
                      ('serial_time', np.float),
                      ('pair_time', np.float),
                      ('api_time', np.float)])

    totN = len(rmax_array) * len(keys) * len(isa) * nrepeats
    runtimes = np.empty(totN, dtype=dtype)
    index = 0
    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for rmax in rmax_array:
            bins = np.logspace(np.log10(rmin),
                               np.log10(rmax),
                               nbins)
            pimax = rmax  # Set to rmax for comparisons between wp and xi
            print("Working on rmax = {0}".format(rmax), file=sys.stderr)
            start_thread_index = index
            if 'DD' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        t0 = time.time()
                        _, api_time = DD(autocorr, nthreads, bins, x, y, z,
                                         verbose=True, c_api_timer=True,
                                         isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DD'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            if 'DDrppi' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        t0 = time.time()
                        _, api_time = DDrppi(autocorr, nthreads, pimax,
                                             bins, x, y, z,
                                             verbose=True, c_api_timer=True,
                                             isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDrppi'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            if 'wp' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        t0 = time.time()
                        _, api_time = wp(boxsize, pimax, nthreads,
                                         bins, x, y, z,
                                         verbose=True, c_api_timer=True,
                                         isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'wp'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            if 'xi' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        t0 = time.time()
                        _, api_time = xi(boxsize, nthreads, bins, x, y, z,
                                         verbose=True, c_api_timer=True,
                                         isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'xi'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            print("{0}".format(runtimes[start_thread_index:index]))
            sys.stdout.flush()
            
    print("index = {0} totN = {1}".format(index, totN))
    return keys, isa, runtimes


def benchmark_mocks_threads_all(rmax_array=[5.0, 10.0, 15.0, 20.0, 25.0,
                                            40.0, 50.0, 60.0, 70.0, 80.0],
                                nrepeats=1,
                                keys=None,
                                isa=None):
    from Corrfunc.mocks import DDrppi_mocks, DDtheta_mocks
    allkeys = ['DDrppi (DD)',
               'DDtheta (DD)',
               'DDrppi (DR)',
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

    rmax_array = np.array(rmax_array)
    print("Benchmarking mocks routines = {0} with isa = {1}".format(keys, isa))
    mocks_file = pjoin(dirname(abspath(Corrfunc.__file__)),
                       "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
    ra, dec, cz = read_catalog(mocks_file)

    rand_file = pjoin(dirname(abspath(Corrfunc.__file__)),
                      "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
    rand_ra, rand_dec, rand_cz = read_catalog(rand_file)
    cosmology = 1
    rmin = 0.1
    nbins = 20
    nthreads = max_threads
    
    dtype = np.dtype([('repeat', np.int),
                      ('name', 'S16'),
                      ('isa', 'S16'),
                      ('rmax', np.float),
                      ('nthreads', np.int),
                      ('runtime', np.float),
                      ('serial_time', np.float),
                      ('pair_time', np.float),
                      ('api_time', np.float)])
    
    totN = len(rmax_array) * len(keys) * len(isa) * nrepeats
    runtimes = np.empty(totN, dtype=dtype)
    index = 0
    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for rmax in rmax_array:
            bins = np.logspace(np.log10(rmin),
                               np.log10(rmax),
                               nbins)
            pimax = rmax  # Set to rmax for comparisons between wp and xi
            print("Working on rmax = {0}".format(rmax), file=sys.stderr)
            start_thread_index = index
            if 'DDtheta (DD)' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 1
                        t0 = time.time()
                        _, api_time = DDtheta_mocks(autocorr, nthreads,
                                                    bins,
                                                    ra, dec,
                                                    verbose=True,
                                                    c_api_timer=True,
                                                    isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDtheta (DD)'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            if 'DDtheta (DR)' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 0
                        t0 = time.time()
                        _, api_time = DDtheta_mocks(autocorr, nthreads,
                                                    bins,
                                                    ra, dec,
                                                    RA2=rand_ra,
                                                    DEC2=rand_dec,
                                                    verbose=True,
                                                    c_api_timer=True,
                                                    isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDtheta (DR)'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1
                        
            if 'DDrppi (DD)' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
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
                        runtimes['name'][index] = 'DDrppi (DD)'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            if 'DDrppi (DR)' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
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
                        runtimes['name'][index] = 'DDrppi (DR)'
                        runtimes['isa'][index] = run_isa
                        runtimes['rmax'][index] = rmax
                        runtimes['nthreads'][index] = nthreads
                        runtimes['runtime'][index] = t1 - t0
                        serial_time, pair_time = _get_times(stderr_filename)
                        runtimes['serial_time'][index] = serial_time
                        runtimes['pair_time'][index] = pair_time
                        runtimes['api_time'][index] = api_time
                        index += 1

            print("{0}".format(runtimes[start_thread_index:index]))
            sys.stdout.flush()
            
    print("index = {0} totN = {1}".format(index, totN))
    return keys, isa, runtimes

if len(sys.argv) == 1:
    print("Running theory benchmarks")
    keys, isa, runtimes = benchmark_theory_threads_all(nrepeats=5)
    np.savez('theory_scaling_rmax.npz', keys=keys, isa=isa,
             runtimes=runtimes)
    print("Theory: runtimes = {0}".format(runtimes))

    print("Running mocks benchmarks")
    keys, isa, runtimes = benchmark_mocks_threads_all(nrepeats=5)
    np.savez('mocks_scaling_rmax.npz', keys=keys, isa=isa,
             runtimes=runtimes)
    print("Mocks: runtimes = {0}".format(runtimes))
    
else:
    timings_file = sys.argv[1]
    print("Loading benchmarks from file = {0}".format(timings_file))
    xx = np.load(timings_file)
    try:
        keys = xx['keys']
        isa = xx['isa']
        try:
            runtimes = xx['runtimes']
        except KeyError:
            # Previous versions of this script used 'all_runtimes'
            runtimes = xx['all_runtimes']
            
    except KeyError:
        print("Error: Invalid timings file = `{0}' passed in the "
              "command-line ".format(timings_file))
        raise
    
    nthreads = set(nthreads for nthreads in runtimes['nthreads'])
    if min(nthreads) > 1:
        msg = "Can not scale to equivalent serial run. Min. nthreads "\
              "must be set to 1"
        raise ValueError(msg)

    if 'theory' in timings_file:
        output_file = pjoin(dirname(__file__), '../tables/',
                            'timings_Mr19_rmax_theory.tex')
    else:
        output_file = pjoin(dirname(__file__), '../tables/',
                            'timings_Mr19_rmax_mocks.tex')

    with open(output_file, 'w') as f:

        print("# Nthreads ", end='')
        print(" Nthreads ", end='', file=f)
        for mod in keys:
            print("{0:16s}{1:9s}{0:16s}".format("", mod), end='')
            print("& {0:16s}{1:9s}{0:16s}".format("", mod), end='', file=f)
        print("")
        print("\\\\", file=f)
        print("# {0:11s} ".format(""), end='')
        print(" {0:11s} ".format(""), end='', file=f)
        for mod in keys:
            for run_isa in isa:
                print("{0:12s}".format(run_isa), end='')
                print("& {0:12s}".format(run_isa), end='', file=f)
            print("     ", end='')
            print("     ", end='', file=f)

        print("")
        print("\\\\", file=f)

        for it in nthreads:
            print(" {0:5d} ".format(it), end='')
            print(" {0:5d} ".format(it), end='', file=f)
            for mod in keys:
                for run_isa in isa:
                    ind = (runtimes['nthreads'] == it) & \
                          (runtimes['name'] == mod) & \
                          (runtimes['isa'] == run_isa)
                    s_ind = (runtimes['nthreads'] == 1) & \
                            (runtimes['name'] == mod) & \
                            (runtimes['isa'] == run_isa)

                    # Note, this is only looking at the pair-counting time
                    # and ignoring the total runtime ~ pair_time + serial_time
                    # For the mocks, serial_time (~0.08 sec) *is* the limiting
                    # factor in efficiency.
                    serial_avg = np.mean(runtimes['pair_time'][s_ind])
                    para_avg = np.mean(runtimes['pair_time'][ind])

                    # If you want the Amdahl's law limited max. effiency,
                    # uncomment the following lines.
                    # serial_avg = np.mean(runtimes['runtime'][s_ind])
                    # para_avg = np.mean(runtimes['runtime'][ind])
                    # serial_para_runs = np.mean(runtimes['serial_time'][ind])
                    # serial_serial_runs = np.mean(runtimes['serial_time'][s_ind])
                    # theoretical_best_time = serial_runtime_in_parallel_runs + \
                    #                         (serial_avg-serial_runtime_in_serial_runs)/it
                    # print("{0:9.1f}({1:3.1f})".format((serial_avg/it)/para_avg*100,
                    #                               (serial_avg/it)/theoretical_best_time*100.0),
                    #       end='')
                    print("{0:12.1f}".format((serial_avg/it)/para_avg*100),
                          end='')
                    print("& {0:12.1f} ".format(serial_avg/it/para_avg*100),
                          end='', file=f)

                print("    |", end='')

            print("")
            print("\\\\", file=f)
    
