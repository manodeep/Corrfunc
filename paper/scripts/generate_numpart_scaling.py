#!/usr/bin/env python

from __future__ import print_function
import numpy as np

import Corrfunc

from Corrfunc.io import read_catalog
from os.path import join as pjoin, abspath, dirname
import time
import sys
import os
import os.path as path
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

    # Can be at most 2, for cross-correlations
    # Once for dataset1, and once for dataset2
    with open(filename, 'r') as f:
        for l in f:
            if 'gridlink' in l and 'better' not in l:
                splits = l.split()
                serial_time += np.float(splits[-2])

            if l.startswith('0%'):
                splits = l.split()
                pair_time = np.float(splits[-2])

    return (serial_time, pair_time)

def benchmark_theory_threads_all(numpart_frac=[0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 
                                               0.6, 0.7, 0.8, 0.9, 1.0],
                                 nrepeats=3,
                                 keys=None,
                                 isa=None):

    from Corrfunc.theory import DD, DDrppi, wp, xi
    allkeys = [#'DD', 'DDrppi',
               'wp', 'xi']
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
    numpart_frac = np.array(numpart_frac)
    print("Benchmarking theory routines {0} for isa = {1}".format(keys, isa))
    allx, ally, allz = read_catalog()

    rmin = 0.1
    rmax = 84.0
    nbins = 20
    bins = np.logspace(np.log10(rmin),
                       np.log10(rmax),
                       nbins)
    pimax = rmax  # Set to rmax for comparisons between wp and xi

    autocorr = 1
    boxsize = 420.0
    nthreads = max_threads
    
    dtype = np.dtype([('repeat', np.int),
                      ('name', 'S16'),
                      ('isa', 'S16'),
                      ('rmax', np.float),
                      ('ndata',np.int),
                      ('nrand',np.int),
                      ('nthreads', np.int),
                      ('runtime', np.float),
                      ('serial_time', np.float),
                      ('pair_time', np.float),
                      ('api_time', np.float)])

    totN = len(numpart_frac) * len(keys) * len(isa) * nrepeats
    runtimes = np.empty(totN, dtype=dtype)
    runtimes['nthreads'][:] = nthreads
    runtimes['rmax'][:] = rmax


    index = 0
    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for frac in numpart_frac:
            npts = np.int(frac * len(allx))
            print("Working with N = {0}".format(npts), file=sys.stderr)

            x = np.random.choice(allx, npts, replace=False)
            y = np.random.choice(ally, npts, replace=False)
            z = np.random.choice(allz, npts, replace=False)

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
                        runtimes['repeat'][index] = repeat
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts
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
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts
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
                        runtimes['repeat'][index] = repeat
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts
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
                        runtimes['repeat'][index] = repeat
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts
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
    # autocorr is always 1 for theory routines -> 'nrand' == 'ndata'
    runtimes['nrand'][:] = (runtimes['ndata'][:]).copy()
    return keys, isa, runtimes


def benchmark_mocks_threads_all(numpart_frac=[0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 
                                              0.6, 0.7, 0.8, 0.9, 1.0],
                                nrepeats=3,
                                keys=None,
                                isa=None):
    from Corrfunc.mocks import DDrppi_mocks, DDtheta_mocks
    allkeys = [#'DDrppi (DD)',
               #'DDtheta (DD)',
               'DDrppi (DR)',
               'DDtheta (DR)'
              ]
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

    print("Benchmarking mocks routines = {0} with isa = {1}".format(keys, isa))
    mocks_file = pjoin(dirname(abspath(Corrfunc.__file__)),
                       "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
    allra, alldec, allcz = read_catalog(mocks_file)

    rand_file = pjoin(dirname(abspath(Corrfunc.__file__)),
                      "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
    allrand_ra, allrand_dec, allrand_cz = read_catalog(rand_file)
    cosmology = 1
    rmin = 0.1
    rmax = 84.0
    angmax = 10.0
    nbins = 20
    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins)

    # set to rmax for easier handling of
    # scaling with number of  particles
    pimax = rmax
    angbins = np.logspace(np.log10(rmin), np.log10(angmax), nbins)

    nthreads = max_threads
    dtype = np.dtype([('repeat', np.int),
                      ('name', 'S16'),
                      ('isa', 'S16'),
                      ('rmax', np.float),
                      ('ndata',np.int),
                      ('nrand',np.int),
                      ('nthreads', np.int),
                      ('runtime', np.float),
                      ('serial_time', np.float),
                      ('pair_time', np.float),
                      ('api_time', np.float)])
    
    totN = len(numpart_frac) * len(keys) * len(isa) * nrepeats
    runtimes = np.empty(totN, dtype=dtype)
    runtimes['nthreads'][:] = nthreads
    runtimes['rmax'][:] = rmax

    index = 0
    stderr_filename = 'stderr.txt'
    for run_isa in isa:
        for frac in numpart_frac:
            npts = np.int(frac * len(allra))
            npts_rand = np.int(frac * len(allrand_ra))
            print("Working with (N, nrand) = {0} {1}".format(npts, npts_rand), file=sys.stderr)

            ra  = np.random.choice(allra, npts, replace=False)
            dec = np.random.choice(alldec, npts, replace=False)
            cz  = np.random.choice(allcz, npts, replace=False)

            rand_ra  = np.random.choice(allrand_ra, npts_rand, replace=False)
            rand_dec = np.random.choice(allrand_dec, npts_rand, replace=False)
            rand_cz  = np.random.choice(allrand_cz, npts_rand, replace=False)

            start_thread_index = index
            if 'DDtheta (DD)' in keys:
                for repeat in range(nrepeats):
                    runtimes['repeat'][index] = repeat
                    with stderr_redirected(to=stderr_filename):
                        autocorr = 1
                        t0 = time.time()
                        _, api_time = DDtheta_mocks(autocorr, nthreads,
                                                    angbins,
                                                    ra, dec,
                                                    verbose=True,
                                                    c_api_timer=True,
                                                    isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDtheta (DD)'
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts
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
                                                    angbins,
                                                    ra, dec,
                                                    RA2=rand_ra,
                                                    DEC2=rand_dec,
                                                    verbose=True,
                                                    c_api_timer=True,
                                                    isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDtheta (DR)'
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts_rand
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
                                                   rbins, ra, dec, cz,
                                                   verbose=True,
                                                   c_api_timer=True,
                                                   isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDrppi (DD)'
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts
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
                                                   rbins, ra, dec, cz,
                                                   RA2=rand_ra,
                                                   DEC2=rand_dec,
                                                   CZ2=rand_cz,
                                                   verbose=True,
                                                   c_api_timer=True,
                                                   isa=run_isa)
                        t1 = time.time()
                        runtimes['name'][index] = 'DDrppi (DR)'
                        runtimes['isa'][index] = run_isa
                        runtimes['ndata'][index] = npts
                        runtimes['nrand'][index] = npts_rand
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
    #print("Running theory benchmarks")
    #keys, isa, runtimes = benchmark_theory_threads_all(nrepeats=3)
    #np.savez('theory_scaling_numpart.npz', keys=keys, isa=isa,
    #         runtimes=runtimes)
    #print("Theory: runtimes = {0}".format(runtimes))

    print("Running mocks benchmarks")
    keys, isa, runtimes = benchmark_mocks_threads_all(nrepeats=3)
    np.savez('mocks_scaling_numpart.npz', keys=keys, isa=isa,
             runtimes=runtimes)
    print("Mocks: runtimes = {0}".format(runtimes))
    
else:
    timings_file = sys.argv[1]
    mock = 'mock' in timings_file
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
        
    all_np = np.array(sorted(list(set(runtimes['ndata']))))
    min_np = min(all_np)
    #assert (runtimes['ndata'] == runtimes['nrand']).all()  # for simplicity
    
    if 'theory' in timings_file:
        output_file = pjoin(dirname(__file__), '../tables/',
                            'timings_Mr19_numpart_theory.tex')
    else:
        output_file = pjoin(dirname(__file__), '../tables/',
                            'timings_Mr19_numpart_mocks.tex')

    with open(output_file, 'w') as f:

        print("# Nparticles ", end='')
        print(" Nparticles ", end='', file=f)
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

        for npart in all_np:
            print(" {0:5d} ".format(npart), end='')
            print(" {0:5d} ".format(npart), end='', file=f)
            for mod in keys:
                for run_isa in isa:
                    ind = (runtimes['ndata'] == npart) & \
                          (runtimes['name'] == mod) & \
                          (runtimes['isa'] == run_isa)
                    s_ind = (runtimes['ndata'] == min_np) & \
                            (runtimes['name'] == mod) & \
                            (runtimes['isa'] == run_isa)

                    # Note, this is only looking at the pair-counting time
                    # and ignoring the total runtime ~ pair_time + serial_time
                    # For the mocks, serial_time (~0.08 sec) *is* the limiting
                    # factor in efficiency.
                    serial_avg = (runtimes['pair_time'][s_ind]).mean()
                    para_avg = (runtimes['pair_time'][ind]).mean()

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
                    print("{0:12.1f} ".format(para_avg),
                          end='')
                    print("& {0:12.1f} ".format(para_avg),
                          end='', file=f)

                print("    |", end='')

            print("")
            print("\\\\", file=f)
    
    # Begin plotting
    plt_scaling = False  # do we ever want to plot scalings for Npart?
    import matplotlib.pyplot as plt
    import seaborn
    seaborn.set_style('ticks')
    seaborn.set_style({"xtick.direction": "in","ytick.direction": "in", 'xtick.top':True, 'ytick.right':True})
    seaborn.set_context('paper')
    seaborn.set_palette('Dark2')
    
    assert len(keys) == 2
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, squeeze=False, figsize=(5,2.5))
    for ax in axes.T:
        ax[-1].set_xlabel(r'$N_\mathrm{particles}$')
    for ax in axes:
        ax[0].set_ylabel(r'Scaling efficiency' if plt_scaling else 'Runtime [sec]')
    fig.subplots_adjust(hspace=0, wspace=0)
    axes = axes.reshape(-1)
    plt_mod = {'DDrppi':r'$\mathrm{DD}(r_p,\pi)$', 'DD':r'$\mathrm{DD}(r)$', 'wp':r'$w_p(r_p)$', 'xi':r'$\xi(r)$',
               'DDrppi (DD)':r'$\mathrm{DD}(r_p,\pi)$', 'DDtheta (DD)':r'$\mathrm{DD}(\theta)$',
               'DDrppi (DR)':r'$\mathrm{DR}(r_p,\pi)$', 'DDtheta (DR)':r'$\mathrm{DR}(\theta)$'}
    for mod,ax in zip(keys,axes):
        ax.set_title(plt_mod[mod], position=(0.1,0.8), loc='left')
        if mock:
            ax.set_xlim(6e1, 2e5)
            ax.set_ylim(1e-3,500.)
        else:
            ax.set_xlim(1e3, 2e6)
            ax.set_ylim(1e-3,500.)
        plt_isa = {'avx':'AVX', 'sse42':'SSE 4.2', 'fallback':'Fallback'}
        for run_isa in isa:
            rt = []
            s_ind = (runtimes['ndata'] == min_np) & \
                    (runtimes['name'] == mod) & \
                    (runtimes['isa'] == run_isa)
            serial_avg = (runtimes['runtime'][s_ind]).mean()
            
            for npart in all_np:
                ind = (runtimes['ndata'] == npart) & \
                      (runtimes['name'] == mod) & \
                      (runtimes['isa'] == run_isa)
                para_avg = (runtimes['runtime'][ind]).mean()
                rt += [serial_avg/it/para_avg if plt_scaling else para_avg]
            
            ax.loglog(all_np, rt, label=plt_isa[run_isa])
            
            if 'avx' in run_isa:
                pltx = np.concatenate([all_np, [all_np[-1]*10]])
                plty = rt[-1]*(.6*pltx/pltx[-2])**2.
                ax.loglog(pltx, plty, ':', c='k')
        
    if mock:
        axes[0].annotate(r'$\propto N^2$', xy=(.4, .1), xycoords='axes fraction')
        axes[1].legend(loc='upper right')
    else:
        axes[0].annotate(r'$\propto N^2$', xy=(.55, .1), xycoords='axes fraction')
    #axes[0].get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))
    
    fig_fn = '{}.pdf'.format('.'.join(path.basename(timings_file).split('.')[:-1]))
    fig.tight_layout()
    fig.savefig(fig_fn)