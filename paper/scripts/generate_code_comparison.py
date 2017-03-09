#!/usr/bin/env python
'''
A script to do timing tests of different codes with Corrfunc.
In each case, we endeavor to perform the same calculations
with Corrfunc that each code does internally, thus producing
a fair comparison.  We also check that the results of each code
match those of Corrfunc.
'''

import time
import argparse
from collections import defaultdict
import numpy as np
import os
import os.path as path
import Corrfunc
import subprocess32
import re

import multiprocessing
max_threads = multiprocessing.cpu_count()

registry = {}  # all serial functions
parallel_registry = {}  # specifically multithreaded functions

# a decorator to register codes to test
def register_code(name, serial=True, multithreaded=False):
    def wrapper(test_func):
        if multithreaded:
            parallel_registry[name] = test_func
        if serial:
            registry[name] = test_func
        return test_func
    return wrapper

# Begin function definitions

import scipy.spatial
@register_code('SciPy cKDTree')
def test_scipy_ckdtree(points, bins):
    t = time.time()
    tree = scipy.spatial.cKDTree(points, leafsize=16, compact_nodes=True, balanced_tree=True)
    results = tree.count_neighbors(tree, bins)
    t = time.time() - t
    
    # cKDTree returns all pairs <= r, not pairs in bins.  Include a bin that goes all the way to r~0.
    bins = [0.] + list(bins)
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=False, nthreads=1, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                                    periodic=False, output_ravg=False, verbose=False)
    cf_t = time.time() - cf_t
    cum_counts = cf_results['npairs'].cumsum()
    
    assert np.all(results == cum_counts), (results,cum_counts)
    
    return cf_t, t


import sklearn.neighbors
@register_code('scikit-learn KDTree')
def test_sklearn_kdtree(points, bins):
    t = time.time()
    tree = sklearn.neighbors.KDTree(points, leaf_size=40)
    results = tree.two_point_correlation(points, bins, dualtree=True)  # dualtree is faster in our tests
    t = time.time() - t
    
    # scikit-learn returns all pairs <= r, not pairs in bins.  Include a bin that goes all the way to r~0.
    bins = [0.] + list(bins)
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=False, nthreads=1, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                                    periodic=False, output_ravg=False, verbose=False)
    cf_t = time.time() - cf_t
    cum_counts = cf_results['npairs'].cumsum()
    assert np.all(results == cum_counts), (results,cum_counts)
    
    return cf_t, t
    

import kdcount.correlate
import kdcount.models
@register_code('kdcount', multithreaded=True)
def test_kdcount(points, bins, nthreads=1):
    kdpoints = kdcount.models.points(points)
    kdbins = kdcount.correlate.RBinning(bins)
    t = time.time()
    counts = kdcount.correlate.paircount(kdpoints, kdpoints, kdbins, np=nthreads)
    t = time.time() - t
    counts = counts.sum1
    
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=False, nthreads=nthreads, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                                    periodic=False, output_ravg=False, verbose=False)
    cf_t = time.time() - cf_t
    cf_npairs = cf_results['npairs']
    assert np.all(counts == cf_npairs), (counts,cf_npairs)
    
    return cf_t, t
    
import halotools.mock_observables
@register_code('halotools', multithreaded=True)
def test_halotools(points, bins, nthreads=1):
    # even though we have to pass randoms, RR_precomputed should avoid evaluating this
    RR_pre = np.ones_like(bins[:-1])
    t = time.time()
    xi = halotools.mock_observables.tpcf(points, bins, num_threads=nthreads, max_sample_size=len(points),
                                         RR_precomputed=RR_pre, NR_precomputed=len(points), randoms=points)
    t = time.time() - t
    paircounts = np.round(xi + 1)  # could be affected by roundoff error
    
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=False, nthreads=nthreads, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                                    periodic=False, output_ravg=False, verbose=False)
    cf_t = time.time() - cf_t
    cf_npairs = cf_results['npairs']
    assert np.all(paircounts == cf_npairs), (paircounts,cf_npairs)
    
    return cf_t, t
    

import treecorr
@register_code('Treecorr', multithreaded=True)
def test_treecorr(points, bins, nthreads=1):
    cat = treecorr.Catalog(x=points[:,0], y=points[:,1], z=points[:,2])
    config = {'nbins':len(bins)-1, 'min_sep':bins[0], 'max_sep':bins[-1], 'num_threads':nthreads, 'bin_slop':0.}
    nn = treecorr.NNCorrelation(config=config)
    t = time.time()
    nn.process(cat)
    t = time.time() - t
    tc_npairs = nn.npairs*2  # Corrfunc gives doubled counts
    
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=True, nthreads=nthreads, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    #X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                                    periodic=False, output_ravg=True, verbose=False)
    cf_t = time.time() - cf_t
    cf_npairs = cf_results['npairs']
    assert np.all(tc_npairs == cf_npairs), (tc_npairs,cf_npairs)
    
    return cf_t, t

@register_code('CUTE_box', multithreaded=True)
def test_cute_box(points, bins, nthreads=1):
    # Cute prefers to use linear binning, so we'll discard the input bins and use our own
    bins = np.linspace(0., bins[-1], len(bins))
    
    np.savetxt('cute_box_data.dat', points)
    
    path_to_cute = '/home/lgarrison/CUTE/CUTE_box/CUTE_box'
    param_fn = 'cute_box_param.txt'
    
    output_fn = 'cute_box_out.txt'
    os.environ['OMP_NUM_THREADS'] = str(nthreads)
    stdout = subprocess32.check_output([path_to_cute, param_fn])
    del os.environ['OMP_NUM_THREADS']
    #print stdout
    outfile = np.loadtxt(output_fn)
    cute_counts = outfile[:,-1]
    
    # parse the output to get the processing time
    # this does not include any catalog reading time
    matches = re.search(r'Total time ellapsed (\d+(?:.\d+)?) ms', stdout)
    t = float(matches.group(1))/1000.
    
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=True, nthreads=nthreads, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    periodic=True, output_ravg=False, verbose=False,
                                    #xbin_refine_factor=7, ybin_refine_factor=7, zbin_refine_factor=1, max_cells_per_dim=1000,
                                    boxsize=1100.)
    cf_t = time.time() - cf_t
    cf_npairs = cf_results['npairs']
    cf_npairs[0] -= len(points)  # cute doesn't count self-pairs
    # CUTE is often off by a few pairs.
    # This could be due to their binning scheme, which relies on floating point math to handle decisions about bin boundaries
    assert np.allclose(cute_counts, cf_npairs, atol=10, rtol=.1), (cute_counts,cf_npairs)
    
    return cf_t, t

@register_code('mlpack RangeSearch')
def test_mlpack(points, bins):
    # range search only supports one bin, and runs out of memory for large Rmax
    bins = [bins[0], bins[-1]/2.5]
    
    mlpack_data_fn = 'mlpack_data.csv'
    np.savetxt(mlpack_data_fn, points)
    mlpack_out_fn = 'mlpack_neighbors.csv'
    mlpack_cmd = 'mlpack_range_search --min={min} --max={max} -v -r {infn} --neighbors_file={outfn}'.format(min=bins[0], max=bins[-1], infn=mlpack_data_fn, outfn=mlpack_out_fn)
    stdout = subprocess32.check_output(mlpack_cmd.split(' '))
    
    # parse the output to get the processing time
    # this does not include any catalog reading time
    # mlpack is often much slower than this suggests.  is the difference purely output time?
    t_rs = float(re.search(r'range_search/computing_neighbors: (\d+(.\d+)?)s', stdout).group(1))
    t_tree = float(re.search(r'tree_building: (\d+(.\d+)?)s', stdout).group(1))
    t = t_rs + t_tree
    
    # Read the neighbors file to count the number of pairs
    mlpack_pairs = 0
    with open(mlpack_out_fn, 'rb') as fp:
        for line in fp:
            line = line.strip()
            neighbors = (line.count(',') + 1) if line else 0
            mlpack_pairs += neighbors
    
    cf_t = time.time()
    cf_results = Corrfunc.theory.DD(autocorr=False, nthreads=1, binfile=bins,
                                    X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                                    X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                                    periodic=False, output_ravg=False, verbose=False)
    cf_t = time.time() - cf_t
    cf_npairs = cf_results['npairs']
    assert np.all(mlpack_pairs == cf_npairs), (mlpack_pairs,cf_npairs)
    
    return cf_t, t


#@register_code('swot', serial=False, multithreaded=True)  # swot is too slow to run at these problem sizes
def test_swot(points, bins, nthreads=1):
    # Swot only supports evenly spaced bins, and we can't go very big
    bins = np.linspace(bins[0]/100, bins[-1]/100., len(bins))
    
    np.savetxt('swot_test_data.txt', points)
    
    mpi = ['mpirun', '-np', str(nthreads)]
    exec_path = os.getenv('HOME') + '/swot/bin/swot'
    
    options = ['-c', 'swot_config.txt']
    output = subprocess32.check_output(mpi + [exec_path] + options, stderr=subprocess32.STDOUT)
    t = float(re.search(r'CPU wall time: (\d+(.\d+)?)s', output).group(1))
    
    swot_results = np.loadtxt('swot_corr.out')
    swot_npairs = swot_results[:,6]
    
    cf_t = time.time()
    cf_DD = Corrfunc.theory.DD(autocorr=True, nthreads=nthreads, binfile=bins,
                               X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                               X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                               periodic=False, output_ravg=True, verbose=False)
    cf_DR = Corrfunc.theory.DD(autocorr=False, nthreads=nthreads, binfile=bins,
                               X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                               X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                               periodic=False, output_ravg=True, verbose=False)
    cf_RR = Corrfunc.theory.DD(autocorr=True, nthreads=nthreads, binfile=bins,
                               X1=points[:,0], Y1=points[:,1], Z1=points[:,2],
                               X2=points[:,0], Y2=points[:,1], Z2=points[:,2],
                               periodic=False, output_ravg=True, verbose=False)
    cf_t = time.time() - cf_t
    cf_npairs = cf_DD['npairs']
    assert np.all(swot_npairs == cf_npairs), (swot_npairs,cf_npairs)
    
    return cf_t, t


# Run all of the codes in the given registry
# nthreads=None specifies to not pass a nthread value
def test_from_registry(reg, rmin=1, rmax=90., nbin=19, nthreads=max_threads, nrepeat=3,
                       np_fracs=[0.01, 0.03, 0.1, 0.25, 0.5, 1.0]
                       #np_fracs=[0.01, 0.03, 0.1]
                      ):
    #points = Corrfunc.io.read_catalog()
    #points = np.array(points).T  #.copy()  # shape (N,3)
    points = np.loadtxt('halos_emulator_1100box_Neff3_00.txt')
    assert (points >= 0).all() and (points < 1100.).all()
    dtype = points.dtype  # float64
    
    bins = np.logspace(np.log(rmin), np.log(rmax), nbin+1, base=np.e)
    
    nthreads_kwarg = {'nthreads':nthreads} if nthreads else {}
    
    # last axis is [corrfunc, other code]
    func_runtimes = {n:np.empty((len(np_fracs),nrepeat,2), dtype=float) for n in reg}
    all_np = [int(frac*len(points)) for frac in np_fracs]
    func_runtimes['all_np'] = all_np
    
    for i,numpart in enumerate(all_np):
        print('\tTesting with {} particles'.format(numpart))
        # generate the random subsample
        points = points.reshape(-1).view(dtype=[('x',dtype,3)])
        subsample = np.random.choice(points, numpart, replace=False)
        subsample = subsample.view(dtype=dtype).reshape(-1,3)
        
        for name,func in reg.iteritems():
            print('\t\t{}'.format(name))
            for j in xrange(nrepeat):
                # Do the timing
                rt = func(subsample, bins, **nthreads_kwarg)
                func_runtimes[name][i,j] = rt
        
    return func_runtimes


import matplotlib.pyplot as plt
import seaborn
from cycler import cycler
from collections import OrderedDict

def plot_results(timings_fn):
    parallel = 'parallel' in timings_fn
    seaborn.set_style('ticks')
    seaborn.set_style({"xtick.direction": "in","ytick.direction": "in", 'xtick.top':True, 'ytick.right':True})
    seaborn.set_context('paper')

    _func_runtimes = dict(np.load(timings_fn))
    all_np = _func_runtimes.pop('all_np')
    
    # Assign a persistent line style to each func
    order = ['halotools', 'kdcount', 'Treecorr', 'CUTE_box', 'scikit-learn KDTree', 'SciPy cKDTree', 'mlpack RangeSearch', 'swot']
    prop_cycle =  cycler(linestyle=['-','--']) * cycler(color=seaborn.color_palette('colorblind',4))
    prop_cycle = prop_cycle()  # make infinite
    func_props = [next(prop_cycle) for func in order]
    func_runtimes = OrderedDict()
    for func,props in zip(order, func_props):
        if func in _func_runtimes:
            func_runtimes[func] = _func_runtimes.pop(func), props
    assert len(_func_runtimes) == 0
    del _func_runtimes
    
    fig, ax = plt.subplots(1,1, figsize=(5,4))
    ax.set_xlabel(r'$N_\mathrm{particles}$')
    ax.set_ylabel(r'Corrfunc speed-up: $t_\mathrm{other}/t_\mathrm{Corrfunc}$')
    #ax.set_ylabel(r'Corrfunc time [sec]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    for func,(runtimes,props) in func_runtimes.iteritems():
        means = runtimes.mean(axis=-2)
        ax.plot(all_np, means[:,1]/means[:,0], marker='', label=func, **props)
        #ax.plot(all_np, means[:,0],label=func, marker='o', **next(prop_cycle))
        print(r'{}: {:.1f} -- {:.1f}\times'.format(func, (means[:,1]/means[:,0]).min(), (means[:,1]/means[:,0]).max()))
    ax.axhline(1., linestyle=':', c='k')
    if parallel:
        ax.legend(loc='upper left')
        ax.annotate('Corrfunc faster', xy=(4e6, 1.1), horizontalalignment='right', verticalalignment='bottom')
        ax.annotate('Corrfunc slower', xy=(4e6, 0.9), horizontalalignment='right', verticalalignment='top')
    else:
        ax.legend(loc=(.6,.13))
        ax.set_ylim(top=1e1)
        ax.annotate('Corrfunc faster', xy=(1e5, 1.05), horizontalalignment='left', verticalalignment='bottom')
        ax.annotate('Corrfunc slower', xy=(1e5, 0.95), horizontalalignment='left', verticalalignment='top')
    
    fig_fn = '{}.pdf'.format('.'.join(path.basename(timings_fn).split('.')[:-1]))
    fig.tight_layout()
    fig.savefig(fig_fn)


# Begin test driver
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results_file', help='If given, plot the saved results.  Otherwise run the tests and save the results to disk.', nargs='?')
    args = parser.parse_args()
    
    if not args.results_file:
        # Run the tests
        #print('Begin testing single-threaded codes')
        #runtimes = test_from_registry(registry, nthreads=None)
        #np.savez('codes_scaling_numpart', **runtimes)

        print('Begin testing multi-threaded codes')
        parallel_runtimes = test_from_registry(parallel_registry)
        np.savez('codes_scaling_numpart_parallel', **parallel_runtimes)
    else:
        plot_results(args.results_file)
        