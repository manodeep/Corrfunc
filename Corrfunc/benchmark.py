from __future__ import print_function
import numpy as np
try:
    import kdcount
except:
    kdcount = None

try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None

try:
    from sklearn.neighbors import KDTree
except ImportError:
    KDTree = None

try:
    from halotools.mock_observables import npairs_3d \
        as halotools_npairs
except ImportError:
    halotools_npairs = None

try:
    from periodic_kdtree import PeriodicCKDTree
except ImportError:
    PeriodicCKDTree = None
    
try:
    from numba.decorators import autojit
except ImportError:
    autojit = None
    
from Corrfunc.theory import DD
from Corrfunc.utils import read_catalog
import time


def kdcount_timing(data1, data2, rbins, period):
    tree1 = kdcount.KDTree(data1, boxsize=period).root
    tree2 = kdcount.KDTree(data2, boxsize=period).root
    return tree1.count(tree2, r=rbins)


def halotools_timing(data1, data2, rbins, period):
    return halotools_npairs(data1, data2, rbins, period)


def cKDTree_timing(data1, data2, rbins, period):
    tree1 = cKDTree(data1, boxsize=period)
    tree2 = cKDTree(data2, boxsize=period)
    return tree1.count_neighbors(tree2, rbins)


def KDTree_timing(data1, data2, rbins, period):
    tree1 = cKDTree(data1, boxsize=period)
    tree2 = cKDTree(data2, boxsize=period)

    return tree1.count_neighbors(tree2, rbins)


def corrfunc_timing(data1, data2, rbins, period):
    x1 = data1[:, 0]
    y1 = data1[:, 1]
    z1 = data1[:, 2]

    x2 = data2[:, 0]
    y2 = data2[:, 1]
    z2 = data2[:, 2]
    nthreads = 1
    autocorr = 0
    periodic = False
    if period is not None:
        periodic = True

    result, api_time = DD(autocorr, nthreads, rbins,
                          x1, y1, z1,
                          X2=x2, Y2=y2, Z2=z2,
                          periodic=periodic,
                          c_api_timer=True,
                          isa='fastest')

    return (result['npairs'], api_time)


def corrfunc_timing_sse42(data1, data2, rbins, period):
    x1 = data1[:, 0]
    y1 = data1[:, 1]
    z1 = data1[:, 2]

    x2 = data2[:, 0]
    y2 = data2[:, 1]
    z2 = data2[:, 2]
    nthreads = 1
    autocorr = 0
    periodic = False
    if period is not None:
        periodic = True

    result, api_time = DD(autocorr, nthreads, rbins,
                          x1, y1, z1,
                          X2=x2, Y2=y2, Z2=z2,
                          periodic=periodic,
                          c_api_timer=True,
                          isa='SSE42')

    return (result['npairs'], api_time)


def corrfunc_timing_fallback(data1, data2, rbins, period):
    x1 = data1[:, 0]
    y1 = data1[:, 1]
    z1 = data1[:, 2]

    x2 = data2[:, 0]
    y2 = data2[:, 1]
    z2 = data2[:, 2]
    nthreads = 1
    autocorr = 0
    periodic = False
    if period is not None:
        periodic = True

    result, api_time = DD(autocorr, nthreads, rbins,
                          x1, y1, z1,
                          X2=x2, Y2=y2, Z2=z2,
                          periodic=periodic,
                          c_api_timer=True,
                          isa='FALLBACK')

    return (result['npairs'], api_time)


def periodic_wrap(dx, period):
    if dx > 0.5*period:
        dx -= period
    elif dx < -0.5*period:
        dx += period
        
    return dx


def pairwise_python(data1, data2, rbins, period):
    """
    These are adapted from from JakeVDP's famous pairwise blog entry
    https://jakevdp.github.io/blog/2013/06/15/numba-vs-cython-take-2/
    """

    x1 = data1[:, 0]
    y1 = data1[:, 1]
    z1 = data1[:, 2]

    x2 = data2[:, 0]
    y2 = data2[:, 1]
    z2 = data2[:, 2]

    N1 = len(x1)
    N2 = len(x2)
    Nbins = len(rbins) - 1
    sqr_rbins = rbins*rbins
    sqr_rmax = sqr_rbins[-1]
    sqr_rmin = sqr_rbins[0]
    counts = np.zeros(Nbins, dtype=np.int32)
    for i in xrange(N1):
        for j in xrange(N2):
            tmp = 0.0
            diffx = x1[i] - x2[j]
            diffy = y1[i] - y2[j]
            diffz = z1[i] = z2[j]
            if period is not None:
                diffx = periodic_wrap(diffx, period)
                diffy = periodic_wrap(diffy, period)
                diffz = periodic_wrap(diffz, period)

            tmp = diffx*diffx + diffy*diffy + diffz*diffz
            if tmp >= sqr_rmax or tmp <= sqr_rmin:
                continue

            for ibin in range(Nbins-1, 0, -1):
                if tmp >= sqr_rbins[ibin] and tmp < sqr_rbins[ibin+1]:
                    counts[ibin] += 1
                    break

    return counts


def main():
    functions_to_test = []
    function_names = []
    if kdcount is not None:
        functions_to_test.append(kdcount_timing)
        function_names.append("kdcount")

    if cKDTree is not None:
        functions_to_test.append(cKDTree_timing)
        function_names.append("cKDTree")

    if KDTree is not None:
        functions_to_test.append(KDTree_timing)
        function_names.append("KDTree")
       
    if halotools_npairs is not None:
        functions_to_test.append(halotools_timing)
        function_names.append("halotools")
    
    functions_to_test.append(corrfunc_timing)
    function_names.append("Corrfunc(AVX)")

    functions_to_test.append(corrfunc_timing_sse42)
    function_names.append("Corrfunc(SSE42)")

    functions_to_test.append(corrfunc_timing_fallback)
    function_names.append("Corrfunc(fallback)")

    # functions_to_test.append(pairwise_python)
    # function_names.append("Python")

    # if autojit is not None:
    #     pairwise_numba = autojit(pairwise_python)
    #     functions_to_test.append(pairwise_numba)
    #     function_names.append("Numba")

    X, Y, Z = read_catalog()
    npts = [100, 500, 1e3, 1e4, 1e5, 5e5]
    npts = [int(i) for i in npts]
    npts[npts > len(X)] = len(X)
    print("## Npts            ", end="")
    for ifunc in function_names:
        print(" {0:16s}".format(ifunc), end="")

    print("")

    for n in npts:
        X_sample = np.random.choice(X, n)
        Y_sample = np.random.choice(Y, n)
        Z_sample = np.random.choice(Z, n)
        data1 = np.array([X_sample, Y_sample, Z_sample]).reshape(n, 3)
        data2 = data1
        rbins = np.logspace(0.1, np.log10(24.0), 10)

        period = None
        print("{0:7d}   ".format(n), end="")
        for ii, func_name in enumerate(function_names):
            t0 = time.time()
            results_and_time = functions_to_test[ii](data1, data2, rbins,
                                                     period)
            try:
                result = np.array([np.long(l[0]) for l in results_and_time])
                runtime = results_and_time[1]
            except (IndexError, TypeError):
                result = results_and_time
                runtime = None
                
            # print("For {1}: result = {0}".format(result, func_name),
            #       file=sys.stderr)
            if runtime is None:
                t1 = time.time()
                runtime = t1-t0
            print(" {0:16.3f}".format(runtime), end="")

        print("")

if __name__ == '__main__':
    main()
