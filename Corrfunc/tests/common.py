from os.path import dirname, abspath, join as pjoin

import numpy as np
import pytest

from Corrfunc.io import read_fastfood_catalog

@pytest.fixture(scope='module')
def gals_Mr19():
    
    filename = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/data", "gals_Mr19.ff")
    x, y, z, w = read_fastfood_catalog(filename, need_weights=True)
    
    return x, y, z, w


@pytest.fixture(scope='module')
def Mr19_mock_northonly():
    filename = pjoin(dirname(abspath(__file__)),
                "../../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
    ra,dec,cz,w = read_fastfood_catalog(filename, need_weights=True)
    
    return ra, dec, cz, w


@pytest.fixture(scope='module')
def Mr19_randoms_northonly():
    filename = pjoin(dirname(abspath(__file__)),
                "../../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
    ra,dec,cz,w = read_fastfood_catalog(filename, need_weights=True)
    
    return ra, dec, cz, w


@pytest.fixture(scope='module')
def nthreads():
    '''Use as many threads as cores that are available to this process'''
    import os
    return len(os.sched_getaffinity(0))


def _check_against_reference(results, filename, results_cols=(-2, -1, 2), ref_cols=(0, 4, 1)):
    # results is output of Python function
    # filename the reference counts
    # ref_cols in order of npairs, weightavg, rpavg
    refs = np.loadtxt(filename, unpack=True, usecols=ref_cols)
    for ii, items in enumerate(results):
        for icol, ref in zip(results_cols, refs):
            assert np.allclose(items[icol], ref[ii])
