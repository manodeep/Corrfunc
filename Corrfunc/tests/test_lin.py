#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
from os.path import dirname, abspath, join as pjoin

import numpy as np
import pytest

from Corrfunc.io import read_fastfood_catalog
from Corrfunc.tests.common import gals_Mr19
from Corrfunc.tests.common import maxthreads

def _allclose(a, b, ravg_name=None):
    toret = True
    for name in ['npairs','weightavg',ravg_name]:
        if name == 'npairs':
            ac = np.all(a[name] == b[name])
        else:
            ac = np.allclose(a[name], b[name])
        if not ac:
            print("Mis-match for {0}: (a, b) = {1} ".format(name, list(zip(a[name], b[name]))))
        toret &= ac

    return toret


def _test_bin_types(pair_counter, *args, ravg_name=None, bin_types = ['custom', 'lin'], **kwargs):
    res = (pair_counter(*args, bin_type=bin_type,
           weight_type="pair_product", **kwargs) for bin_type in bin_types)
    assert _allclose(*res, ravg_name=ravg_name)
    

@pytest.fixture
def nthreads():
    return maxthreads()


@pytest.fixture(scope='module', params=['gals_Mr19',])
def points(request, npts=10**4):
    if request.param == 'gals_Mr19':
        filename = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/data", "gals_Mr19.ff")
        data = read_fastfood_catalog(filename, need_weights=True)
        data = np.array(data)
        np.random.seed(42)
        i = np.random.choice(data.shape[-1], size=npts, replace=False)
        data = data[:,i]

        return dict(data=data, boxsize=420.)
    
    raise ValueError(kind)
    
@pytest.fixture(scope='module', params=['Mr19_mock_northonly',])
def points_mock(request, npts=10**4):
    if request.param == 'Mr19_mock_northonly':
        filename = pjoin(dirname(abspath(__file__)),
                    "../../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
        data = read_fastfood_catalog(filename, need_weights=True)
        data = np.array(data)
        np.random.seed(42)
        i = np.random.choice(data.shape[-1], size=npts, replace=False)
        data = data[:,i]

        return dict(data=data)
    
    raise ValueError(kind)


@pytest.fixture(scope='module')
def linbins():
    binfile = np.linspace(2.0, 20.1, 12)
    return binfile


@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_mocks_DDrppi(points_mock, linbins, isa, nthreads, autocorr=True, pimax = 40.0, cosmology=1):
    from Corrfunc.mocks import DDrppi_mocks
    ra,dec,cz,w = points_mock['data']
    
    _test_bin_types(DDrppi_mocks, autocorr, cosmology, nthreads,
                  pimax, linbins,
                  ra, dec, cz, weights1=w,
                  output_rpavg=True, verbose=True,
                  isa=isa, ravg_name='rpavg')
    
@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_mocks_DDsmu(points_mock, linbins, isa, nthreads, autocorr=True, periodic=True, cosmology=1, nmu_bins = 10, mu_max = 1.0):
    from Corrfunc.mocks import DDsmu_mocks
    ra,dec,cz,w = points_mock['data']
    
    _test_bin_types(DDsmu_mocks, autocorr, cosmology, nthreads,
                  mu_max, nmu_bins, linbins,
                  ra, dec, cz, weights1=w,
                  output_savg=True, verbose=True,
                  isa=isa, ravg_name='savg')


@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_mocks_DDtheta(points_mock, linbins, isa, nthreads, autocorr=True):
    from Corrfunc.mocks import DDtheta_mocks
    ra,dec,cz,w = points_mock['data']
    
    _test_bin_types(DDtheta_mocks, autocorr, nthreads, linbins,
                  ra, dec, RA2=ra, DEC2=dec,
                  weights1=w,
                  weights2=w,
                  output_thetaavg=True, fast_acos=False,
                  verbose=True, isa=isa, ravg_name='thetaavg')


@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_theory_DD(points, linbins, isa, nthreads, autocorr=True, periodic=True):
    from Corrfunc.theory import DD
    x,y,z,w = points['data']
    boxsize = points['boxsize']
    
    _test_bin_types(DD, autocorr, nthreads, linbins, x, y, z, weights1=w,
                    verbose=True, periodic=periodic, boxsize=boxsize, isa=isa,
                    output_ravg=True, ravg_name='ravg')
    

@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_theory_DDrppi(points, linbins, isa, nthreads, autocorr=True, periodic=True, pimax=40.):
    from Corrfunc.theory import DDrppi
    x,y,z,w = points['data']
    boxsize = points['boxsize']

    _test_bin_types(DDrppi, autocorr, nthreads, pimax,
                  linbins, x, y, z, weights1=w,
                  verbose=True, periodic=periodic,
                  boxsize=boxsize, isa=isa,
                  output_rpavg=True, ravg_name='rpavg')


@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_theory_DDsmu(points, linbins, isa, nthreads, autocorr=True, periodic=True,
                                      mu_max=0.5, nmu_bins=10):
    from Corrfunc.theory import DDsmu
    x,y,z,w = points['data']
    boxsize = points['boxsize']

    _test_bin_types(DDsmu, autocorr, nthreads, linbins,
                  mu_max, nmu_bins,
                  x, y, z, weights1=w,
                  verbose=True, periodic=periodic,
                  boxsize=boxsize, output_savg=True, isa=isa,
                  ravg_name='savg')

    
@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_theory_wp(points, linbins, isa, nthreads, autocorr=True, periodic=True, pimax=40.):
    from Corrfunc.theory import wp
    x,y,z,w = points['data']
    boxsize = points['boxsize']
    
    _test_bin_types(wp, boxsize, pimax, nthreads,
                  linbins, x, y, z, weights=w,
                  verbose=True, isa=isa,
                  output_rpavg=True, ravg_name='rpavg')

    
@pytest.mark.parametrize('isa', ['fallback','sse42','avx','avx512f'])
def test_linear_binning_theory_xi(points, linbins, isa, nthreads, autocorr=True, periodic=True):
    from Corrfunc.theory import xi
    x,y,z,w = points['data']
    boxsize = points['boxsize']
    
    _test_bin_types(xi, boxsize, nthreads, linbins,
                  x, y, z, weights=w,
                  verbose=True, isa=isa,
                  output_ravg=True, ravg_name='ravg')
