#!/usr/bin/env python

from os.path import dirname, abspath, join as pjoin

import pytest
import numpy as np

from Corrfunc.tests.common import Mr19_mock_northonly, Mr19_randoms_northonly
from Corrfunc.tests.common import (check_against_reference,
                                   check_vpf_against_reference)
from Corrfunc.tests.common import generate_isa_and_nthreads_combos
    

@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_DDrppi_mocks(Mr19_mock_northonly, isa, nthreads):
    from Corrfunc.mocks import DDrppi_mocks
    
    pimax = 40.0
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../mocks/tests/", "bins")
    autocorr = 1
    cosmology = 1
    
    ra,dec,cz,w = Mr19_mock_northonly
    results_DDrppi_mocks = DDrppi_mocks(autocorr, cosmology, nthreads,
                              pimax, binfile,
                              ra, dec, cz, weights1=w,
                              weight_type='pair_product',
                              output_rpavg=True, verbose=True,
                              isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../mocks/tests/", "Mr19_mock.DD")
    check_against_reference(results_DDrppi_mocks, file_ref, ravg_name='rpavg', ref_cols=(0, 4, 1))

    
@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_DDsmu_mocks(Mr19_randoms_northonly, isa, nthreads):
    from Corrfunc.mocks import DDsmu_mocks
    
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../mocks/tests/", "bins")
    autocorr = 1
    mu_max = 1.0
    nmu_bins = 10
    cosmology = 1
    
    ra,dec,cz,w = Mr19_randoms_northonly
    results_DDsmu_mocks = DDsmu_mocks(autocorr, cosmology, nthreads,
                              mu_max, nmu_bins, binfile,
                              ra, dec, cz, weights1=w,
                              weight_type='pair_product',
                              output_savg=True, verbose=True,
                              isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../mocks/tests/", "Mr19_mock_DDsmu.RR")
    check_against_reference(results_DDsmu_mocks, file_ref, ravg_name='savg', ref_cols=(0, 4, 1))
    
    
@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_DDtheta_mocks(Mr19_mock_northonly, isa, nthreads):
    from Corrfunc.mocks import DDtheta_mocks
    
    autocorr = 1
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../mocks/tests/", "angular_bins")
    
    ra,dec,cz,w = Mr19_mock_northonly
    results_DDtheta_mocks = DDtheta_mocks(autocorr, nthreads, binfile,
                              ra, dec,
                              weights1=w,
                              weight_type='pair_product',
                              output_thetaavg=True, fast_acos=False,
                              verbose=True, isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../mocks/tests/", "Mr19_mock_wtheta.DD")
    check_against_reference(results_DDtheta_mocks, file_ref, ravg_name='thetaavg', ref_cols=(0, 4, 1))


@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_vpf_mocks(Mr19_mock_northonly, isa, nthreads):
    from Corrfunc.mocks import vpf_mocks
    
    print("Beginning the VPF")
    # Max. sphere radius of 10 Mpc
    rmax = 10.0
    # 10 bins..so counts in spheres of radius 1, 2, 3, 4...10 Mpc spheres
    nbin = 10
    num_spheres = 10000
    num_pN = 6
    threshold_neighbors = 1  # does not matter since we have the centers
    centers_file = pjoin(dirname(abspath(__file__)),
                         "../../mocks/tests/data/",
                         "Mr19_centers_xyz_forVPF_rmax_10Mpc.txt")
    cosmology = 1
    
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../mocks/tests/", "angular_bins")
    
    ra,dec,cz,w = Mr19_mock_northonly
    results_vpf_mocks = vpf_mocks(rmax, nbin, num_spheres, num_pN,
                                    threshold_neighbors, centers_file,
                                    cosmology,
                                    ra, dec, cz, ra, dec, cz, 
                                    verbose=True, isa=isa,)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../mocks/tests/", "Mr19_mock_vpf")
    check_vpf_against_reference(results_vpf_mocks, file_ref)
    