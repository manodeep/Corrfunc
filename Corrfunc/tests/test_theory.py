#!/usr/bin/env python

from os.path import dirname, abspath, join as pjoin

import pytest
import numpy as np

from Corrfunc.tests.common import gals_Mr19
from Corrfunc.tests.common import (check_against_reference,
                                   check_vpf_against_reference)
from Corrfunc.tests.common import (generate_isa_and_nthreads_combos,
                                   maxthreads)


@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_DD(gals_Mr19, isa, nthreads):
    from Corrfunc.theory import DD
    
    boxsize = 420.
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../theory/tests/", "bins")
    autocorr = 1
    periodic = 1
    
    x,y,z,w = gals_Mr19
    results_DD = DD(autocorr, nthreads, binfile, x, y, z,
                            weights1=w, weight_type='pair_product',
                            periodic=periodic, boxsize=boxsize,
                            output_ravg=True, verbose=True,
                            isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "Mr19_DD_periodic")
    check_against_reference(results_DD, file_ref,
                            ravg_name='ravg', ref_cols=(0, 4, 1))


@pytest.mark.parametrize('funcname', ['DD', 'DDrppi', 'DDsmu'])
def test_boxsize(gals_Mr19, funcname, isa='fastest', nthreads=maxthreads()):
    '''Test the non-cubic and periodic boxsize features
    '''
    import Corrfunc.theory
    func = getattr(Corrfunc.theory, funcname)

    boxsize = 420.
    binfile = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "bins")
    periodic_ref = pjoin(dirname(abspath(__file__)),
                         "../../theory/tests/",
                         "Mr19_{}_periodic".format(funcname))
    nonperiodic_ref = pjoin(dirname(abspath(__file__)),
                            "../../theory/tests/",
                            "Mr19_{}_nonperiodic".format(funcname))
    autocorr = 1
    periodic = 1
    pimax = 40.0
    mu_max = 0.5
    nmu_bins = 10

    x, y, z, w = gals_Mr19
    args = [autocorr, nthreads, binfile, x, y, z]
    kwargs = dict(weights1=w, weight_type='pair_product',
                  periodic=periodic,
                  verbose=True, isa=isa)

    if funcname == 'DDrppi':
        args.insert(2, pimax)
        kwargs['output_rpavg'] = True
        ravg_name = 'rpavg'
    elif funcname == 'DDsmu':
        kwargs['output_savg'] = True
        args[3:3] = (mu_max, nmu_bins)
        ravg_name = 'savg'
    else:
        kwargs['output_ravg'] = True
        ravg_name = 'ravg'

    # scalar periodic
    results_DD = func(*args, boxsize=boxsize, **kwargs)
    check_against_reference(results_DD, periodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))

    # 0-d array, periodic
    results_DD = func(*args, boxsize=np.array(boxsize), **kwargs)
    check_against_reference(results_DD, periodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))

    # 1-d array, len 1, periodic
    results_DD = func(*args, boxsize=np.array([boxsize]), **kwargs)
    check_against_reference(results_DD, periodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))

    # 1-d array, len 3, periodic
    results_DD = func(*args, boxsize=np.array([boxsize, boxsize, boxsize]),
                      **kwargs)
    check_against_reference(results_DD, periodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))

    # 3-tuple periodic cube
    results_DD = func(*args, boxsize=(boxsize, boxsize, boxsize), **kwargs)
    check_against_reference(results_DD, periodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))

    # scalar non-periodic
    results_DD = func(*args, boxsize=(-1., -1., -1.), **kwargs)
    check_against_reference(results_DD, nonperiodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))

    # non-periodic via oversized box
    results_DD = func(*args, boxsize=(2000., 2000., 2000.), **kwargs)
    check_against_reference(results_DD, nonperiodic_ref,
                            ravg_name=ravg_name, ref_cols=(0, 4, 1))


@pytest.mark.parametrize('N', [1, 2])
def test_narrow_extent(N, isa='fastest', nthreads=maxthreads()):
    '''Test that a narrow particle distribution in a large box
    does not throw an error.
    Regression test for
    https://github.com/manodeep/Corrfunc/pull/276#issuecomment-1164872508
    '''
    from Corrfunc.theory import DD
    boxsize = (3., 3., 3.)
    seed = 42
    autocorr = 1
    r_bins = [0.2, 0.6, 1.0]
    if N == 2:
        pos = np.array([[0., 0.],
                        [0., 0.],
                        [0., 0.5],
                        ])
    elif N == 1:
        pos = np.array([[0.1],
                        [0.2],
                        [0.3],
                        ])
    else:
        raise NotImplemented(N)

    results = DD(autocorr, nthreads, r_bins, pos[0], pos[1], pos[2],
                 boxsize=boxsize, periodic=True,
                 verbose=True)

    if N == 2:
        assert np.all(results['npairs'] == [2, 0])
    elif N == 1:
        assert np.all(results['npairs'] == [0, 0])


@pytest.mark.parametrize('autocorr', [0, 1], ids=['cross', 'auto'])
@pytest.mark.parametrize('binref', [1, 2, 3])
@pytest.mark.parametrize('min_sep_opt', [False, True])
@pytest.mark.parametrize('maxcells', [1, 2, 3])
def test_duplicate_cellpairs(autocorr, binref, min_sep_opt, maxcells,
                             isa='fastest', nthreads=maxthreads()):
    '''A test to implement Manodeep's example from
    https://github.com/manodeep/Corrfunc/pull/277#issuecomment-1190921894
    where a particle pair straddles the wrap.

    Also tests the large Rmax case, where Rmax is almost as large as L/2.
    '''
    from Corrfunc.theory import DD
    boxsize = 432.
    r_bins = np.array([0.01, 0.4]) * boxsize
    pos = np.array([[0.02, 0.98],
                    [0., 0.],
                    [0., 0.0],
                    ]) * boxsize

    results = DD(autocorr, nthreads, r_bins, pos[0], pos[1], pos[2],
                 X2=pos[0], Y2=pos[1], Z2=pos[2],
                 boxsize=boxsize, periodic=True, isa=isa,
                 verbose=True, max_cells_per_dim=maxcells,
                 xbin_refine_factor=binref, ybin_refine_factor=binref,
                 zbin_refine_factor=binref, enable_min_sep_opt=min_sep_opt,
                 )

    assert np.all(results['npairs'] == [2])

    r_bins = np.array([0.2, 0.3, 0.49]) * boxsize
    pos = np.array([[0., 0.],
                    [0., 0.],
                    [0., 0.48],
                    ]) * boxsize

    results = DD(autocorr, nthreads, r_bins, pos[0], pos[1], pos[2],
                 X2=pos[0], Y2=pos[1], Z2=pos[2],
                 boxsize=boxsize, periodic=True, isa=isa,
                 verbose=True, max_cells_per_dim=maxcells,
                 xbin_refine_factor=binref, ybin_refine_factor=binref,
                 zbin_refine_factor=binref, enable_min_sep_opt=min_sep_opt,
                 )

    assert np.all(results['npairs'] == [0, 2])


@pytest.mark.parametrize('autocorr', [0, 1], ids=['cross', 'auto'])
@pytest.mark.parametrize('binref', [1, 2, 3], ids=['ref1', 'ref2', 'ref3'])
@pytest.mark.parametrize('min_sep_opt', [False, True], ids=['noMSO', 'MSO'])
@pytest.mark.parametrize('maxcells', [1, 2, 3], ids=['max1', 'max2', 'max3'])
@pytest.mark.parametrize('boxsize', [123., (51., 75., 123.)],
                         ids=['iso', 'aniso'])
@pytest.mark.parametrize('funcname', ['DD', 'DDrppi', 'DDsmu'])
@pytest.mark.parametrize('periodic', [False, True], ids=['nowrap', 'wrap'])
def test_brute(autocorr, binref, min_sep_opt, maxcells, boxsize, funcname,
               periodic, isa='fastest', nthreads=maxthreads()):
    '''Generate two small point clouds near particles near (0,0,0)
    and (L/2,L/2,L/2) and compare against the brute-force answer.

    In the periodic case, use close to the max allowable Rmax,
    0.49*min(boxsize).

    Tests both isotropic and anisotropic boxes.
    '''
    import Corrfunc.theory

    np.random.seed(1234)
    npts = 100
    eps = 0.2  # as a fraction of boxsize
    boxsize = np.array(boxsize)
    if periodic:
        bins = np.linspace(0.01, 0.49*boxsize.min(), 20)
    else:
        # non-periodic has no upper limit on Rmax
        bins = np.linspace(0.01, 2*boxsize.max(), 20)
    pimax = np.floor(0.49*boxsize.min())
    mu_max = 0.5
    nmu_bins = 10
    func = getattr(Corrfunc.theory, funcname)

    # two clouds of width eps*boxsize
    pos = np.random.uniform(low=-eps, high=eps,
                            size=(npts, 3))*boxsize
    pos[npts//2:] += boxsize/2.  # second cloud is in center of box
    pos %= boxsize

    # Compute the pairwise distance between particles with broadcasting.
    # Broadcasting (npts,1,3) against (npts,3) yields (npts,npts,3),
    # which is the array of all npts^2 (x,y,z) differences.
    pdiff = np.abs(pos[:, np.newaxis] - pos)
    if periodic:
        mask = pdiff >= boxsize/2
        pdiff -= mask*boxsize

    args = [autocorr, nthreads, bins, pos[:, 0], pos[:, 1], pos[:, 2]]
    kwargs = dict(periodic=periodic, isa=isa, boxsize=boxsize,
                  X2=pos[:, 0], Y2=pos[:, 1], Z2=pos[:, 2],
                  max_cells_per_dim=maxcells,
                  xbin_refine_factor=binref, ybin_refine_factor=binref,
                  zbin_refine_factor=binref, enable_min_sep_opt=min_sep_opt,
                  )

    if funcname == 'DDrppi':
        # Compute rp^2 = dx^2 + dy^2, and pi = abs(dz)
        args.insert(2, pimax)
        sqr_rpdiff = (pdiff[:, :, :2]**2).sum(axis=-1).reshape(-1)
        pidiff = np.abs(pdiff[:, :, 2]).reshape(-1)
        pibins = np.linspace(0., pimax, int(pimax)+1)
        brutecounts, _, _ = np.histogram2d(sqr_rpdiff, pidiff,
                                           bins=(bins**2, pibins))
        brutecounts = brutecounts.reshape(-1)  # corrfunc flattened convention
    elif funcname == 'DDsmu':
        # Compute s^2 = dx^2 + dy^2 + dz^2, mu = |dz| / s
        args[3:3] = (mu_max, nmu_bins)
        sdiff = np.sqrt((pdiff**2).sum(axis=-1).reshape(-1))
        sdiff[sdiff == 0.] = np.inf  # don't divide by 0
        mu = np.abs(pdiff[:, :, 2]).reshape(-1) / sdiff
        mubins = np.linspace(0, mu_max, nmu_bins+1)
        brutecounts, _, _ = np.histogram2d(sdiff, mu,
                                           bins=(bins, mubins))
        brutecounts = brutecounts.reshape(-1)  # corrfunc flattened convention
    elif funcname == 'DD':
        # Compute dist^2 = x^2 + y^2 + z^2, flattening because we don't
        # care which dist came from which particle pair
        sqr_pdiff = (pdiff**2).sum(axis=-1).reshape(-1)
        brutecounts, _ = np.histogram(sqr_pdiff, bins=bins**2)
    else:
        raise NotImplementedError(funcname)

    # spot-check that we have non-zero counts
    assert np.any(brutecounts > 0)

    pos = pos.T  # for indexing convenience
    results = func(*args, **kwargs)

    assert np.all(results['npairs'] == brutecounts)



@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_DDrppi(gals_Mr19, isa, nthreads):
    from Corrfunc.theory import DDrppi
    
    boxsize = 420.
    pimax = 40.0
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../theory/tests/", "bins")
    autocorr = 1
    periodic = 1
    
    x,y,z,w = gals_Mr19
    results_DDrppi = DDrppi(autocorr, nthreads, pimax, binfile, x, y, z,
                            weights1=w, weight_type='pair_product',
                            periodic=periodic, boxsize=boxsize,
                            output_rpavg=True, verbose=True,
                            isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "Mr19_DDrppi_periodic")
    check_against_reference(results_DDrppi, file_ref, ravg_name='rpavg', ref_cols=(0, 4, 1))

    
@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_DDsmu(gals_Mr19, isa, nthreads):
    from Corrfunc.theory import DDsmu
    
    boxsize = 420.
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../theory/tests/", "bins")
    autocorr = 1
    periodic = 1
    mu_max = 0.5
    nmu_bins = 10
    
    x,y,z,w = gals_Mr19
    results_DDsmu = DDsmu(autocorr, nthreads, binfile,
                            mu_max, nmu_bins,
                            x, y, z,
                            weights1=w, weight_type='pair_product',
                            periodic=periodic, boxsize=boxsize,
                            output_savg=True, verbose=True,
                            isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "Mr19_DDsmu_periodic")
    check_against_reference(results_DDsmu, file_ref, ravg_name='savg', ref_cols=(0, 4, 1))
    
    
@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos(extra_isa=['AVX2']))
def test_wp(gals_Mr19, isa, nthreads):
    from Corrfunc.theory import wp
    
    boxsize = 420.
    pimax = 40.
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../theory/tests/", "bins")
    
    x,y,z,w = gals_Mr19
    results_wp = wp(boxsize, pimax, nthreads, binfile,
                            x, y, z,
                            weights=w, weight_type='pair_product',
                            output_rpavg=True, verbose=True,
                            isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "Mr19_wp")
    check_against_reference(results_wp, file_ref, ravg_name='rpavg', cf_name='wp',
                            ref_cols=(4, 5, 1, 0))
    

@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_xi(gals_Mr19, isa, nthreads):
    from Corrfunc.theory import xi
    
    boxsize = 420.
    binfile = pjoin(dirname(abspath(__file__)),
                     "../../theory/tests/", "bins")
    
    x,y,z,w = gals_Mr19
    results_xi = xi(boxsize, nthreads, binfile,
                            x, y, z,
                            weights=w, weight_type='pair_product',
                            output_ravg=True, verbose=True,
                            isa=isa)
    
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "Mr19_xi")
    check_against_reference(results_xi, file_ref, ravg_name='ravg', cf_name='xi', ref_cols=(4, 5, 1, 0))
    

@pytest.mark.parametrize('isa,nthreads', generate_isa_and_nthreads_combos())
def test_vpf(gals_Mr19, isa, nthreads):
    from Corrfunc.theory import vpf
    
    boxsize = 420.
    rmax = 10.0
    nbin = 10
    nspheres = 10000
    num_pN = 6
    seed = -1234
    periodic = 1
    
    x,y,z,w = gals_Mr19
    results_vpf  = vpf(rmax, nbin, nspheres, num_pN,
                          seed, x, y, z, verbose=True, periodic=periodic,
                          boxsize=boxsize)
    #results_vpf = results_vpf.view(dtype=np.float64).reshape(nbin,-1)  # flatten to same shape as results
    file_ref = pjoin(dirname(abspath(__file__)),
                    "../../theory/tests/", "Mr19_vpf_periodic")
    check_vpf_against_reference(results_vpf, file_ref)
