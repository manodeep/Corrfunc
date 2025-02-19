# Corrfunc forked:
I rewrite the utils.py, so that the `convert_rp_pi_counts_to_wp` and `convert_3d_counts_to_cf` works more properly. The formal call is still valid, and I added that if you provide weight in the pair-counting, it will account for the weighting.

```python
# load some catalogue...

# Code 1
wei_norm = galaxy['w'] / (galaxy['w'].mean())
wei_norm_r = random['w'] / (random['w'].mean())

dd = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                   RA1=galaxy['ra'], DEC1=galaxy['dec'], CZ1=galaxy['distance'], weights1=wei_norm, is_comoving_dist=True, weight_type='pair_product')
dr = DDrppi_mocks(
    autocorr=False, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin, 
    RA1=galaxy['ra'], DEC1=galaxy['dec'], CZ1=galaxy['distance'], weights1=wei_norm, 
    RA2=random['ra'], DEC2=random['dec'], CZ2=random['distance'], weights2=wei_norm_r, 
    is_comoving_dist=True, weight_type='pair_product')
rr = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                   RA1=random['ra'], DEC1=random['dec'], CZ1=random['distance'], weights1=wei_norm_r, is_comoving_dist=True, weight_type='pair_product')

Nd = len(galaxy)
Nr = len(random)

wp_1 = convert_rp_pi_counts_to_wp(Nd, Nd, Nr, Nr, dd, dr, dr, rr, pimax=pimax, nrpbins=Nbins)


# Code 2
dd = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                   RA1=galaxy['ra'], DEC1=galaxy['dec'], CZ1=galaxy['distance'], weights1=galaxy['w'], is_comoving_dist=True, weight_type='pair_product')
dr = DDrppi_mocks(
    autocorr=False, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin, 
    RA1=galaxy['ra'], DEC1=galaxy['dec'], CZ1=galaxy['distance'], weights1=galaxy['w'], 
    RA2=random['ra'], DEC2=random['dec'], CZ2=random['distance'], weights2=random['w'], 
    is_comoving_dist=True, weight_type='pair_product')
rr = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                   RA1=random['ra'], DEC1=random['dec'], CZ1=random['distance'], weights1=random['w'], is_comoving_dist=True, weight_type='pair_product')

Nd = galaxy['w'].sum()
Nr = random['w'].sum()

wp_2 = convert_rp_pi_counts_to_wp(Nd, Nd, Nr, Nr, dd, dr, dr, rr, pimax=pimax, nrpbins=Nbins)
assert np.isclose(wp_1, wp_2).all()
```
This code will work.

Note that, for simplicity, I didn't add new parameters to the function. 
Instead you can 
- normalize the weight first, e.g. `weight_normal = weight / weight.mean()`, and pass the parameter in the old way.
- pass the sum of weight of dataset1 to `ND1`, sum of weight of dataset2 to `ND2`, etc. This is reasonable, because if you assume no weighting(that is weight of every point is 1), then the sum of weights equals to number of points.

Both way will work. 