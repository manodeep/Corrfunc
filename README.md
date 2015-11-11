/* Author: Manodeep Sinha <manodeep@gmail.com>
	 Date: At some point in 2015. 
	 LICENSE: MIT 
*/

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/manodeep/Corrfunc/master/LICENSE)
[![Coverity Scan](https://img.shields.io/coverity/scan/6982.svg)](https://scan.coverity.com/projects/manodeep-corrfunc)
[![Travis Build]https://travis-ci.org/manodeep/Corrfunc.svg?branch=python3)](https://travis-ci.org/manodeep/Corrfunc)

#Description

This repo contains a set of codes to measure the following OpenMP parallelized clustering 
measures in a cosmological box (co-moving XYZ) or on a mock (RA, DEC, CZ). Also, 
contains the associated paper to be published in Astronomy & Computing Journal (at some point). 

## Clustering Measures on a Cosmological box

All codes that work on cosmological boxes with co-moving positions are located in the 
``xi_theory`` directory. The various clustering measures are:

1. ``xi_of_r`` -- Measures auto/cross-correlations between two boxes. The boxes do not need to be cubes.

2. ``xi`` -- Measures 3-d auto-correlation in a cubic cosmological box. Assumes PERIODIC boundary conditions.

3. ``wp`` -- Measures auto 2-d point projected correlation function in a cubic cosmological box. Assumes PERIODIC boundary conditions. 

4. ``xi_rp_pi`` -- Measures the auto/cross correlation function between two boxes. The boxes do not need to be cubes. 

5. ``vpf`` -- Measures the void probability function + counts-in-cells. 

## Clustering measures on a Mock

All codes that work on mock catalogs (RA, DEC, CZ) are located in the ``xi_mocks`` directory. The
various clustering measures are:

1. ``DDrppi`` -- The standard auto/cross correlation between two data sets. The outputs, DD, DR and RR
can be combined using ``wprp`` to produce the Landy-Szalay estimator for $w_p(r_p)$. 

2. ``wtheta`` -- Computes angular correlation function between two data sets. The outputs from 
``DDtheta_mocks`` need to be combined with ``wtheta`` to get the full $\omega(\theta)$

3. ``vpf`` -- Computes the void probability function on mocks. 

# Science options

1. PERIODIC (ignored in case of wp/xi) -- switches PERIODIC boundary
conditions on/off. Enabled by default. 

2. OUTPUT_RPAVG -- switches on output of <rp> in each ``rp`` bin. Can be
a massive performance hit (~ 2.2x in case of wp). Disabled by default.
Needs code option DOUBLE_PREC to be enabled as well. For the mocks, 
OUTPUT_RPAVG causes only a mild increase in runtime and is enabled by 
default.

3. OUTPUT_THETAAVG -- switches on output of <theta> in each theta bin. 
Can be extremely slow (~5x) depending on compiler, and CPU capabilities. 
Disabled by default. 


## Mocks

1. LINK_IN_DEC -- creates binning in declination for mocks. Please check that for 
your desired binning in $r_p$/$\theta$, this binning does not produce incorrect 
results (due to numerical precision). 

2. LINK_IN_RA -- creates binning in RA once binning in DEC has been enabled. Same 
numerical issues as LINK_IN_DEC

3. FAST_DIVIDE --  Divisions are slow but required $DD(r_p,\pi)$. This Makefile
option (in mocks.options) replaces the divisions to a reciprocal followed by a 
Newton-Raphson. The code will run ~20% faster at the expense of some numerical precision. 
Please check that the loss of precision is not important for your use-case. Also, note 
that the mocks tests for $DD(r_p, \pi)$*will fail* if you enable FAST_DIVIDE. 

# Common Code options for both Mocks and Cosmological Boxes

1. DOUBLE_PREC -- does the calculations in double precision. Disabled
by default. 

2. USE_AVX -- uses the AVX instruction set found in Intel/AMD CPUs >= 2011
(Intel: Sandy Bridge or later; AMD: Bulldozer or later). Enabled by
default - code will not compile if the CPU does not support AVX instructions.
On Linux, check for "avx" in /proc/cpuinfo under flags. If you do not have
AVX, but have a SSE4 system instead, email me - I will send you a copy of
the code with SSE4 intrinsics. Or, take the relevant SSE code from the public repo at 
[pairwise](https://bitbucket.org/manodeep/pairwise).

3. USE_OMP -- uses OpenMP parallelization. Scaling is great for DD (perfect scaling
up to 12 threads in my tests) and okay (runtime becomes constant ~6-8 threads in
my tests) for DDrppi and wp. 


*Optimization for your architecture*

1. The values of bin_refine_factor and/or zbin_refine_factor in the countpairs_*.c
files control the cache-misses, and consequently, the runtime. In my trial-and-error
methods, I have seen any values larger than 3 are always slower. But some different
combination of 1/2 for (z)bin_refine_factor might be faster on your platform. 

2. If you have AVX2/AVX-512/KNC, you will need to rewrite the entire AVX section.

# Author

Pairwise is written/maintained by Manodeep Sinha. Please contact the [author](mailto:manodeep@gmail.com) in
case of any issues.

# LICENSE

Corrfunc is released under the MIT license. Basically, do what you want
with the code including using it in commercial application.

# Project URL
 
* version control (https://bitbucket.org/manodeep/corrfunc)
