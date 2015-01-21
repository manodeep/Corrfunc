/* Author: Manodeep Sinha <manodeep@gmail.com>
	 Date: At some point in early 2015. 
	 LICENSE: MIT 
*/

This repo contains a set of codes to measure the following three OpenMP parallelized correlation measures in a cosmological box. Also, contains the associated paper to be published
in Astronomy & Computing Journal (hopefully). 

1. Measures auto/cross-correlations between two boxes. 
Can be used to make full 3-d correlation function --  $\xi(r)$
Codes are in xi_of_r sub-directory.

2. Measures the auto/cross-correlations between two boxes.
Can be used to make projected correlation function -- $\xi(r_p,\pi)$
Codes are in the xi_rp_pi sub-directory.

3. Measures the full 2-point projected auto-correlation function
in a periodic cosmological box. Codes are in the wp sub-directory.


*Science options*

1. PERIODIC (ignored in case of wp) -- switches PERIODIC boundary
conditions on/off. Enabled by default. 

2. OUTPUT_RPAVG -- switches on output of <rp> in each bin. Can be
a massive performance hit (~ 2.2x in case of wp). Disabled by default.
Needs code option DOUBLE_PREC to be enabled as well. 

*Code options*

1. DOUBLE_PREC -- does the calculations in double precision. Disabled
by default. 

2. USE_AVX -- uses the AVX instruction set found in Intel/AMD CPUs >= 2011
(Intel: Sandy Bridge or later; AMD: Bulldozer or later). Enabled by
default - code will not compile if the CPU does not support AVX instructions.
On Linux, check for "avx" in /proc/cpuinfo under flags. If you do not have
AVX, but have a SSE4 system instead, email me - I will send you a copy of
the code with SSE4 intrinsics. 

3. USE_OMP -- uses OpenMP parallelization. Scaling is great for DD (perfect scaling
up to 12 threads in my tests) and okay (runtime becomes constant ~6-8 threads in
my tests) for DDrppi and wp. 

*Optimization for your architecture*

1. The values of bin_refine_factor and/or zbin_refine_factor in the countpairs_*.c
files control the cache-misses, and consequently, the runtime. In my trial-and-error
methods, I have seen any values larger than 3 are always slower. But some different
combination of 1/2 for (z)bin_refine_factor might be faster on your platform. 

2. If you have AVX2/AVX-512/KNC, you will need to rewrite the entire AVX section.


