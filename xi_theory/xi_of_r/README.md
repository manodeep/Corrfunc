/* Distributed under the MIT LICENSE
   Author: Manodeep Sinha
   Email : manodeep@gmail.com
   Year  : 2014
   Location: Aspen, CO, USA

   Created to complement xi(rp, pi) DD/DR/RR pair-counting 
   for CMASS galaxies. Specifically for Alexie + Shun's
   project. 
*/

Code to compute xi(r) on a rectangular parallelopiped.

Inputs: 

* file1     - string. filename containing X Y Z (co-moving, including redshift distortions) positions in Mpc/h
* format1   - char.   format for file1. 'a' denotes ascii, 'f' denotes fast-food (fortran binary) format.
* file2     - string. filename containing X Y Z (co-moving, including redshift distortions) positions in Mpc/h
* format2   - char.   format for file2. 'a' denotes ascii, 'f' denotes fast-food (fortran binary) format.
* bins      - string.  filename containing the bins in the [rmin, rmax] format (units of Mpc/h). 
* nthreads  - integer. Only required if compiled with openmp support.

Typical command-line would be:

1. Auto-correlations : ./DD cmassmock_AbM_swot_scatter0.10dex_All_V_PEAK_Zspace.ff cmassmock_AbM_swot_scatter0.10dex_All_V_PEAK_Zspace.ff f bins [nthreads] > DD

2. Cross-correlations: ./DD cmassmock_AbM_swot_scatter0.10dex_All_V_PEAK_Zspace.ff f random_AbM_swot_scatter0.10dex_All_V_PEAK_Zspace.ff f bins [nthreads] > DR
	
Outputs: Pair counts in r bins. Printed to stdout. The columns are [Npairs ravg   rmin rmax] (written out by countpairs.c)

Makefile option:  

1. DOUBLE_PREC 			computes in double-precision. 
2. PERIODIC 				enables periodic boundary conditions
3. OUTPUT_RPAVG			computes ravg explicitly. Disabled by default (the ravg term in the output file is identically 0.0)
4. USE_AVX  				enables the hand-vectorized AVX intrinsics - will only work if cpu >= Sandy Bridge
5. USE_OMP          Compile with openmp support - requires an additional nthreads parameter on the command-line. The 
                    openmp scaling is fantastic up to ~ 12 threads for a compute-intense work (i.e., > 6-7s, AVX, PERIODIC).
										YMMV.


Compile  options: Use -xhost and -opt-prefetch -ipo  (for icc)


Notes:

1) Fast-food format is a fortran binary format that contains a header, followed by the xyz data.
The header contains an integer array of 5 elements (idat), floating array with 9 elements (fdat) and a float with
the redshift corresponding to that file. idat[1] contains the number of particles in the data file.
fdat[0] contains the box-size. This header is followed by sizeof(double/float)*Npart x, then y, then z.
Don't forget the padding bytes when reading in C (all of these are implemented in ftread.c). There is a
version of ftwrite.c floating around on the web (comes with Rockstar if nothing else).