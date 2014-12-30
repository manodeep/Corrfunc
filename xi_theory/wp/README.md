/* Distributed under the MIT LICENSE
   Author: Manodeep Sinha
   Email : manodeep@gmail.com
   Year  : 2012
   Location: Nashville, TN
	 Created for the MCMC chains to compute wp
	 on a periodic, cosmological box
*/

Code to compute wp(rp) on a periodic, cosmological box. Ignores the Makefile option PERIODIC.  
 
Inputs:
	boxsize - double. Boxsize for the periodic cube that contains all the particles. 
	file    - string. filename containing X Y Z (co-moving, including redshift distortions) positions in Mpc/h
	format  - char.   format for file2. 'a' denotes ascii, 'f' denotes fast-food (fortran binary) format.
  bins    - string.  filename containing the bins in the [rmin, rmax] format (units of Mpc/h). 
	pimax   - float/double. line of sight integration distance in co-moving Mpc/h. Default 40 Mpc/h.
	Nthreads - Number of OpenMP threads to use, if compiled with OpenMP support. Scaling is poor so anything more
					 	 than 8 threads is pointless. Not required unless compiled with OpenMP. 
	
Outputs: Prints to stdout. The columns are [ wp <rp> rlow rupp npairs ] where <rp> is identically 0.0
unless the Makefile option OUTPUT_RPAVG is defined. rlow is the lower edge of the bin, while rupp is the
upper edge. npairs are the number of pairs in the bin (not double-counted). 


Makefile option:  

1. DOUBLE_PREC 			computes in double-precision. 
2. PERIODIC 				IS IGNORED. PERIODIC is *always* ON. 
3. OUTPUT_RPAVG			computes ravg explicitly. Disabled by default (the ravg term in the output file is identically 0.0)
4. USE_AVX  				enables the hand-vectorized AVX intrinsics - will only work if cpu >= Sandy Bridge
5. USE_OMP          Compile with openmp support - requires an additional nthreads parameter on the command-line. The 
                    openmp scaling is fairly poor - so use up to 4-6 threads. 

*Notes*

1) The binning in the pi-direction is currently fixed to 1 Mpc/h bins. This can be changed by modifying the
line in DDrppi.c that sets npibin = (int) pimax. 

2) Fast-food format is a fortran binary format that contains a header, followed by the xyz data.
The header contains an integer array of 5 elements (idat), floating array with 9 elements (fdat) and a float with
the redshift corresponding to that file. idat[1] contains the number of particles in the data file.
fdat[0] contains the box-size. This header is followed by sizeof(double/float)*Npart x, then y, then z.
Don't forget the padding bytes when reading in C (all of these are implemented in ftread.c). There is a
version of ftwrite.c floating around on the web (comes with Rockstar if nothing else).

*Possible runtime optimizations*

1) Change bin_refine_factor and/or zbin_refine_factor in countpairs_rp_pi.c and profile. 








