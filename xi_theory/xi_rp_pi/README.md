/* Distributed under the MIT LICENSE
   Author: Manodeep Sinha
   Email : manodeep@gmail.com
   Year  : 2014
   Location: Aspen, CO, USA
   Created to run DD/DR/RR pair-counting for
   CMASS galaxies. Specifically for Alexie + Shun's
   project. 
*/

Code to compute xi(rp, pi) on a rectangular parallelopiped.

Inputs:
	file1    - string. filename containing X Y Z (co-moving, including redshift distortions) positions in Mpc/h
	format1  - char.   format for file1. 'a' denotes ascii, 'f' denotes fast-food (fortran binary) format.
	file2    - string. filename containing X Y Z (co-moving, including redshift distortions) positions in Mpc/h
	format2  - char.   format for file2. 'a' denotes ascii, 'f' denotes fast-food (fortran binary) format.
  bins     - string.  filename containing the bins in the [rmin, rmax] format (units of Mpc/h). 
	pimax    - float/double. line of sight integration distance in co-moving Mpc/h. Default 40 Mpc/h.
	Nthreads - Number of OpenMP threads to use (if compiled with OpenMP support). Scaling is poor so anything more
					 	 than 8 threads is pointless. 
	
Outputs: Pair counts in (rp,pi) bins. Printed to stdout. The columns are [Npairs rpavg log10(rpbin_outer_edge)  pibin] (written out by countpairs.c)

Makefile option:  DOUBLE_PREC computes in double-precision. 
Compile  options: Use -xhost and -opt-prefetch. 

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








