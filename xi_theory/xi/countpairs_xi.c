/* File: countpairs.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "gridlink.h"//function proto-type for gridlink
#include "countpairs_xi.h" //function proto-type
#include "cellarray.h" //definition of struct cellarray
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar
#include "sglib.h"

#ifdef USE_AVX
#include "avx_calls.h"
#endif

#define BLOCK_SIZE     NVEC

#ifdef USE_OMP
#include <omp.h>
#endif

#ifndef PERIODIC
#warning "xi is only valid for PERIODIC boundary conditions. Ignoring the Makefile (non)-definition of PERIODIC"
#endif


void free_results_xi(results_countpairs_xi **results) 
{
	if(results == NULL)
		return;
	if(*results == NULL)
		return;

	results_countpairs_xi *tmp = *results;

	free(tmp->rupp);
	free(tmp->xi);
	free(tmp->npairs);
	free(tmp->rpavg);
	free(tmp);
	tmp = NULL;
}

results_countpairs_xi *countpairs_xi(const int64_t ND1, DOUBLE * restrict X1, DOUBLE * restrict Y1, DOUBLE * restrict Z1,
									 const double boxsize, 
#ifdef USE_OMP
									 const int numthreads,
#endif
									 const char *binfile)
{
	
  int bin_refine_factor=2;
	
  /***********************
   *initializing the  bins
   ************************/
  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin >= 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins must be > 0");
  

  /*---Create 3-D lattice--------------------------------------*/
  int nmesh_x=0,nmesh_y=0,nmesh_z=0;
  const DOUBLE xmin = 0.0, xmax=boxsize;
  const DOUBLE ymin = 0.0, ymax=boxsize;
  const DOUBLE zmin = 0.0, zmax=boxsize;
		
  cellarray *lattice = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
  if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
	fprintf(stderr,"countpairs> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n");
	bin_refine_factor *=2;
	const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;  		
	for(int64_t i=0;i<totncells;i++) {
	  free(lattice[i].x);
	  free(lattice[i].y);
	  free(lattice[i].z);
	}
	free(lattice);
	lattice = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
  }
  const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;  		
#ifdef USE_OMP
  omp_set_num_threads(numthreads);
#pragma omp parallel for schedule(dynamic) 
#endif
	for(int64_t icell=0;icell<totncells;icell++) {
		const cellarray *first=&(lattice[icell]);
		DOUBLE *x = first->x;
		DOUBLE *y = first->y;
		DOUBLE *z = first->z;
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,x,i,j); \
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,y,i,j);										\
			SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,z,i,j) }
		
		SGLIB_ARRAY_QUICK_SORT(DOUBLE, z, first->nelements, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
	}

	const DOUBLE side=boxsize;
  
#ifndef USE_OMP
	uint64_t npairs[nrpbin];	
	for(int i=0; i < nrpbin;i++) npairs[i] = 0;
#else
	uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, nrpbin);
#endif

  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
    rupp_sqr[i] = rupp[i]*rupp[i];
	}

#ifdef OUTPUT_RPAVG  
#ifndef USE_OMP	
  DOUBLE rpavg[nrpbin];
  for(int i=0; i < nrpbin;i++) {
    rpavg[i] = 0.0;
  }
#else	
	DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE),numthreads,nrpbin);
#endif//USE_OMP
#endif//OUTPUT_RPAVG

	const DOUBLE sqr_rpmax=rupp_sqr[nrpbin-1];
  const DOUBLE sqr_rpmin=rupp_sqr[0];

#ifdef USE_AVX
  AVX_FLOATS m_rupp_sqr[nrpbin];
  for(int i=0;i<nrpbin;i++) {
    m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
  }
#ifdef OUTPUT_RPAVG
  AVX_FLOATS m_kbin[nrpbin];
  for(int i=0;i<nrpbin;i++) {
    m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
  }
#endif//RPAVG  
#endif//AVX


  int interrupted=0;
  int64_t numdone=0;
  init_my_progressbar(totncells,&interrupted);

  /*---Loop-over-Data1-particles--------------------*/
#ifdef USE_OMP
#pragma omp parallel shared(numdone) 
  {
		int tid = omp_get_thread_num();
		uint64_t npairs[nrpbin];
		for(int i=0;i<nrpbin;i++) npairs[i] = 0;
#ifdef OUTPUT_RPAVG		
		DOUBLE rpavg[nrpbin];
		for(int i=0; i < nrpbin;i++) {
			rpavg[i] = 0.0;
		}
#endif		


#pragma omp for  schedule(dynamic) 
#endif
		for(int64_t index1=0;index1<totncells;index1++) {

#ifdef USE_OMP
		  if (omp_get_thread_num() == 0)
#endif
			my_progressbar(numdone,&interrupted);


#ifdef USE_OMP
		  #pragma omp atomic
#endif
		  numdone++;


		  cellarray *first = &(lattice[index1]);
		  const int iz = index1 % nmesh_z ;
		  const int ix = index1 / (nmesh_z * nmesh_y) ;
		  const int iy = (index1 - iz - ix*nmesh_z*nmesh_y)/nmesh_z ;
		  assert( ((iz + nmesh_z*iy + nmesh_z*nmesh_y*ix) == index1) && "Index reconstruction is wrong");
		  for(int iix=-bin_refine_factor;iix<=bin_refine_factor;iix++){
				const int iiix=(ix+iix+nmesh_x)%nmesh_x;
				DOUBLE off_xwrap=0.0;
				if(ix + iix >= nmesh_x) {
				  off_xwrap = -side;
				} else if (ix + iix < 0) {
				  off_xwrap = side;
				}

				
				for(int iiy=-bin_refine_factor;iiy<=bin_refine_factor;iiy++){
				  const int iiiy=(iy+iiy+nmesh_y)%nmesh_y;
				  DOUBLE off_ywrap = 0.0;
				  if(iy + iiy >= nmesh_y) {
						off_ywrap = -side;
				  } else if (iy + iiy < 0) {
						off_ywrap = side;
				  }
				  
				  for(int iiz=0;iiz<=bin_refine_factor;iiz++){
						const int iiiz=(iz+iiz+nmesh_z)%nmesh_z;
						DOUBLE off_zwrap = 0.0;
						if(iz + iiz >= nmesh_z) {
							off_zwrap = -side;
						} else if (iz + iiz < 0) {
							off_zwrap = side;
						}

						assert(iiix >= 0 && iiix < nmesh_x && iiiy >= 0 && iiiy < nmesh_y && iiiz >= 0 && iiiz < nmesh_z && "Checking that the second pointer is in range");
						const int64_t index2 = iiix*nmesh_y*nmesh_z + iiiy*nmesh_z + iiiz;
						const cellarray * second = &(lattice[index2]);
						const DOUBLE *x1 = first->x;
						const DOUBLE *y1 = first->y;
						const DOUBLE *z1 = first->z;
					
						const DOUBLE *x2 = second->x;
						const DOUBLE *y2 = second->y;
						const DOUBLE *z2 = second->z;

						for(int64_t i=0;i<first->nelements;i++) {
						  const int64_t start = (index1==index2) ? (i+1):0;					
							DOUBLE x1pos=x1[i] + off_xwrap;
							DOUBLE y1pos=y1[i] + off_ywrap;;
							DOUBLE z1pos=z1[i] + off_zwrap;
					  
#ifndef USE_AVX
							for(int64_t j=start;j<second->nelements;j+=BLOCK_SIZE) {
								int block_size=second->nelements - j;
								if(block_size > BLOCK_SIZE) block_size=BLOCK_SIZE;
						
								for(int jj=0;jj<block_size;jj++) {
									const DOUBLE dx = x1pos - x2[j+jj];
									const DOUBLE dy = y1pos - y2[j+jj];
									const DOUBLE dz = z1pos - z2[j+jj];
									const DOUBLE r2 = (dx*dx + dy*dy + dz*dz);
									const DOUBLE sqr_dz = dz*dz;
									if(dz < 0) {
										continue;
									}	else if(sqr_dz >= sqr_rpmax) {
										j = second->nelements;
									  break;
									}
									if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
										continue;
									}
#ifdef OUTPUT_RPAVG
									const DOUBLE r = SQRT(r2);
#endif					
									for(int kbin=nrpbin-1;kbin>=1;kbin--){
										if(r2 >= rupp_sqr[kbin-1]) {
											npairs[kbin]++;
#ifdef OUTPUT_RPAVG						
											rpavg[kbin] += r;
#endif						
											break;
										}
									}//searching for kbin loop
								}
							}//end of j loop
					  
#else //beginning of AVX section
		
#ifdef OUTPUT_RPAVG
							union int8 {
								AVX_INTS m_ibin;
								int ibin[NVEC];
							};
							union int8 union_rpbin;
					  
							union float8{
								AVX_FLOATS m_Dperp;
								DOUBLE Dperp[NVEC];
							};
							union float8 union_mDperp;
#endif				
					  
							const AVX_FLOATS m_x1pos = AVX_SET_FLOAT(x1pos);
							const AVX_FLOATS m_y1pos = AVX_SET_FLOAT(y1pos);
							const AVX_FLOATS m_z1pos = AVX_SET_FLOAT(z1pos);
					  
							int64_t j;
							for(j=start;j<=(second->nelements-NVEC);j+=NVEC) {
							  //Load the x/y/z arrays (NVEC at a time)
								const AVX_FLOATS x2pos = AVX_LOAD_FLOATS_UNALIGNED(&x2[j]);
								const AVX_FLOATS y2pos = AVX_LOAD_FLOATS_UNALIGNED(&y2[j]);
								const AVX_FLOATS z2pos = AVX_LOAD_FLOATS_UNALIGNED(&z2[j]);
						
								//x1-x2
								const AVX_FLOATS m_xdiff = AVX_SUBTRACT_FLOATS(m_x1pos,x2pos);
								//y1-y2
								const AVX_FLOATS m_ydiff = AVX_SUBTRACT_FLOATS(m_y1pos,y2pos);
								//z1-z2
								const AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_z1pos,z2pos);
								
								//set constant := sqr_rpmax
								const AVX_FLOATS m_sqr_rpmax = AVX_SET_FLOAT(sqr_rpmax);
								//set constant := sqr_rpmin
								const AVX_FLOATS m_sqr_rpmin = AVX_SET_FLOAT(sqr_rpmin);

								//set constant m_zero to 0.0
								const AVX_FLOATS m_zero  = AVX_SET_FLOAT((DOUBLE) 0.0);

								//(x1-x2)^2
								const AVX_FLOATS m_xdiff_sqr = AVX_SQUARE_FLOAT(m_xdiff);

								//(y1-y2)^2
								const AVX_FLOATS m_ydiff_sqr = AVX_SQUARE_FLOAT(m_ydiff);

								//(z1-z2)^2
								const AVX_FLOATS m_zdiff_sqr = AVX_SQUARE_FLOAT(m_zdiff);

								//(x1-x2)^2 + (y1-y2)^2
								const AVX_FLOATS m_xydiff_sqr_sum = AVX_ADD_FLOATS(m_xdiff_sqr,m_ydiff_sqr);
						
								//r2 now will contain (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 
								AVX_FLOATS r2 = AVX_ADD_FLOATS(m_zdiff_sqr,m_xydiff_sqr_sum);
								AVX_FLOATS m_mask_left;
						
								{

								  //Check if dz^2 >= sqr_rpmax. If so, break. 
								  AVX_FLOATS m_mask_pimax = AVX_COMPARE_FLOATS(m_zdiff_sqr, m_sqr_rpmax,_CMP_LT_OS);
								  const int test = AVX_TEST_COMPARISON(m_mask_pimax);
								  if(test == 0) {
										//If the execution reaches here -> then none of the NVEC zdiff values
										//are smaller than sqr_rpmax. We can terminate the j-loop now.
										
										//set j so that the remainder loop does not run
										j = second->nelements;
										//break out of the j-loop
										break;
								  }

									m_mask_pimax = AVX_BITWISE_AND(AVX_COMPARE_FLOATS(m_zdiff,m_zero,_CMP_GE_OS),m_mask_pimax);
									
									//Set r2 with sqr_rpmax where the mask is false.
									r2 = AVX_BLEND_FLOATS_WITH_MASK(m_rupp_sqr[nrpbin-1],r2,m_mask_pimax);
									
								  //Check if any of the NVEC distances are less than sqr_rpmax
									m_mask_left = AVX_COMPARE_FLOATS(r2,m_sqr_rpmax,_CMP_LT_OS);
									//If all points are >= sqr_rpmax, continue with the j-loop
									if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
										continue;
									}

						  
									//Create a mask for the NVEC distances that fall within sqr_rpmin and sqr_rpmax (sqr_rpmin <= dist < sqr_rpmax)
									const AVX_FLOATS m_mask = AVX_BITWISE_AND(m_mask_left, AVX_COMPARE_FLOATS(r2, m_sqr_rpmin, _CMP_GE_OS));
									if(AVX_TEST_COMPARISON(m_mask) == 0) {
										continue;
									}
									
									//Update r2 such that all distances that do not satisfy sqr_rpmin <= r2 < sqr_rpmax, get set to sqr_rpmax
									r2 = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2, m_mask);

									//Update the mask that now only contains points that need to be added to the npairs histogram
									m_mask_left = AVX_COMPARE_FLOATS(r2,m_sqr_rpmax,_CMP_LT_OS);
								}

								{ 
						  
#ifdef OUTPUT_RPAVG					  
								  //first do the sqrt since r2 contains squared distances
								  union_mDperp.m_Dperp = AVX_SQRT_FLOAT(r2);
								  AVX_FLOATS m_rpbin = AVX_SET_FLOAT((DOUBLE) 0.0);
#endif 					  
						  
									/* AVX_FLOATS m_all_ones  = AVX_CAST_INT_TO_FLOAT(AVX_SET_INT(-1));//-1 is 0xFFFF... and the cast just reinterprets (i.e., the cast is a no-op) */

								  //Loop over the histogram bins backwards. Most pairs will fall into the outer bins -> more efficient to loop backwards
								  //Remember that rupp[kbin-1] contains upper limit of previous bin -> lower radial limit of kbin
									for(int kbin=nrpbin-1;kbin>=1;kbin--) {
									  //Create a mask of pairwise separations that are greater than the lower radial limit of this bin (kbin)
										const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_GE_OS);
										//Do a bitwise AND to get the mask for separations that fall into this bin 
										const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
										//Create the mask for the remainder. This comparison should be exclusive with the comparison used for the m1 variable.
										m_mask_left = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_LT_OS);
										/* m_mask_left = AVX_XOR_FLOATS(m1, m_all_ones);//XOR with 0xFFFF... gives the bins that are smaller than m_rupp_sqr[kbin] (and is faster than cmp_p(s/d) in theory) */

										//Check the mask 
										const int test2  = AVX_TEST_COMPARISON(m_bin_mask);
										
										//Do a pop-count to add the number of bits. This is somewhat wasteful, since 
										//only 4 bits are set in DOUBLE_PREC mode (8 bits in regular float) but we 
										//are adding up all 32 bits in the integer. However, in my massive amount of 
										//testing with all sorts of faster implementations of popcount and table lookups,
										//builtin hardware popcnt always outperformed everything else. Thanks to NSA
										//for requiring a hardware popcnt I suppose. 
										npairs[kbin] += AVX_BIT_COUNT_INT(test2);
										
										//Add the kbin variable (as float) into the m_rpbin variable. 
										//This would be so much better implemented in AVX2 with support for integers
#ifdef OUTPUT_RPAVG
										m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
#endif					  
										//Check if there are any more valid points left. Break out of the kbin histogram loop if none are left
										const int test3 = AVX_TEST_COMPARISON(m_mask_left);
										if(test3 == 0) {
											break;
										}
									}

									//Since the m_rpbin is an AVX float, I have to truncate to an int to get the bin numbers. 
									//Only required when OUTPUT_RPAVG is enabled (i.e., the next jj-loop with the pragma unroll is in effect)
#ifdef OUTPUT_RPAVG
									union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
#endif					
								}
						
								//All these ops can be avoided (and anything leading to these) if the CPU
								//supports AVX 512 mask_add operation
#ifdef OUTPUT_RPAVG

//protect the unroll pragma in case compiler is not icc.								
#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
								for(int jj=0;jj<NVEC;jj++) {
									const int kbin = union_rpbin.ibin[jj];
									const DOUBLE r = union_mDperp.Dperp[jj];
									rpavg[kbin] += r;
								}
#endif//OUTPUT_RPAVG				  
							}//end of j-loop with AVX intrinsics
		  
							//Now take care of the remainder. 
							for(;j<second->nelements;j++) {
								const DOUBLE dz = z1pos - z2[j];
								const DOUBLE sqr_dz = dz*dz;
								if(sqr_dz >= sqr_rpmax ) break;
								const DOUBLE dx = x1pos - x2[j];
								const DOUBLE dy = y1pos - y2[j];
						
								const DOUBLE r2 = (dx*dx + dy*dy + dz*dz);
								if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
									continue;
								}
#ifdef OUTPUT_RPAVG
								const DOUBLE r = SQRT(r2);
#endif				  
								for(int kbin=nrpbin-1;kbin>=1;kbin--){
									if(r2 >= rupp_sqr[kbin-1]) {
										npairs[kbin]++;
#ifdef OUTPUT_RPAVG
										rpavg[kbin] += r;
#endif					  
										break;
									}
								}//searching for kbin loop
							}//end of remainder j loop
					  
#endif//end of AVX section
					  
						}//end of i loop
				  }//iiz loop over bin_refine_factor
				}//iiy loop over bin_refine_factor
			}//iix loop over bin_refine_factor
			  
		}//index1 loop over totncells
#ifdef USE_OMP
		for(int j=0;j<nrpbin;j++) {
		  all_npairs[tid][j] = npairs[j];
		}
#ifdef OUTPUT_RPAVG
		for(int j=0;j<nrpbin;j++) {
			all_rpavg[tid][j] = rpavg[j];
		}
#endif

	}//close the omp parallel region
#endif
  finish_myprogressbar(&interrupted);
  

#ifdef USE_OMP
  uint64_t npairs[nrpbin];
  for(int i=0;i<nrpbin;i++) npairs[i] = 0;
	
  for(int i=0;i<numthreads;i++) {
		for(int j=0;j<nrpbin;j++) {
			npairs[j] += all_npairs[i][j];
		}
  }
#ifdef OUTPUT_RPAVG
	DOUBLE rpavg[nrpbin];
	for(int i=0;i<nrpbin;i++) rpavg[i] = 0.0;

  for(int i=0;i<numthreads;i++) {
		for(int j=0;j<nrpbin;j++) {
			rpavg[j] += all_rpavg[i][j];
		}
	}
#endif

#endif


#ifdef OUTPUT_RPAVG  
  for(int i=0;i<nrpbin;i++) {
    if(npairs[i] > 0) {
      rpavg[i] /= (DOUBLE) npairs[i] ;
    }
  }
#endif

	//Pack in the results
	results_countpairs_xi *results = my_malloc(sizeof(*results), 1);
	results->nbin = nrpbin;
	results->npairs = my_malloc(sizeof(uint64_t), nrpbin);
	results->xi     = my_malloc(sizeof(DOUBLE)  , nrpbin);
	results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
	results->rpavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);

    const DOUBLE avgweight2 = 1.0, avgweight1 = 1.0;
    const DOUBLE density=0.5*avgweight2*ND1/(boxsize*boxsize*boxsize);//pairs are not double-counted                                                                                                                                           
    DOUBLE rlow=0.0 ;
    DOUBLE prefac_density_DD=avgweight1*ND1*density;

	for(int i=0;i<nrpbin;i++) {
		results->npairs[i] = npairs[i];
		results->rupp[i] = rupp[i];
#ifdef OUTPUT_RPAVG
		results->rpavg[i] = rpavg[i];
#else
		results->rpavg[i] = 0.0;
#endif

        const DOUBLE weight0 = (DOUBLE) results->npairs[i];
        /* compute xi, dividing summed weight by that expected for a random set */
        const DOUBLE vol=M_PI*(results->rupp[i]*results->rupp[i]*results->rupp[i]-rlow*rlow*rlow);
        const DOUBLE weightrandom = prefac_density_DD*vol;
        results->xi[i] = (weight0/weightrandom-1);
	}

  for(int64_t i=0;i<totncells;i++) {
		free(lattice[i].x);
		free(lattice[i].y);
		free(lattice[i].z);

  }
  
  free(lattice);
  free(rupp);
	
#ifdef USE_OMP
  matrix_free((void **) all_npairs, numthreads);
#ifdef OUTPUT_RPAVG
	matrix_free((void **) all_rpavg, numthreads);
#endif
	
#endif

	return results;

}
