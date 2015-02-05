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
#include "countpairs.h" //function proto-type
#include "cellarray.h" //definition of struct cellarray
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar

#ifdef USE_AVX
#include "avx_calls.h"
#endif


#ifdef USE_OMP
#include <omp.h>
#endif

void free_results(results_countpairs **results) 
{
	if(results == NULL)
		return;
	if(*results == NULL)
		return;

	results_countpairs *tmp = *results;

	free(tmp->rupp);
	free(tmp->npairs);
	free(tmp->rpavg);
	free(tmp);
	tmp = NULL;
}

results_countpairs * countpairs(const int64_t ND1, const DOUBLE * const X1, const DOUBLE * const Y1, const DOUBLE  * const Z1,
								const int64_t ND2, const DOUBLE * const X2, const DOUBLE * const Y2, const DOUBLE  * const Z2,
#ifdef USE_OMP
								const int numthreads,
#endif
								const int autocorr,
								const char *binfile)
{
	
  int bin_refine_factor=1;
  if(autocorr==1) {
    bin_refine_factor=2;
  } else {
    bin_refine_factor=1;
  }
#ifdef USE_OMP
	if(numthreads > 1) {

		//I have written it this way to maintain
		//compatibility with the previous chunk.
		// For instance, DR calculations might be faster
		// for bin_refine=2...
		if(autocorr==1) {
			bin_refine_factor=1;//benchmarked -> gives super-scalar performance
		} else {
			bin_refine_factor=1;//not benchmarked -> so may be changed in the future
		}
	}
#endif


  /***********************
   *initializing the  bins
   ************************/
  double *rupp;
	int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins is valid");
  
	//Find the min/max of the data
	DOUBLE xmin,xmax,ymin,ymax,zmin,zmax;
	xmin=1e10;ymin=1e10;zmin=1e10;
	xmax=0.0;ymax=0.0;zmax=0.0;
	get_max_min(ND1, X1, Y1, Z1, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax);
	
	if(autocorr==0) {
#ifndef SILENT				
		fprintf(stderr,"ND1 = %12"PRId64" [xmin,ymin,zmin] = [%lf,%lf,%lf], [xmax,ymax,zmax] = [%lf,%lf,%lf]\n",ND1,xmin,ymin,zmin,xmax,ymax,zmax);
#endif		
		get_max_min(ND2, X2, Y2, Z2, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax);
#ifndef SILENT		
		fprintf(stderr,"ND2 = %12"PRId64" [xmin,ymin,zmin] = [%lf,%lf,%lf], [xmax,ymax,zmax] = [%lf,%lf,%lf]\n",ND2,xmin,ymin,zmin,xmax,ymax,zmax);
#endif		
	}

#ifndef SILENT	
	fprintf(stderr,"Running with [xmin,xmax] = %lf,%lf\n",xmin,xmax);
	fprintf(stderr,"Running with [ymin,ymax] = %lf,%lf\n",ymin,ymax);
	fprintf(stderr,"Running with [zmin,zmax] = %lf,%lf\n",zmin,zmax);
#endif    

  /*---Create 3-D lattice--------------------------------------*/
  int nmesh_x=0,nmesh_y=0,nmesh_z=0;
		
  cellarray *lattice1 = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
  if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
	fprintf(stderr,"countpairs> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n");
	bin_refine_factor *=2;
	int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;  		
	for(int64_t i=0;i<totncells;i++) {
	  free(lattice1[i].x);
	  free(lattice1[i].y);
	  free(lattice1[i].z);
	}
	free(lattice1);
	lattice1 = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
  }
  

  cellarray *lattice2 = NULL;
  if(autocorr==0) {
	int ngrid2_x=0,ngrid2_y=0,ngrid2_z=0;
	lattice2 = gridlink(ND2, X2, Y2, Z2, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &ngrid2_x, &ngrid2_y, &ngrid2_z);
	assert(nmesh_x == ngrid2_x && "Both lattices have the same number of X bins");
	assert(nmesh_y == ngrid2_y && "Both lattices have the same number of Y bins");
	assert(nmesh_z == ngrid2_z && "Both lattices have the same number of Z bins");
  } else {
	lattice2 = lattice1;
  }
#ifdef PERIODIC
	const DOUBLE xdiff = (xmax-xmin);
	const DOUBLE ydiff = (ymax-ymin);
	const DOUBLE zdiff = (zmax-zmin);
#endif	

#ifndef USE_OMP
	uint64_t npairs[nrpbin];	
	for(int i=0; i < nrpbin;i++) npairs[i] = 0;
#else
	omp_set_num_threads(numthreads);
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

  DOUBLE sqr_rpmax=rupp_sqr[nrpbin-1];
  DOUBLE sqr_rpmin=rupp_sqr[0];

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


  int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;  
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


		  cellarray *first = &(lattice1[index1]);
		  const int iz = index1 % nmesh_z ;
		  const int ix = index1 / (nmesh_z * nmesh_y) ;
		  const int iy = (index1 - iz - ix*nmesh_z*nmesh_y)/nmesh_z ;
		  assert( ((iz + nmesh_z*iy + nmesh_z*nmesh_y*ix) == index1) && "Index reconstruction is wrong");
		  for(int iix=-bin_refine_factor;iix<=bin_refine_factor;iix++){
				int iiix;
#ifdef PERIODIC
				DOUBLE off_xwrap=0.0;
				if(ix + iix >= nmesh_x) {
				  off_xwrap = -xdiff;
				} else if (ix + iix < 0) {
				  off_xwrap = xdiff;
				}
				iiix=(ix+iix+nmesh_x)%nmesh_x;
#else	  
				iiix = iix+ix;
				if(iiix < 0 || iiix >= nmesh_x) {
				  continue;
				}
#endif	  
				
				for(int iiy=-bin_refine_factor;iiy<=bin_refine_factor;iiy++){
				  int iiiy;
#ifdef PERIODIC
				  DOUBLE off_ywrap = 0.0;
				  if(iy + iiy >= nmesh_y) {
						off_ywrap = -ydiff;
				  } else if (iy + iiy < 0) {
						off_ywrap = ydiff;
				  }
				  iiiy=(iy+iiy+nmesh_y)%nmesh_y;
#else	  
				  iiiy = iiy+iy;
				  if(iiiy < 0 || iiiy >= nmesh_y) {
						continue;
				  }
#endif	  
				  
				  for(int iiz=-bin_refine_factor;iiz<=bin_refine_factor;iiz++){
						int iiiz;
#ifdef PERIODIC
						DOUBLE off_zwrap = 0.0;
						if(iz + iiz >= nmesh_z) {
							off_zwrap = -zdiff;
						} else if (iz + iiz < 0) {
							off_zwrap = zdiff;
						}
						iiiz=(iz+iiz+nmesh_z)%nmesh_z;
#else	  
						iiiz = iiz+iz;
						if(iiiz < 0 || iiiz >= nmesh_z) {
							continue;
						}
#endif	  
						assert(iiix >= 0 && iiix < nmesh_x && iiiy >= 0 && iiiy < nmesh_y && iiiz >= 0 && iiiz < nmesh_z && "Checking that the second pointer is in range");
						const int64_t index2 = iiix*nmesh_y*nmesh_z + iiiy*nmesh_z + iiiz;
						const cellarray * second = &(lattice2[index2]);
						const DOUBLE *x1 = first->x;
						const DOUBLE *y1 = first->y;
						const DOUBLE *z1 = first->z;
					
						const DOUBLE *x2 = second->x;
						const DOUBLE *y2 = second->y;
						const DOUBLE *z2 = second->z;
					
						for(int64_t i=0;i<first->nelements;i++) {
							DOUBLE x1pos=x1[i];
							DOUBLE y1pos=y1[i];
							DOUBLE z1pos=z1[i];
#ifdef PERIODIC
							x1pos += off_xwrap;
							y1pos += off_ywrap;
							z1pos += off_zwrap;
#endif		
					  
#ifndef USE_AVX
							for(int64_t j=0;j<second->nelements;j+=NVEC) {
								int block_size=second->nelements - j;
								if(block_size > NVEC) block_size=NVEC;
						
								for(int jj=0;jj<block_size;jj++) {
									const DOUBLE dx = x1pos - x2[j+jj];
									const DOUBLE dy = y1pos - y2[j+jj];
									const DOUBLE dz = z1pos - z2[j+jj];
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
							for(j=0;j<=(second->nelements-NVEC);j+=NVEC) {
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
										if(test3 == 0)
											break;
									}

									//Since the m_rpbin is an AVX float, I have to truncate to an int to get the bin numbers. 
									//Only required when OUTPUT_RPAVG is enabled (i.e., the next jj-loop with the pragma unroll is in effect)
#ifdef OUTPUT_RPAVG
									union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
						
/* 								//All these ops can be avoided (and anything leading to these) if the CPU */
/* 								//supports AVX 512 mask_add operation */

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
								const DOUBLE dx = x1pos - x2[j];
								const DOUBLE dy = y1pos - y2[j];
								const DOUBLE dz = z1pos - z2[j];
						
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
	results_countpairs *results = my_malloc(sizeof(*results), 1);
	results->nbin = nrpbin;
	results->npairs = my_malloc(sizeof(uint64_t), nrpbin);
	results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
	results->rpavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);

	for(int i=0;i<nrpbin;i++) {
		results->npairs[i] = npairs[i];
		results->rupp[i] = rupp[i];
#ifdef OUTPUT_RPAVG
		results->rpavg[i] = rpavg[i];
#else
		results->rpavg[i] = 0.0;
#endif
	}

  for(int64_t i=0;i<totncells;i++) {
		free(lattice1[i].x);
		free(lattice1[i].y);
		free(lattice1[i].z);
		if(autocorr==0) {
			free(lattice2[i].x);
			free(lattice2[i].y);
			free(lattice2[i].z);
		}
  }
  
  free(lattice1);
  if(autocorr==0) {
		free(lattice2);
  }
	free(rupp);
	
#ifdef USE_OMP
  matrix_free((void **) all_npairs, numthreads);
#ifdef OUTPUT_RPAVG
	matrix_free((void **) all_rpavg, numthreads);
#endif
	
#endif

	return results;

}
