#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

#include "function_precision.h"
#include "cellarray.h" //definition of struct cellarray
#include "utils.h" //all of the utilities
#include "gridlink.h"//function proto-type for gridlink
#include "countpairs_wp.h" //function proto-type
#include "progressbar.h" //for the progressbar

#include "sglib.h"

#ifdef USE_AVX
#include "avx_calls.h"
#endif

#ifdef USE_OMP
#include <omp.h>
#endif


#ifndef PERIODIC
#warning "wp is only valid for PERIODIC boundary conditions. Ignoring the Makefile (non)-definition of PERIODIC"
#endif

void free_results_wp(results_countpairs_wp **results)
{
	if(results == NULL)
		return;

	if(*results == NULL)
		return;

	results_countpairs_wp *tmp = *results;
	free(tmp->npairs);
	free(tmp->rupp);
	free(tmp->wp);
	free(tmp->rpavg);
	free(tmp);
	tmp = NULL;
}


results_countpairs_wp *countpairs_wp(const int64_t ND1, DOUBLE * restrict X1, DOUBLE * restrict Y1, DOUBLE * restrict Z1,
																		 const double boxsize, 
#ifdef USE_OMP
																		 const int numthreads,
#endif
																		 const char *binfile,
																		 const double pimax)
{

	int bin_refine_factor=2,zbin_refine_factor=1;
	int nmesh_x, nmesh_y, nmesh_z;
	

#ifdef USE_OMP
	if(numthreads == 1) {
		bin_refine_factor=2;
		zbin_refine_factor=1;
	} else {
		//seems redundant - but I can not discount the possibility that some
		//combination of these refine factors will be faster on a different architecture.
		bin_refine_factor=2;
		zbin_refine_factor=1;
	}
#endif

  /***********************
   *initializing the  bins
   ************************/
	double *rupp;
	double rpmin,rpmax;
	int nbin;
	setup_bins(binfile,&rpmin,&rpmax,&nbin,&rupp);
	assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
	assert(nbin > 0 && "Number of rp bins is valid");

	uint64_t npair[nbin];
	for(int i=0;i<nbin;i++) npair[i] = 0;

#ifdef OUTPUT_RPAVG
	DOUBLE rpavg[nbin];
	for(int i=0;i<nbin;i++) rpavg[i] = 0.0;
#endif	
	

	const DOUBLE xmin = 0.0, xmax=boxsize;
	const DOUBLE ymin = 0.0, ymax=boxsize;
	const DOUBLE zmin = 0.0, zmax=boxsize;
	DOUBLE rupp_sqr[nbin];
  for(int i=0;i<nbin;i++) {
    rupp_sqr[i] = rupp[i]*rupp[i];
  }

	const DOUBLE sqr_rpmin = rupp_sqr[0];
	const DOUBLE sqr_rpmax = rupp_sqr[nbin-1];
	
  //set up the 3-d grid structure. Each element of the structure contains a
  //pointer to the cellarray structure that itself contains all the points
  cellarray *lattice = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
	if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
    fprintf(stderr,"countpairs_wp> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n");
		bin_refine_factor *=2;
		zbin_refine_factor *=2;
		const int64_t totncells = nmesh_x*nmesh_y*(int64_t) nmesh_z;
		for(int64_t i=0;i<totncells;i++) {
			free(lattice[i].x);
			free(lattice[i].y);
			free(lattice[i].z);
		}
		free(lattice);
		lattice = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
	}
  const int64_t totncells = nmesh_x*nmesh_y*(int64_t) nmesh_z;

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

	
#ifdef USE_AVX
  AVX_FLOATS m_rupp_sqr[nbin];
	for(int i=0;i<nbin;i++) {
    m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
	}
#ifdef OUTPUT_RPAVG
	AVX_FLOATS m_kbin[nbin];
	for(int i=0;i<nbin;i++) {
		m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
	}
#endif//RPAVG
#endif
	

#ifdef USE_OMP
  uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, nbin);
#ifdef OUTPUT_RPAVG
	DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE), numthreads, nbin);
#endif//OUTPUT_RPAVG

#else//USE_OMP
  uint64_t local_npair[nbin];
  for(int i=0;i<nbin;i++) {
		local_npair[i]=0;
  }
#ifdef OUTPUT_RPAVG
	DOUBLE local_rpavg[nbin];
  for(int i=0;i<nbin;i++) {
		local_rpavg[i]=0;
  }
#endif//OUTPUT_RPAVG
#endif// USE_OMP

  const DOUBLE side=boxsize;
  int interrupted=0;
  int64_t numdone=0;
  init_my_progressbar(totncells,&interrupted);
  
#ifdef USE_OMP
#pragma omp parallel shared(numdone)
  {
		const int tid = omp_get_thread_num();
		uint64_t local_npair[nbin];
		for(int i=0;i<nbin;i++) {
			local_npair[i]=0;
		}
#ifdef OUTPUT_RPAVG
		DOUBLE local_rpavg[nbin];
		for(int i=0;i<nbin;i++) {
			local_rpavg[i]=0.0;
		}
#endif


#pragma omp for schedule(dynamic) nowait 
#endif
		for(int index1=0;index1<totncells;index1++) {

#ifdef USE_OMP
		  if (omp_get_thread_num() == 0)
#endif
			my_progressbar(numdone,&interrupted);
		  
		  
#ifdef USE_OMP
#pragma omp atomic
#endif
		  numdone++;
		  

			const int iz = index1 % nmesh_z;
			const int ix = index1 / (nmesh_y * nmesh_z );
			const int iy = (index1 - iz - ix*nmesh_z*nmesh_y)/nmesh_z;
			assert( ((iz + nmesh_z*iy + nmesh_z*nmesh_y*ix) == index1) && "Index reconstruction is wrong");
	  
			//Loop in the neighbouring cells in X. 
			//iiix accounts for periodic boundary conditions
			for(int iix=-bin_refine_factor;iix<=bin_refine_factor;iix++) {
				const int iiix=(ix+iix+nmesh_x)%nmesh_x;
				DOUBLE off_xwrap=0.0,off_ywrap=0.0,off_zwrap=0.0;
				off_xwrap=0.0;
				if(ix + iix < 0) {
					off_xwrap = side;
				} else if(ix + iix >= nmesh_x){
					off_xwrap = -side;
				}
				
				//Loop in the neighbouring cells in Y. 
				//iiiy accounts for periodic boundary conditions
				for(int iiy=-bin_refine_factor;iiy<=bin_refine_factor;iiy++) {
					const int iiiy=(iy+iiy+nmesh_y)%nmesh_y;
					off_ywrap=0.0;
					if(iy + iiy < 0) {
						off_ywrap = side;
					} else if(iy + iiy >= nmesh_y){
						off_ywrap = -side;
					}

					//Loop in the neighbouring cells in Z. 
					//iiiz accounts for periodic boundary conditions
					//However, since we are always doing an auto-correlation
					//we can avoid double-counting by only looping forward
					//in Z -> i.e., iiz loop starts from 0 rather than -zbin_refine_factor.
					for(int iiz=0;iiz<=zbin_refine_factor;iiz++) {
						const int iiiz=(iz+iiz+nmesh_z)%nmesh_z ;
						off_zwrap = 0.0;
						if(iz + iiz >= nmesh_z){
							off_zwrap = -side;
						}
			
						const int index2 = iiix*nmesh_y*nmesh_z + iiiy*nmesh_z + iiiz;
						assert(index2 >=0 && index2 < totncells && "index is within limits");
						const cellarray *second = &(lattice[index2]);
			
						const DOUBLE *x2 = second->x;
						const DOUBLE *y2 = second->y;
						const DOUBLE *z2 = second->z;

						const cellarray *first=&(lattice[index1]);
						const DOUBLE *x = first->x;
						const DOUBLE *y = first->y;
						const DOUBLE *z = first->z;
						
						for(int64_t i=0;i<first->nelements;i++) {
							const DOUBLE x1pos = x[i] + off_xwrap;
							const DOUBLE y1pos = y[i] + off_ywrap;
							const DOUBLE z1pos = z[i] + off_zwrap;

#ifndef USE_AVX							
							for(int64_t j=0;j<second->nelements;j++) {
								const DOUBLE dx = x2[j]-x1pos;
								const DOUBLE dy = y2[j]-y1pos;
								const DOUBLE dz = z2[j]-z1pos;
								if(dz < 0 ) {
									continue;
								} else if(dz >= pimax) {
									break;
								}

								const DOUBLE r2 = dx*dx + dy*dy;
								if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
									continue;
								}

#ifdef OUTPUT_RPAVG
								const DOUBLE r = SQRT(r2);
#endif
								for(int kbin=nbin-1;kbin>=1;kbin--) {
									if(r2 >= rupp_sqr[kbin-1]) {
										local_npair[kbin]++;
#ifdef OUTPUT_RPAVG
										local_rpavg[kbin]+=r;
#endif
										break;
									}
								}
								
							}
#else //AVX section
							const AVX_FLOATS m_xpos    = AVX_SET_FLOAT(x1pos);
							const AVX_FLOATS m_ypos    = AVX_SET_FLOAT(y1pos);
							const AVX_FLOATS m_zpos    = AVX_SET_FLOAT(z1pos);

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
							
							int64_t j;
							for(j=0;j<= (second->nelements - NVEC);j+=NVEC) {
								const AVX_FLOATS m_pimax = AVX_SET_FLOAT(pimax);
								const AVX_FLOATS m_zero  = AVX_SET_FLOAT((DOUBLE) 0.0);
					
								const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(&x2[j]);
								const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(&y2[j]);
								const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(&z2[j]);
					
								const AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_z2,m_zpos);//z2[j:j+NVEC-1] - z1
								const AVX_FLOATS m_xdiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_xpos,m_x2));//(x0 - x[j])^2
								const AVX_FLOATS m_ydiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_ypos,m_y2));//(y0 - y[j])^2
								AVX_FLOATS m_dist  = AVX_ADD_FLOATS(m_xdiff,m_ydiff);
					
								AVX_FLOATS m_mask_left;
					
								//Do all the distance cuts using masks here in new scope
								{
								  //the z2 arrays are sorted in increasing order. which means 
								  //the z2 value will increase in any future iteration of j. 
								  //that implies the zdiff values are also monotonically increasing
								  //Therefore, if none of the zdiff values are less than pimax, then 
								  //no future iteration in j can produce a zdiff value less than pimax.
								  //The code terminates the j-loop early in that case (and also sets 
								  //j equal to second->nelements to ensure that the remainder loop
								  //does not run either.
								  AVX_FLOATS m_mask_pimax = AVX_COMPARE_FLOATS(m_zdiff,m_pimax,_CMP_LT_OS);
									const int test = AVX_TEST_COMPARISON(m_mask_pimax);
									if(test == 0) {
									  //If the execution reaches here -> then none of the NVEC zdiff values 
									  //are smaller than pimax. We can terminate the j-loop now. 

									  //set j so that the remainder loop does not run
									  j = second->nelements;
									  //break out of the j-loop
									  break;
									}
									
									//Create a mask with true bits when  0 <= zdiff < pimax.
									m_mask_pimax = AVX_BITWISE_AND(AVX_COMPARE_FLOATS(m_zdiff,m_zero,_CMP_GE_OS),m_mask_pimax);
									
									//Set m_dist with sqr_rpmax where the mask is false. 
									m_dist = AVX_BLEND_FLOATS_WITH_MASK(m_rupp_sqr[nbin-1],m_dist,m_mask_pimax);		  

									//Create a mask for m_dist >= sqr_rpmin
									const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[0],_CMP_GE_OS);

									//Create a mask for m_dist < sqr_rpmax
									m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[nbin-1],_CMP_LT_OS);//will get utilized in the next section

									//Create a combined mask by bitwise and of m1 and m_mask_left. 
									//This gives us the mask for all sqr_rpmin <= m_dist < sqr_rpmax
									const AVX_FLOATS m_mask = AVX_BITWISE_AND(m1,m_mask_left);
									
									//Now check if any pair separations are within range
									const int test1 = AVX_TEST_COMPARISON(m_mask);
									
									//If not, continue with the next iteration of j-loop
									if(test1 == 0) {
										continue;
									}
								}
					
#ifdef OUTPUT_RPAVG
								union_mDperp.m_Dperp = AVX_SQRT_FLOAT(m_dist);
								AVX_FLOATS m_rpbin = AVX_SET_FLOAT((DOUBLE) 0.0);
#endif
								
								//Loop backwards through nbins. m_mask_left contains all the points that are less than rpmax
								for(int kbin=nbin-1;kbin>=1;kbin--) {
									const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_GE_OS);
									const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
									const int test2  = AVX_TEST_COMPARISON(m_bin_mask);
									local_npair[kbin] += AVX_BIT_COUNT_INT(test2);
#ifdef OUTPUT_RPAVG
									m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
#endif
									m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_LT_OS);
									const int test3 = AVX_TEST_COMPARISON(m_mask_left);
									if(test3 == 0)
										break;
								}

#ifdef OUTPUT_RPAVG
								union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
								//protect the unroll pragma in case compiler is not icc.
#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
								for(int jj=0;jj<NVEC;jj++) {
									const int kbin = union_rpbin.ibin[jj];
									const DOUBLE r = union_mDperp.Dperp[jj];
									local_rpavg[kbin] += r;
								}
#endif//OUTPUT_RPAVG

							}
			  
							//Now take care of the rest 
							for(;j<second->nelements;j++){
								const DOUBLE dz = z2[j] - z1pos;
								const DOUBLE dx = x2[j] - x1pos;
								const DOUBLE dy = y2[j] - y1pos;
								
								if(dz < 0 ) {
									continue;
								} else if(dz >= pimax) {
									break;
								}
								const DOUBLE r2 = dx*dx + dy*dy;
#ifdef OUTPUT_RPAVG
								const DOUBLE r = SQRT(r2);
#endif								
								if(r2 < sqr_rpmax && r2 >= sqr_rpmin) {
									for(int kbin=nbin-1;kbin>=1;kbin--) {
										if(r2 >= rupp_sqr[kbin-1] && r2 < rupp_sqr[kbin]) {
											local_npair[kbin]++;
#ifdef OUTPUT_RPAVG
											local_rpavg[kbin] += r;
#endif											
											break;
										}
									}
								}
							}//remainder loop over cellstruct second

#endif//end of AVX section
						}//loop over cellstruct first
					}//iiz loop over adjacent cells
				}//iiy loop over adjacent cells
			}//iix loop over adjacent cells
		}//index1 loop
#ifdef USE_OMP
		for(int j=0;j<nbin;j++) {
			all_npairs[tid][j] = local_npair[j];
		}
#ifdef OUTPUT_RPAVG
		for(int j=0;j<nbin;j++) {
			all_rpavg[tid][j] = local_rpavg[j];
		}
#endif		
  }//omp parallel
  finish_myprogressbar(&interrupted);	


  for(int i=0;i<numthreads;i++) {
		for(int j=0;j<nbin;j++) {
			npair[j] += all_npairs[i][j];
		}
  }
  matrix_free((void **) all_npairs,numthreads);
#ifdef OUTPUT_RPAVG
  for(int i=0;i<numthreads;i++) {
		for(int j=0;j<nbin;j++) {
			rpavg[j] += all_rpavg[i][j];
		}
  }
	matrix_free((void **) all_rpavg, numthreads);
#endif //OUTPUT_RPAVG

	
#else//end of OMP
	for(int i=0;i<nbin;i++) {
		npair[i] = local_npair[i];
	}
#ifdef OUTPUT_RPAVG
	for(int i=0;i<nbin;i++) {
		rpavg[i] = local_rpavg[i];
	}
#endif//OUTPUT_RPAVG
#endif//NO OMP


#ifdef OUTPUT_RPAVG
	for(int i=0;i<nbin;i++) {
		if(npair[i] > 0) {
			rpavg[i] /= (DOUBLE) npair[i];
		}
	}
#endif
	
	for(int64_t i=0;i<totncells;i++) {
		free(lattice[i].x);
		free(lattice[i].y);
		free(lattice[i].z);
	}
	free(lattice);

	//Pack in the results
	results_countpairs_wp *results = my_malloc(sizeof(*results), 1);
	results->nbin  = nbin;
	results->pimax = pimax;
	results->npairs = my_malloc(sizeof(uint64_t), results->nbin);
	results->wp = my_malloc(sizeof(DOUBLE), results->nbin);
	results->rupp   = my_malloc(sizeof(DOUBLE)  , results->nbin);
	results->rpavg  = my_malloc(sizeof(DOUBLE)  , results->nbin);


	const DOUBLE avgweight2 = 1.0, avgweight1 = 1.0;
	const DOUBLE density=0.5*avgweight2*ND1/(boxsize*boxsize*boxsize);//pairs are not double-counted
	DOUBLE rlow=0.0 ;
	DOUBLE prefac_density_DD=avgweight1*ND1*density;
	DOUBLE twice_pimax = 2.0*pimax;

	for(int i=0;i<results->nbin;i++) {
		results->npairs[i] = npair[i];
		results->rupp[i] = rupp[i];
#ifdef OUTPUT_RPAVG
		results->rpavg[i] = rpavg[i];
#else
		results->rpavg[i] = 0.0;
#endif
		const DOUBLE weight0 = (DOUBLE) results->npairs[i];
		/* compute xi, dividing summed weight by that expected for a random set */
		const DOUBLE vol=M_PI*(results->rupp[i]*results->rupp[i]-rlow*rlow)*twice_pimax;
		const DOUBLE weightrandom = prefac_density_DD*vol;
		results->wp[i] = (weight0/weightrandom-1)*twice_pimax;
		rlow=results->rupp[i];
	}
	free(rupp);
	return results;
}




