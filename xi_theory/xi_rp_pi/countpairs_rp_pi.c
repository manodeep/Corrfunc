#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "countpairs_rp_pi.h" //function proto-type
#include "gridlink.h"//function proto-type for gridlink
#include "cellarray.h" //definition of struct cellarray
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar

#ifdef USE_AVX
#include "avx_calls.h"
#endif

#ifdef USE_OMP
#include <omp.h>
#endif


void free_results_rp_pi(results_countpairs_rp_pi **results)
{
	if(results==NULL)
		return;

	if(*results==NULL)
		return;
	
	results_countpairs_rp_pi *tmp = *results;

	free(tmp->npairs);
	free(tmp->rupp);
	free(tmp->rpavg);
	free(tmp);
	tmp = NULL;
}


results_countpairs_rp_pi * countpairs_rp_pi(const int64_t ND1, const DOUBLE *X1, const DOUBLE *Y1, const DOUBLE *Z1,
											const int64_t ND2, const DOUBLE *X2, const DOUBLE *Y2, const DOUBLE *Z2,
#ifdef USE_OMP
											const int numthreads,
#endif
											const int autocorr,
											const char *binfile,
											const double pimax)
{
  int bin_refine_factor=1;
	int zbin_refine_factor=2;
	if(autocorr==1) {
		bin_refine_factor=2;
		zbin_refine_factor=2;
	} else {
		bin_refine_factor=1;
		zbin_refine_factor=1;
	}
#ifdef USE_OMP
	if(numthreads > 1) {
		if(autocorr==1) {
			bin_refine_factor=2;
			zbin_refine_factor=2;
		} else {
			bin_refine_factor=1;
			zbin_refine_factor=1;
		}
	}
#endif
	const int npibin = (int) pimax;
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
		fprintf(stderr,"ND1 = %12"PRId64" [xmin,ymin,zmin] = [%lf,%lf,%lf], [xmax,ymax,zmax] = [%lf,%lf,%lf]\n",ND1,xmin,ymin,zmin,xmax,ymax,zmax);
		get_max_min(ND2, X2, Y2, Z2, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax);
		fprintf(stderr,"ND2 = %12"PRId64" [xmin,ymin,zmin] = [%lf,%lf,%lf], [xmax,ymax,zmax] = [%lf,%lf,%lf]\n",ND2,xmin,ymin,zmin,xmax,ymax,zmax);
	}
	fprintf(stderr,"Running with [xmin,xmax] = %lf,%lf\n",xmin,xmax);
	fprintf(stderr,"Running with [ymin,ymax] = %lf,%lf\n",ymin,ymax);
	fprintf(stderr,"Running with [zmin,zmax] = %lf,%lf\n",zmin,zmax);
	

	/*---Create 3-D lattice--------------------------------------*/
	int nmesh_x=0,nmesh_y=0,nmesh_z=0;

	cellarray *lattice1 = gridlink(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
	cellarray *lattice2 = NULL;
	if(autocorr==0) {
		int ngrid2_x=0,ngrid2_y=0,ngrid2_z=0;
		lattice2 = gridlink(ND2, X2, Y2, Z2, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &ngrid2_x, &ngrid2_y, &ngrid2_z);
		assert(nmesh_x == ngrid2_x && "Both lattices have the same number of X bins");
		assert(nmesh_y == ngrid2_y && "Both lattices have the same number of Y bins");
		assert(nmesh_z == ngrid2_z && "Both lattices have the same number of Z bins");
	} else {
		lattice2 = lattice1;
	}
	const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;

#ifdef PERIODIC
	const DOUBLE xdiff = (xmax-xmin);
	const DOUBLE ydiff = (ymax-ymin);
	const DOUBLE zdiff = (zmax-zmin);
#endif
	
	/* DOUBLE logrpmax,logrpmin,dlogrp; */
  const DOUBLE dpi = pimax/(DOUBLE)npibin ;
  const DOUBLE inv_dpi = 1.0/dpi;

  DOUBLE rupp_sqr[nrpbin];
	const int64_t totnbins = (npibin+1)*(nrpbin+1);
	for(int i=0; i < nrpbin;i++) {
		rupp_sqr[i] = rupp[i]*rupp[i];
	}

	const DOUBLE sqr_rpmax=rupp_sqr[nrpbin-1];
	const DOUBLE sqr_rpmin=rupp_sqr[0];
	
#ifndef USE_OMP
	uint64_t npairs[totnbins];
#ifdef OUTPUT_RPAVG
	DOUBLE rpavg[totnbins];
#endif	
	for(int ibin=0;ibin<totnbins;ibin++) {
		npairs[ibin]=0;
#ifdef OUTPUT_RPAVG		
		rpavg[ibin] = 0.0;
#endif		
	}
#else
	omp_set_num_threads(numthreads);
	uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, totnbins);
#ifdef OUTPUT_RPAVG
	DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE),numthreads,totnbins);
#endif	
#endif

	
#ifdef USE_AVX
  AVX_FLOATS m_rupp_sqr[nrpbin];
  AVX_FLOATS m_kbin[nrpbin];
  for(int i=0;i<nrpbin;i++) {
    m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
    m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
  }
#endif
  
  int interrupted=0;
  int64_t numdone=0;
  init_my_progressbar(totncells,&interrupted);

	
#ifdef USE_OMP
#pragma omp parallel shared(numdone)
	{
		const int tid = omp_get_thread_num();
		uint64_t npairs[totnbins];
		for(int i=0;i<totnbins;i++) npairs[i] = 0;
#ifdef OUTPUT_RPAVG		
		DOUBLE rpavg[totnbins];
		for(int i=0;i<totnbins;i++) rpavg[i] = 0.0;
#endif			

#pragma omp for  schedule(dynamic)
#endif
		/*---Loop-over-lattice1--------------------*/
		for(int index1=0;index1<totncells;index1++) {
		  
#ifdef USE_OMP
		  if (omp_get_thread_num() == 0)
#endif
			my_progressbar(numdone,&interrupted);
		  
		  
#ifdef USE_OMP
#pragma omp atomic
#endif
		  numdone++;
		  
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
					
					for(int iiz=-zbin_refine_factor;iiz<=zbin_refine_factor;iiz++){
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
						const cellarray *second = &(lattice2[index2]);
						const cellarray *first  = &(lattice1[index1]);
						
						const DOUBLE *x1 = first->x;
						const DOUBLE *y1 = first->y;
						const DOUBLE *z1 = first->z;
						
						const DOUBLE *x2 = second->x;
						const DOUBLE *y2 = second->y;
						const DOUBLE *z2 = second->z;
						
						for(int64_t i=0;i<first->nelements;i++){
							DOUBLE x1pos = x1[i];
							DOUBLE y1pos = y1[i];
							DOUBLE z1pos = z1[i];
#ifdef PERIODIC
							x1pos += off_xwrap;
							y1pos += off_ywrap;
							z1pos += off_zwrap;
#endif
							
#ifndef USE_AVX	//Beginning of NO AVX section
							for(int64_t j=0;j<second->nelements;j++) {
								const DOUBLE dx = x2[j]-x1pos;
								const DOUBLE dy = y2[j]-y1pos;
								const DOUBLE dz = FABS(z2[j]-z1pos);
								
								const DOUBLE r2 = dx*dx + dy*dy;
								if(r2 >= sqr_rpmax || r2 < sqr_rpmin || dz >= pimax) {
									continue;
								}

#ifdef OUTPUT_RPAVG								
								const DOUBLE r = SQRT(r2);
#endif								
								int pibin = (int) (dz*inv_dpi);
								pibin = pibin > npibin ? npibin:pibin;
								for(int kbin=nrpbin-1;kbin>=1;kbin--) {
									if(r2 >= rupp_sqr[kbin-1]) {
										const int ibin = kbin*(npibin+1) + pibin;
										npairs[ibin]++;
#ifdef OUTPUT_RPAVG										
										rpavg[ibin]+=r;
#endif										
										break;
									}
								}
							}
							
#else //beginning of AVX section
							union int8 {
								AVX_INTS m_ibin;
								int ibin[NVEC];
							};
							union int8 union_rpbin;
							union int8 union_pibin;

#ifdef OUTPUT_RPAVG							
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
								const AVX_FLOATS x2pos = AVX_LOAD_FLOATS_UNALIGNED(&x2[j]);
								const AVX_FLOATS y2pos = AVX_LOAD_FLOATS_UNALIGNED(&y2[j]);
								const AVX_FLOATS z2pos = AVX_LOAD_FLOATS_UNALIGNED(&z2[j]);
								const AVX_FLOATS m_sqr_rpmax = AVX_SET_FLOAT(sqr_rpmax);
								const AVX_FLOATS m_sqr_rpmin = AVX_SET_FLOAT(sqr_rpmin);
								const AVX_FLOATS m_pimax = AVX_SET_FLOAT(pimax);
								const AVX_FLOATS m_zero  = AVX_SET_FLOAT((DOUBLE) 0.0);
								const AVX_FLOATS m_inv_dpi    = AVX_SET_FLOAT(inv_dpi);
								
								AVX_FLOATS m_zdiff       = AVX_SUBTRACT_FLOATS(z2pos,m_z1pos);
								const AVX_FLOATS m_xdiff = AVX_SUBTRACT_FLOATS(x2pos,m_x1pos);
								const AVX_FLOATS m_ydiff = AVX_SUBTRACT_FLOATS(y2pos,m_y1pos);
								
								m_zdiff = AVX_MAX_FLOATS(m_zdiff,AVX_SUBTRACT_FLOATS(m_zero,m_zdiff));//dz = fabs(dz) => dz = max(dz, -dz);
								
								AVX_FLOATS m_dist  = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_xdiff),AVX_SQUARE_FLOAT(m_ydiff));
								AVX_FLOATS m_mask_left;
								
								//Do all the distance cuts using masks here in new scope
								{
									const AVX_FLOATS m_mask_pimax = AVX_COMPARE_FLOATS(m_zdiff,m_pimax,_CMP_LT_OS);
									const int test = AVX_TEST_COMPARISON(m_mask_pimax);
									if(test == 0) {
										continue;
									}

									const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_dist,m_sqr_rpmin,_CMP_GE_OS);
									m_dist = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax,m_dist,m_mask_pimax);
									
									m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_sqr_rpmax,_CMP_LT_OS);//will get utilized in the next section
									const AVX_FLOATS m_mask = AVX_BITWISE_AND(m1,m_mask_left);
									int test1 = AVX_TEST_COMPARISON(m_mask);
									if(test1 == 0) {
										continue;
									}
									
									//So there's at least one point that is in range - let's find the bin
									m_zdiff = AVX_BLEND_FLOATS_WITH_MASK(m_pimax, m_zdiff, m_mask);
#ifdef OUTPUT_RPAVG
									union_mDperp.m_Dperp = AVX_SQRT_FLOAT(m_dist);
#endif									
									union_pibin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(AVX_MULTIPLY_FLOATS(m_zdiff,m_inv_dpi));
								}
								
								{
									AVX_FLOATS m_rpbin     = AVX_SET_FLOAT((DOUBLE) nrpbin);
									//AVX_FLOATS m_all_ones  = AVX_CAST_INT_TO_FLOAT(AVX_SET_INT(-1));
									for(int kbin=nrpbin-1;kbin>=1;kbin--) {
										const AVX_FLOATS m_mask_low = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_GE_OS);
										const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m_mask_low,m_mask_left);
										m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
										m_mask_left = AVX_COMPARE_FLOATS(m_dist, m_rupp_sqr[kbin-1],_CMP_LT_OS);
										//m_mask_left = AVX_XOR_FLOATS(m_mask_low, m_all_ones);//XOR with 0xFFFF... gives the bins that are smaller than m_rupp_sqr[kbin] (and is faster than cmp_p(s/d) in theory)
										const int test = AVX_TEST_COMPARISON(m_mask_left);
										if(test==0) {
											break;
										}
									}
									union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
								}
								
								//update the histograms
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma unroll(NVEC)
#endif		  
								for(int jj=0;jj<NVEC;jj++) {
									int rpbin = union_rpbin.ibin[jj];
									int pibin = union_pibin.ibin[jj];
									int ibin = rpbin*(npibin+1) + pibin;
									npairs[ibin]++;
#ifdef OUTPUT_RPAVG
									rpavg [ibin] += union_mDperp.Dperp[jj];
#endif									
								}
							}

							//remainder loop
							for(;j<second->nelements;j++) {
								const DOUBLE dx = x2[j]-x1pos;
								const DOUBLE dy = y2[j]-y1pos;
								const DOUBLE dz = FABS(z2[j]-z1pos);
								
								const DOUBLE r2 = dx*dx + dy*dy;
								if(r2 >= sqr_rpmax || r2 < sqr_rpmin || dz >= pimax) {
									continue;
								}
#ifdef OUTPUT_RPAVG								
								const DOUBLE r = SQRT(r2);
#endif								
								int pibin = (int) (dz*inv_dpi);
								pibin = pibin > npibin ? npibin:pibin;
								for(int kbin=nrpbin-1;kbin>=1;kbin--) {
									if(r2 >= rupp_sqr[kbin-1]) {
										int ibin = kbin*(npibin+1) + pibin;
										npairs[ibin]++;
#ifdef OUTPUT_RPAVG
										rpavg [ibin] += r;
#endif										
										break;
									}
								}
							} //end of j-remainder loop
							
#endif //end of AVX section
							
						}//end of i-loop over first
					}//iiz loop over zbin_refine_factor
				}//iiy loop over bin_refine_factor
			}//iix loop over bin_refine_factor
		}//index1 loop over totncells
#ifdef USE_OMP
		for(int i=0;i<totnbins;i++) {
			all_npairs[tid][i] = npairs[i];
#ifdef OUTPUT_RPAVG			
			all_rpavg[tid][i] = rpavg[i];
#endif
		}
	}//close the omp parallel region
#endif
	finish_myprogressbar(&interrupted);


	
#ifdef USE_OMP
	uint64_t npairs[totnbins];
#ifdef OUTPUT_RPAVG	
	DOUBLE rpavg[totnbins];
#endif	
	for(int i=0;i<totnbins;i++) {
		npairs[i] = 0;
#ifdef OUTPUT_RPAVG		
		rpavg[i] = 0.0;
#endif		
	}

	for(int i=0;i<numthreads;i++) {
		for(int j=0;j<totnbins;j++) {
			npairs[j] += all_npairs[i][j];
#ifdef OUTPUT_RPAVG
			rpavg[j] += all_rpavg[i][j];
#endif			
		}
	}
#endif

#ifdef OUTPUT_RPAVG	
	for(int i=0;i<totnbins;i++){
		if(npairs[i] > 0) {
			rpavg[i] /= ((DOUBLE) npairs[i] );
		}
	}
#endif


  //Pack in the results
	results_countpairs_rp_pi *results = my_malloc(sizeof(*results), 1);
	results->nbin   = nrpbin;
	results->npibin = npibin;
	results->pimax  = pimax;
	results->npairs = my_malloc(sizeof(uint64_t), totnbins);
	results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
	results->rpavg  = my_malloc(sizeof(DOUBLE)  , totnbins);

	for(int i=0;i<nrpbin;i++) {
		results->rupp[i] = rupp[i];
		for(int j=0;j<npibin;j++) {
			int index = i*(npibin+1) + j;
			results->npairs[index] = npairs[index];
#ifdef OUTPUT_RPAVG
			results->rpavg[index] = rpavg[index];
#else
			results->rpavg[index] = 0.0;
#endif
		}
	}

	free(rupp);
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
#ifdef USE_OMP
	matrix_free((void **) all_npairs, numthreads);
#ifdef OUTPUT_RPAVG	
	matrix_free((void **) all_rpavg, numthreads);
#endif	
#endif
	
	return results;
}
