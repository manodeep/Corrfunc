#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

#include "gridlink.h"//function proto-type for gridlink
#include "countspheres.h" //function proto-type
#include "cellarray.h" //definition of struct cellarray
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar

#ifdef USE_AVX
#include "avx_calls.h"
#endif


void free_results_countspheres(results_countspheres **results)
{
	if(results == NULL)
		return;
	if(*results == NULL)
		return;

	results_countspheres *tmp = *results;
	matrix_free((void **) tmp->pN, tmp->nbin);
	free(tmp);
	tmp = NULL;
}



results_countspheres * countspheres(const int64_t np, const DOUBLE * restrict X, const DOUBLE * restrict Y, const DOUBLE * restrict Z,
																		const double rmax, const unsigned int nbin, const unsigned int nc,
																		const unsigned int num_pN,
																		unsigned long seed)

{
	struct timeval t0,t1;
  const gsl_rng_type * T = gsl_rng_mt19937;
  gsl_rng * rng;
	int bin_refine_factor=1;
	
	rng = gsl_rng_alloc (T);
  gsl_rng_set(rng, seed);

	uint64_t *counts = my_calloc(sizeof(*counts),nbin);
	DOUBLE **pN = (DOUBLE **) matrix_calloc(sizeof(DOUBLE), nbin, num_pN);
	
  const DOUBLE rstep = rmax/(DOUBLE)nbin ;
  const DOUBLE inv_rstep = ((DOUBLE) 1.0)/rstep;
  /* const DOUBLE inv_rcube = ((DOUBLE) 1.0)/rcube; */
  const DOUBLE rmax_sqr = rmax*rmax;
#ifdef USE_AVX  
  AVX_FLOATS m_rmax_sqr = AVX_SET_FLOAT(rmax_sqr);
  AVX_FLOATS m_rupp_sqr[nbin];
  for(unsigned int k=0;k<nbin;k++) {
    m_rupp_sqr[k] = AVX_SET_FLOAT((k+1)*rstep*rstep*(k+1));
  }
#endif

	//Find the min/max of the data
	DOUBLE xmin,xmax,ymin,ymax,zmin,zmax;
	xmin=1e10;ymin=1e10;zmin=1e10;
	xmax=0.0;ymax=0.0;zmax=0.0;
	get_max_min(np, X, Y, Z, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax);
	
  //First create the 3-d linklist
  int nmesh_x=0,nmesh_y=0,nmesh_z=0;
  gettimeofday(&t0,NULL);
	cellarray *lattice = gridlink(np, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rmax, rmax, rmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
  gettimeofday(&t1,NULL);
	int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z; 

#ifndef SILENT
	fprintf(stderr,"Running with [xmin,xmax] = %lf,%lf\n",xmin,xmax);
	fprintf(stderr,"Running with [ymin,ymax] = %lf,%lf\n",ymin,ymax);
	fprintf(stderr,"Running with [zmin,zmax] = %lf,%lf\n",zmin,zmax);
#endif

	const DOUBLE xdiff = (xmax-xmin);
	const DOUBLE ydiff = (ymax-ymin);
	const DOUBLE zdiff = (zmax-zmin);
	const DOUBLE inv_xdiff = ((DOUBLE) 1.0)/xdiff;
	const DOUBLE inv_ydiff = ((DOUBLE) 1.0)/ydiff;
	const DOUBLE inv_zdiff = ((DOUBLE) 1.0)/zdiff;

  int interrupted=0;
  init_my_progressbar(nc,&interrupted);

  /* loop through centers, placing each randomly */
	unsigned int ic=0;
	while(ic < nc) {
    my_progressbar(ic,&interrupted);
    const DOUBLE xc = xdiff*gsl_rng_uniform (rng) + xmin;
    const DOUBLE yc = ydiff*gsl_rng_uniform (rng) + ymin;
    const DOUBLE zc = zdiff*gsl_rng_uniform (rng) + zmin;

#ifndef PERIODIC
		//Check that the biggest sphere will not intersect
		//with the box edges iff non-periodic conditions are set
		if((xc - xmin) < rmax || (xmax - xc) < rmax ||
			 (yc - ymin) < rmax || (ymax - yc) < rmax ||
			 (zc - zmin) < rmax || (zmax - zc) < rmax) {
		
			continue;
		}
#endif		
		ic++;
		
    int ix = (int)(nmesh_x*(xc-xmin)*inv_xdiff);
    int iy = (int)(nmesh_y*(yc-ymin)*inv_ydiff);
    int iz = (int)(nmesh_z*(zc-zmin)*inv_zdiff);
    if(ix > nmesh_x-1) ix--;
    if(iy > nmesh_y-1) iy--;
    if(iz > nmesh_z-1) iz--;

    assert(ix >= 0 && ix < nmesh_x && "x-position is inside limits");
    assert(iy >= 0 && iy < nmesh_y && "y-position is inside limits");
    assert(iz >= 0 && iz < nmesh_z && "z-position is inside limits");
    
    for(unsigned int ibin=0;ibin<nbin;ibin++) {
      counts[ibin]=0;
    }

    gettimeofday(&t0,NULL);
    for(int iix=-bin_refine_factor;iix<=bin_refine_factor;iix++) {
			int iiix;
#ifdef PERIODIC
			DOUBLE off_xwrap=0.0;
			if(ix + iix >= nmesh_x) {
				off_xwrap = -xdiff;
			} else if (ix + iix < 0) {
				off_xwrap = xdiff;
			}
			iiix=(ix+iix+nmesh_x)%nmesh_x;
			const DOUBLE newxpos = xc + off_xwrap;
#else
			iiix = iix+ix;
			if(iiix < 0 || iiix >= nmesh_x) {
				continue;
			}
			const DOUBLE newxpos = xc;
#endif  

      for(int iiy=-bin_refine_factor;iiy<=bin_refine_factor;iiy++) {
				int iiiy;
#ifdef PERIODIC
				DOUBLE off_ywrap = 0.0;
				if(iy + iiy >= nmesh_y) {
					off_ywrap = -ydiff;
				} else if (iy + iiy < 0) {
					off_ywrap = ydiff;
				}
				iiiy=(iy+iiy+nmesh_y)%nmesh_y;
				const DOUBLE newypos = yc + off_ywrap;
#else
				iiiy = iiy+iy;
				if(iiiy < 0 || iiiy >= nmesh_y) {
					continue;
				}
				const DOUBLE newypos = yc;
#endif

				for(int iiz=-bin_refine_factor;iiz<=bin_refine_factor;iiz++) {
					int iiiz;
#ifdef PERIODIC
					DOUBLE off_zwrap = 0.0;
					if(iz + iiz >= nmesh_z) {
						off_zwrap = -zdiff;
					} else if (iz + iiz < 0) {
						off_zwrap = zdiff;
					}
					iiiz=(iz+iiz+nmesh_z)%nmesh_z;
					const DOUBLE newzpos = zc + off_zwrap;
#else
					iiiz = iiz+iz;
					if(iiiz < 0 || iiiz >= nmesh_z) {
						continue;
					}
					const DOUBLE newzpos = zc;
#endif

					const int64_t index=iiix*nmesh_y*nmesh_z + iiiy*nmesh_z + iiiz;
					const cellarray *first = &(lattice[index]);
					const DOUBLE *x2 = first->x;
					const DOUBLE *y2 = first->y;
					const DOUBLE *z2 = first->z;
#ifndef USE_AVX
					for(int64_t j=0;j<first->nelements;j+=NVEC) {
						int block_size=first->nelements - j;
						if(block_size > NVEC) block_size=NVEC;
						for(int jj=0;jj<block_size;jj++) {
							int ibin;
							DOUBLE dx=x2[j+jj]-newxpos;
							DOUBLE dy=y2[j+jj]-newypos;
							DOUBLE dz=z2[j+jj]-newzpos;
							DOUBLE r2 = dx*dx + dy*dy + dz*dz;
							if(r2 >= rmax_sqr) continue;
							ibin = (int) (SQRT(r2)*inv_rstep);
							counts[ibin]++;
						}
					}
#else //beginning of AVX section

					const AVX_FLOATS m_xc    = AVX_SET_FLOAT(newxpos);
					const AVX_FLOATS m_yc    = AVX_SET_FLOAT(newypos);
					const AVX_FLOATS m_zc    = AVX_SET_FLOAT(newzpos);

					int64_t j;
					for(j=0;j<=(first->nelements-NVEC);j+=NVEC) {

						const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(&x2[j]);
						const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(&y2[j]);
						const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(&z2[j]);
						
						const AVX_FLOATS m_dx = AVX_SUBTRACT_FLOATS(m_xc,m_x1);
						const AVX_FLOATS m_dy = AVX_SUBTRACT_FLOATS(m_yc,m_y1);
						const AVX_FLOATS m_dz = AVX_SUBTRACT_FLOATS(m_zc,m_z1);
						
						const AVX_FLOATS m_r2 = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dx),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dy),AVX_SQUARE_FLOAT(m_dz)));
						const AVX_FLOATS m_mask = AVX_COMPARE_FLOATS(m_r2,m_rmax_sqr,_CMP_LT_OS);
						if(AVX_TEST_COMPARISON(m_mask) == 0) {
							continue;
						}
						
						for(int k=nbin-1;k>=1;k--){
							AVX_FLOATS m_mask1 = AVX_COMPARE_FLOATS(m_r2,m_rupp_sqr[k],_CMP_LT_OS);
							AVX_FLOATS m_mask2 = AVX_COMPARE_FLOATS(m_r2,m_rupp_sqr[k-1],_CMP_GE_OS);
							AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m_mask1,m_mask2);
							int test2 = AVX_TEST_COMPARISON(m_bin_mask);
							counts[k] += AVX_BIT_COUNT_INT(test2);
							AVX_FLOATS m_mask_left = AVX_COMPARE_FLOATS(m_r2,m_rupp_sqr[k-1],_CMP_LT_OS);
							int test3 = AVX_TEST_COMPARISON(m_mask_left);
							if(test3 == 0) {
								break;
							} else if(k==1){
								counts[0] += AVX_BIT_COUNT_INT(test3);
							}
						}
					}


					//Take care of the rest
					for(;j<first->nelements;j++) {
						const DOUBLE dx=x2[j]-newxpos;
						const DOUBLE dy=y2[j]-newypos;
						const DOUBLE dz=z2[j]-newzpos;
						const DOUBLE r2 = (dx*dx + dy*dy + dz*dz);
						if(r2 >= rmax_sqr) continue;
						const unsigned int ibin = (unsigned int) (SQRT(r2)*inv_rstep);
						if(ibin < nbin) counts[ibin]++;
					}
#endif // end of AVX section
				}
      }
    }
	
		/* compute cumulative counts, i.e. counts changes from the number of galaxies
			 in shell ibin to  the number of galaxies in shell ibin or any smaller shell */
		for(unsigned int ibin=1;ibin<nbin;ibin++){
			counts[ibin]+=counts[ibin-1];
    }
    
    /* compute pN's */
    for(unsigned int ibin=0;ibin<nbin;ibin++) {
			for(unsigned int i=0;i<num_pN;i++) {
				if(counts[ibin] == i) {
					pN[ibin][i] += (DOUBLE) 1.0;
				}
			}
    }
    gettimeofday(&t1,NULL);
  }//loop over number of spheres

	free(counts);
  gsl_rng_free (rng);
  for(int i=0;i<totncells;i++) {
    free(lattice[i].x);
		free(lattice[i].y);
		free(lattice[i].z);
  }
  free(lattice);
	finish_myprogressbar(&interrupted);

	//prepare the results
	results_countspheres *results = my_malloc(sizeof(*results),1);
	results->rmax = rmax;
	results->nbin = nbin;
	results->nc   = nc;
	results->num_pN = num_pN;
	results->pN = (DOUBLE **) matrix_malloc(sizeof(DOUBLE), nbin, num_pN);
	const DOUBLE inv_nc = ((DOUBLE) 1.0)/(DOUBLE) nc;
	for(unsigned int i=0;i<num_pN;i++) {
		for(unsigned int ibin=0;ibin<nbin;ibin++) {
			(results->pN)[ibin][i] = pN[ibin][i] * inv_nc;
		}
	}
	matrix_free((void **) pN, nbin);

	return results;
}
