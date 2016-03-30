/* File: countspheres_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>

#include "defs.h"
#include "countspheres_mocks.h" //function proto-type
#include "gridlink_mocks.h"//function proto-type for gridlink (which is a copy of the theory gridlink in order to avoid name-space collisions)
#include "utils.h" //all of the utilities
#include "set_cosmo_dist.h"//cosmological distance calculations
#include "cosmology_params.h"//init_cosmology

#ifndef SILENT
#include "progressbar.h" //for the progressbar
#endif

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif


int count_neighbors(const DOUBLE xcen,const DOUBLE ycen,const DOUBLE zcen,const DOUBLE smin,const DOUBLE inv_rcube,const DOUBLE rmax,
                    const int ngrid,const cellarray *lattice, const int nthreshold, const int bin_refine_factor)
{
    int numngb=0;
    const cellarray *cellstruct;
    const DOUBLE rmax_sqr = (DOUBLE) (rmax*rmax);
    int ix = (int)(ngrid*(xcen-smin)*inv_rcube);
    int iy = (int)(ngrid*(ycen-smin)*inv_rcube);
    int iz = (int)(ngrid*(zcen-smin)*inv_rcube);
    if(ix > ngrid-1) ix--;
    if(iy > ngrid-1) iy--;
    if(iz > ngrid-1) iz--;

    assert(ix >= 0 && ix < ngrid && "x-position is inside limits");
    assert(iy >= 0 && iy < ngrid && "y-position is inside limits");
    assert(iz >= 0 && iz < ngrid && "z-position is inside limits");

    const int min_ix = ix - bin_refine_factor < 0 ?             0:ix - bin_refine_factor;
    const int max_ix = ix + bin_refine_factor > ngrid-1 ? ngrid-1:ix + bin_refine_factor;
    for(int iix=min_ix;iix<=max_ix;iix++) {
        const DOUBLE newxpos = xcen;
        const int min_iy = iy - bin_refine_factor < 0 ?             0:iy - bin_refine_factor;
        const int max_iy = iy + bin_refine_factor > ngrid-1 ? ngrid-1:iy + bin_refine_factor;

        for(int iiy=min_iy;iiy<=max_iy;iiy++) {
            const DOUBLE newypos = ycen;
            const int min_iz = iz - bin_refine_factor < 0 ?             0:iz - bin_refine_factor;
            const int max_iz = iz + bin_refine_factor > ngrid-1 ? ngrid-1:iz + bin_refine_factor;

            for(int iiz=min_iz;iiz<=max_iz;iiz++) {
                DOUBLE newzpos = zcen;
                const int64_t index=iix*ngrid*ngrid + iiy*ngrid + iiz;
                cellstruct = &(lattice[index]);
                DOUBLE *x2 = cellstruct->pos;
                DOUBLE *y2 = cellstruct->pos + NVEC;
                DOUBLE *z2 = cellstruct->pos + 2*NVEC;

                for(int i=0;i<cellstruct->nelements;i+=NVEC) {
                    int block_size = cellstruct->nelements - i ;
                    if(block_size > NVEC) block_size = NVEC;
                    for(int ii=0;ii<block_size;ii++) {
                        const DOUBLE dx=x2[ii]-newxpos;
                        const DOUBLE dy=y2[ii]-newypos;
                        const DOUBLE dz=z2[ii]-newzpos;
                        const DOUBLE r2 = dx*dx + dy*dy + dz*dz;
                        if (r2 < rmax_sqr) numngb++;
                    }
                    if(numngb > nthreshold) return numngb;

                    x2 += 3*NVEC;
                    y2 += 3*NVEC;
                    z2 += 3*NVEC;

                }
            }
        }
    }
    return numngb;
}


void free_results_countspheres_mocks(results_countspheres_mocks **results)
{
    if(results == NULL)
        return;
    if(*results == NULL)
        return;

    results_countspheres_mocks *tmp = *results;
    matrix_free((void **) tmp->pN, tmp->nbin);
    free(tmp);
    tmp = NULL;
}





results_countspheres_mocks * countspheres_mocks(const int64_t Ngal, DOUBLE *xgal, DOUBLE *ygal, DOUBLE *zgal,
                                                const int64_t Nran, DOUBLE *xran, DOUBLE *yran, DOUBLE *zran,
                                                const int threshold_neighbors,
                                                const DOUBLE rmax, const int nbin, const int nc,
                                                const int num_pN,
                                                const char *centers_file,
                                                const int cosmology)

{
    /*---Gridlink-variables----------------*/
    int ngrid;


    //*---Measurement-----------------------*/
    int itry,isucceed;

    init_cosmology(cosmology);

    //Input validation
    assert(rmax > 0.0  && "rmax has to be positive");
    assert(nbin >=1    && "Number of bins has to be at least 1");
    assert(nc >=1      && "Number of spheres has to be at least 1");
    assert(num_pN >= 1 && "Number of pN's requested must be at least 1");

    int need_randoms=0;
    int64_t num_centers_in_file=0;
    FILE *fpcen = fopen(centers_file,"r");
    if(fpcen != NULL) {
        double rr = 0.0;
        int num_read = fscanf(fpcen,"%*f %*f %*f %lf",&rr);
        assert(num_read == 1 && "Could not read max. sphere radius from the centers file");
        num_centers_in_file = getnumlines(centers_file,'#');
        if( rr >= rmax && num_centers_in_file >= nc) {
            need_randoms = 0;
            rewind(fpcen);
        } else {
            fclose(fpcen);
            num_centers_in_file = 0;
            need_randoms = 1;
        }
    } else {
        num_centers_in_file = 0;
        need_randoms = 1;
    }
    if(need_randoms==1) {
        fpcen = my_fopen(centers_file,"w");
    }
#ifndef SILENT
    fprintf(stderr,"%s> found %"PRId64" centers (need %d centers) - need randoms = %d\n",__FUNCTION__,num_centers_in_file,nc,need_randoms);
#endif


    //set up the interpolation for comoving distances
    double *redshifts,*comoving_distance;
    redshifts=my_calloc(sizeof(*redshifts),COSMO_DIST_SIZE);
    comoving_distance=my_calloc(sizeof(*comoving_distance),COSMO_DIST_SIZE);
    int Nzdc = set_cosmo_dist(MAX_REDSHIFT_FOR_COSMO_DIST, COSMO_DIST_SIZE, redshifts, comoving_distance, cosmology);
    const DOUBLE inv_speed_of_light=1.0/SPEED_OF_LIGHT;

    gsl_interp *interpolation;
    gsl_interp_accel *accelerator;
    accelerator =  gsl_interp_accel_alloc();
    interpolation = gsl_interp_alloc (gsl_interp_linear,Nzdc);
    gsl_interp_init(interpolation, redshifts, comoving_distance, Nzdc);

    DOUBLE rcube=0.0 ;
    for(int i=0;i<Ngal;i++) {
        const DOUBLE new_phi   = xgal[i];
        const DOUBLE new_theta = ygal[i];
        const DOUBLE new_cz    = zgal[i];
        const DOUBLE dc = gsl_interp_eval(interpolation, redshifts, comoving_distance, new_cz*inv_speed_of_light, accelerator);
        if(dc>rcube) rcube = dc;

        xgal[i] = dc*COSD(new_theta)*COSD(new_phi) ;
        ygal[i] = dc*COSD(new_theta)*SIND(new_phi) ;
        zgal[i] = dc*SIND(new_theta) ;
    }

    if (need_randoms == 1) {
        for(int i=0;i<Nran;i++) {
            const DOUBLE new_phi   = xran[i];
            const DOUBLE new_theta = yran[i];
            const DOUBLE new_cz    = zran[i];
            const DOUBLE dc = gsl_interp_eval(interpolation, redshifts, comoving_distance, new_cz*inv_speed_of_light, accelerator);

            xran[i] = dc*COSD(new_theta)*COSD(new_phi) ;
            yran[i] = dc*COSD(new_theta)*SIND(new_phi) ;
            zran[i] = dc*SIND(new_theta);
        }
    }
    free(redshifts);free(comoving_distance);
    gsl_interp_free(interpolation);
    gsl_interp_accel_free(accelerator);


    /*---Shift-coordinates--------------------------------*/
#ifndef SILENT
    fprintf(stderr,"%s> maximum distance = %f. ",__FUNCTION__,rcube) ;
#endif
    rcube = rcube + 1. ; //add buffer

    for(int i=0;i<Ngal;i++) {
        xgal[i] += rcube ;
        ygal[i] += rcube ;
        zgal[i] += rcube ;
    }

    if(need_randoms == 1) {
        for(int i=0;i<Nran;i++) {
            xran[i] += rcube ;
            yran[i] += rcube ;
            zran[i] += rcube ;
        }
    }
    rcube = 2.0*rcube ;
    const DOUBLE inv_rcube = 1.0/rcube;
    const DOUBLE rmax_sqr = rmax*rmax;
#ifndef SILENT
    fprintf(stderr," Bounding cube size = %f\n",rcube) ;
#endif

    /*---Construct-grid-to-speed-up-neighbor-searching----*/
    //First create the 3-d linklist
    int bin_refine_factor=1;
    cellarray *lattice=NULL;//pointer to the full 3-d volume for galaxies
    cellarray *randoms_lattice=NULL;//pointer to the full 3-d volume for randoms
    const DOUBLE xmin=0.0,xmax=rcube;
    const DOUBLE ymin=0.0,ymax=rcube;
    const DOUBLE zmin=0.0,zmax=rcube;
    const DOUBLE smin=0.0;

    {
        //new scope -> no need for the nmesh_x/y/z variables since the data is in a cube
        int nmesh_x,nmesh_y,nmesh_z;
        lattice = gridlink(Ngal, xgal, ygal, zgal, xmin, xmax, ymin, ymax, zmin, zmax, rmax, rmax, rmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
        assert(nmesh_x == nmesh_y && nmesh_x == nmesh_z && "The number of grid cells should be identical");
        ngrid = nmesh_x;//could have been nmesh_y/nmesh_z -> all three are equal
    }

    if(need_randoms == 1) {
        int nran_x,nran_y,nran_z;
        randoms_lattice = gridlink(Nran, xran, yran, zran, xmin, xmax, ymin, ymax, zmin, zmax, rmax, rmax, rmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nran_x, &nran_y, &nran_z);
        assert(nran_x == nran_y && nran_x == nran_z && "The number of (randoms) grid cells should be identical");
        assert(nran_x == ngrid && "The number of grid cells for randoms should be identical to that in the data");
    }


    /*---Prepare-radial-arrays----------------------------*/
    int *counts = my_calloc(sizeof(*counts),nbin);
    DOUBLE **pN = (DOUBLE **) matrix_calloc(sizeof(DOUBLE), nbin, num_pN);

    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    const DOUBLE inv_rstep = ((DOUBLE) 1.0)/rstep;

#if defined(USE_AVX) && defined(__AVX__)
    AVX_FLOATS m_rupp_sqr[nbin];
    AVX_FLOATS m_rmax_sqr = AVX_SET_FLOAT(rmax_sqr);
    for(int k=0;k<nbin;k++) {
        m_rupp_sqr[k] = AVX_SET_FLOAT((k+1)*rstep*(k+1)*rstep);
    }
#endif


    itry=0 ;
    isucceed=0 ;
    int ncenters_written=0;

#ifndef SILENT    
    int interrupted;
    init_my_progressbar(nc, &interrupted);
#endif    
    while(isucceed < nc && itry < Nran) {

#ifndef SILENT      
      my_progressbar(isucceed,&interrupted);
#endif      

      DOUBLE xcen,ycen,zcen;
      int Nnbrs_ran=0;
      if((need_randoms == 1 && isucceed > num_centers_in_file) || num_centers_in_file == 0) {
            xcen = xran[itry] ;
            ycen = yran[itry] ;
            zcen = zran[itry] ;
            Nnbrs_ran = count_neighbors(xcen,ycen,zcen,smin,inv_rcube,rmax,ngrid,randoms_lattice, threshold_neighbors, bin_refine_factor);
        } else {
            double rr=0.0;
            const int MAXBUFSIZE=10000;
            char buffer[MAXBUFSIZE];
            assert( fgets(buffer,MAXBUFSIZE,fpcen) != NULL && "ERROR: Could not read-in centers co-ordinates");
            int nitems = sscanf(buffer,"%"DOUBLE_FORMAT" %"DOUBLE_FORMAT" %"DOUBLE_FORMAT" %lf",&xcen,&ycen,&zcen,&rr);
            if(nitems != 4) {
                fprintf(stderr,"ERROR: nitems = %d xcen = %lf ycen = %lf zcen %lf rr = %lf\n",
                        nitems,xcen,ycen,zcen,rr);
                fprintf(stderr,"buffer = `%s' \n",buffer);
            }
            assert(nitems == 4 && "Read the centers from the centers file");
            assert(rr >= rmax && "Rmax from the center file is >= rmax");
            Nnbrs_ran = threshold_neighbors + 1;
        }

        if(Nnbrs_ran > threshold_neighbors) {  //ignore if sphere overlaps edge
            for(int k=0;k<nbin;k++) {  //initialize counts
                counts[k] = 0 ;
            }

            int ix = (int)(ngrid*(xcen-smin)*inv_rcube);
            int iy = (int)(ngrid*(ycen-smin)*inv_rcube);
            int iz = (int)(ngrid*(zcen-smin)*inv_rcube);
            if(ix > ngrid-1) ix--;
            if(iy > ngrid-1) iy--;
            if(iz > ngrid-1) iz--;

            assert(ix >= 0 && ix < ngrid && "x-position is inside limits");
            assert(iy >= 0 && iy < ngrid && "y-position is inside limits");
            assert(iz >= 0 && iz < ngrid && "z-position is inside limits");

            const int min_ix = ix - bin_refine_factor < 0 ?             0:ix - bin_refine_factor;
            const int max_ix = ix + bin_refine_factor > ngrid-1 ? ngrid-1:ix + bin_refine_factor;
            for(int iix=min_ix;iix<=max_ix;iix++) {
                const DOUBLE newxpos = xcen;
#if defined(USE_AVX) && defined(__AVX__)
                const AVX_FLOATS m_newxpos = AVX_SET_FLOAT(newxpos);
#endif

                const int min_iy = iy - bin_refine_factor < 0 ?             0:iy - bin_refine_factor;
                const int max_iy = iy + bin_refine_factor > ngrid-1 ? ngrid-1:iy + bin_refine_factor;

                for(int iiy=min_iy;iiy<=max_iy;iiy++) {
                    const DOUBLE newypos = ycen;
#if defined(USE_AVX) && defined(__AVX__)
                    const AVX_FLOATS m_newypos = AVX_SET_FLOAT(newypos);
#endif

                    const int min_iz = iz - bin_refine_factor < 0 ?             0:iz - bin_refine_factor;
                    const int max_iz = iz + bin_refine_factor > ngrid-1 ? ngrid-1:iz + bin_refine_factor;

                    for(int iiz=min_iz;iiz<=max_iz;iiz++) {
                        const DOUBLE newzpos = zcen;
#if defined(USE_AVX) && defined(__AVX__)
                        const AVX_FLOATS m_newzpos = AVX_SET_FLOAT(newzpos);
#endif
                        const int index=iix*ngrid*ngrid + iiy*ngrid + iiz;
                        const cellarray *cellstruct = &(lattice[index]);
                        DOUBLE *x2 = cellstruct->pos;
                        DOUBLE *y2 = cellstruct->pos + NVEC;
                        DOUBLE *z2 = cellstruct->pos + 2*NVEC;
                        int ipart;
                        for(ipart=0;ipart<=(cellstruct->nelements-NVEC);ipart+=NVEC) {
#if !(defined(USE_AVX) && defined(__AVX__))
                            int ibin[NVEC];
#if  __INTEL_COMPILER
#pragma simd vectorlengthfor(DOUBLE)
#endif
                            for(int k=0;k<NVEC;k++) {
                                const DOUBLE dx=x2[k]-newxpos;
                                const DOUBLE dy=y2[k]-newypos;
                                const DOUBLE dz=z2[k]-newzpos;
                                const DOUBLE r = SQRT(dx*dx + dy*dy + dz*dz);
                                ibin[k] = (int) (r*inv_rstep);
                            }
                            x2 += 3*NVEC;
                            y2 += 3*NVEC;
                            z2 += 3*NVEC;

#ifdef  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
                            for(int k=0;k<NVEC;k++) {
                                if(ibin[k] < nbin) counts[ibin[k]]++;
                            }


                            //Here is the AVX part
#else
                            const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(x2);
                            const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(y2);
                            const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(z2);

                            x2 += 3*NVEC;
                            y2 += 3*NVEC;
                            z2 += 3*NVEC;

                            const AVX_FLOATS m_xdiff = AVX_SUBTRACT_FLOATS(m_x2,m_newxpos);
                            const AVX_FLOATS m_ydiff = AVX_SUBTRACT_FLOATS(m_y2,m_newypos);
                            const AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_z2,m_newzpos);
                            AVX_FLOATS m_dist  = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_xdiff),AVX_SQUARE_FLOAT(m_ydiff));
                            m_dist = AVX_ADD_FLOATS(m_dist,AVX_SQUARE_FLOAT(m_zdiff));
                            AVX_FLOATS m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_rmax_sqr,_CMP_LT_OS);
                            const int test = AVX_TEST_COMPARISON(m_mask_left);
                            if(test == 0)
                                continue;

                            for(int kbin=nbin-1;kbin>=1;kbin--) {
                                const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_GE_OS);
                                const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
                                const int test2  = AVX_TEST_COMPARISON(m_bin_mask);
                                counts[kbin] += AVX_BIT_COUNT_INT(test2);
                                m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_LT_OS);
                                const int test3 = AVX_TEST_COMPARISON(m_mask_left);
                                if(test3 == 0) {
                                    break;
                                } else if(kbin==1){
                                    counts[0] += AVX_BIT_COUNT_INT(test3);
                                }
                            }

#endif //endof AVX section
                        }

                        //Take care of the rest
                        for(int ipos=0;ipart < cellstruct->nelements;ipart++,ipos++) {
                            const DOUBLE dx=x2[ipos]-newxpos;
                            const DOUBLE dy=y2[ipos]-newypos;
                            const DOUBLE dz=z2[ipos]-newzpos;
                            const DOUBLE r2 = (dx*dx + dy*dy + dz*dz);
                            if(r2 >= rmax_sqr) continue;
                            const int ibin = (int) (SQRT(r2)*inv_rstep);
                            counts[ibin]++;
                        }
                    }
                }
            }
            //Output the center into the file -> either
            if((need_randoms == 1 && isucceed > num_centers_in_file) || num_centers_in_file == 0) {
                fprintf(fpcen,"%lf \t %lf \t %lf \t %lf\n",xcen,ycen,zcen,rmax);
                ncenters_written++;
            }


            /* compute cumulative counts, i.e. n1 changes from the number of galaxies
               in shell ibin to  the number of galaxies in shell ibin or any smaller shell */
            for(int ibin=1;ibin<nbin;ibin++){
                counts[ibin]+=counts[ibin-1];
            }

            for(int ibin=0;ibin<nbin;ibin++) { //compute statistics
                for(int i=0;i<num_pN;i++) {
                    if(counts[ibin] == i) {
                        pN[ibin][i] += (DOUBLE) 1.0;
                    }
                }
            }
            isucceed++ ;
        }
        itry++ ;
    }
    fclose(fpcen);

#ifndef SILENT
    finish_myprogressbar(&interrupted);
    fprintf(stderr,"%s> Placed %d centers out of %d trials.\n",__FUNCTION__,isucceed,itry);
    fprintf(stderr,"%s> num_centers_in_file = %"PRId64" ncenters_written = %d\n",__FUNCTION__,num_centers_in_file,ncenters_written);
#endif
    assert(isucceed > 0 && "Placed > 0 spheres in the volume");
    if(isucceed < nc) {
        fprintf(stderr,"WARNING: Could only place `%d' out of requested `%d' spheres. Increase the random-sample size might improve the situation\n",isucceed,nc);
    }

    /* free(Navg);free(Nvar);free(P0);free(P1);free(P2); */
    free(counts);
    int64_t totncells = ngrid*ngrid*ngrid;
    for(int64_t icell=0;icell < totncells;icell++) {
        free(lattice[icell].pos);
        /* free(lattice[icell].y); */
        /* free(lattice[icell].z); */
        if(need_randoms == 1) {
            free(randoms_lattice[icell].pos);
            /* free(randoms_lattice[icell].y); */
            /* free(randoms_lattice[icell].z); */
        }
    }

    free(lattice);
    if(need_randoms == 1) {
        free(randoms_lattice);
    }
    //prepare the results
    results_countspheres_mocks *results = my_malloc(sizeof(*results),1);
    results->rmax = rmax;
    results->nbin = nbin;
    results->nc   = nc;
    results->num_pN = num_pN;
    results->pN = (DOUBLE **) matrix_malloc(sizeof(DOUBLE), nbin, num_pN);
    const DOUBLE inv_nc = ((DOUBLE) 1.0)/(DOUBLE) isucceed;//actual number of spheres placed
    for(int i=0;i<num_pN;i++) {
        for(int ibin=0;ibin<nbin;ibin++) {
            (results->pN)[ibin][i] = pN[ibin][i] * inv_nc;
        }
    }
    matrix_free((void **) pN, nbin);

    return results;
}
