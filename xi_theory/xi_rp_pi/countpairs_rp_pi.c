/* File: countpairs_rp_pi.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "countpairs_rp_pi.h" //function proto-type
#include "gridlink.h"//function proto-type for gridlink
#include "cellarray.h" //definition of struct cellarray
#include "sglib.h"//sorting
#include "utils.h" //all of the utilities

#ifndef SILENT
#include "progressbar.h" //for the progressbar
#endif

#include "countpairs_rp_pi_driver.h"

#if defined(USE_OMP) && defined(_OPENMP)
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


results_countpairs_rp_pi * countpairs_rp_pi(const int64_t ND1, DOUBLE *X1, DOUBLE *Y1, DOUBLE *Z1,
                                            const int64_t ND2, DOUBLE *X2, DOUBLE *Y2, DOUBLE *Z2,
#if defined(USE_OMP) && defined(_OPENMP)
                                            const int numthreads,
#endif
                                            const int autocorr,
                                            const char *binfile,
                                            const DOUBLE pimax)
{
    int bin_refine_factor=2;
    int zbin_refine_factor=1;
    const int npibin = (int) pimax;
#ifdef PERIODIC
    const int periodic = 1;
#else
    const int periodic = 0;
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
    const DOUBLE xdiff = (xmax-xmin);
    const DOUBLE ydiff = (ymax-ymin);
    const DOUBLE zdiff = (zmax-zmin);
    if(rpmax < 0.05*xdiff) bin_refine_factor = 1;

    /*---Create 3-D lattice--------------------------------------*/
    int nmesh_x=0,nmesh_y=0,nmesh_z=0;
    cellarray_index *lattice1 = gridlink_index(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax,pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
        fprintf(stderr,"countpairs> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n");
        bin_refine_factor *=2;
        free(lattice1);
        lattice1 = gridlink_index(ND1, X1, Y1, Z1, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    }

    
    cellarray_index *lattice2 = NULL;
    if(autocorr==0) {
        int ngrid2_x=0,ngrid2_y=0,ngrid2_z=0;
        lattice2 = gridlink_index(ND2, X2, Y2, Z2, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &ngrid2_x, &ngrid2_y, &ngrid2_z);
        assert(nmesh_x == ngrid2_x && "Both lattices have the same number of X bins");
        assert(nmesh_y == ngrid2_y && "Both lattices have the same number of Y bins");
        assert(nmesh_z == ngrid2_z && "Both lattices have the same number of Z bins");
    } else {
        lattice2 = lattice1;
    }
    const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;

    //Generate the unique set of neighbouring cells to count over. 
    assign_ngb_cells(lattice1, lattice2, totncells, bin_refine_factor, bin_refine_factor, zbin_refine_factor, nmesh_x, nmesh_y, nmesh_z, xdiff, ydiff, zdiff, autocorr, periodic);


#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,x,i,j); \
            SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,y,i,j);               \
            SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,z,i,j) }

#if defined(USE_OMP) && defined(_OPENMP)
    omp_set_num_threads(numthreads);
#pragma omp parallel for schedule(dynamic)
#endif
    for(int64_t icell=0;icell<totncells;icell++) {
        const cellarray_index *first=&(lattice1[icell]);
        if(first->nelements > 0) {
          const int64_t start = first->start;
          DOUBLE *x = X1 + start;
          DOUBLE *y = Y1 + start;
          DOUBLE *z = Z1 + start;
          
          SGLIB_ARRAY_QUICK_SORT(DOUBLE, z, first->nelements, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
        }

    }
    if(autocorr == 0) {
#if defined(USE_OMP) && defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
        for(int64_t icell=0;icell<totncells;icell++) {
            const cellarray_index *first=&(lattice2[icell]);
            if(first->nelements > 0) {
                const int64_t start = first->start;
                DOUBLE *x = X2 + start;
                DOUBLE *y = Y2 + start;
                DOUBLE *z = Z2 + start;
                
                SGLIB_ARRAY_QUICK_SORT(DOUBLE, z, first->nelements, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
            }
#undef MULTIPLE_ARRAY_EXCHANGER
        }
    }

    DOUBLE rupp_sqr[nrpbin];
    const int64_t totnbins = (npibin+1)*(nrpbin+1);
    for(int i=0; i < nrpbin;i++) {
        rupp_sqr[i] = rupp[i]*rupp[i];
    }

    const DOUBLE sqr_rpmax=rupp_sqr[nrpbin-1];
    const DOUBLE sqr_rpmin=rupp_sqr[0];

#if !(defined(USE_OMP) && defined(_OPENMP))
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
#else//Use OMP
    omp_set_num_threads(numthreads);
    uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, totnbins);
#ifdef OUTPUT_RPAVG
    DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE),numthreads,totnbins);
#endif
#endif//OMP


#ifndef SILENT    
    int interrupted=0;
    int64_t numdone=0;
    init_my_progressbar(totncells,&interrupted);
#endif

#if defined(USE_OMP) && defined(_OPENMP)
#ifndef SILENT    
#pragma omp parallel shared(numdone)
#else
#pragma omp parallel    
#endif    
    {
        const int tid = omp_get_thread_num();
        uint64_t npairs[totnbins];
        for(int i=0;i<totnbins;i++) npairs[i] = 0;
#ifdef OUTPUT_RPAVG
        DOUBLE rpavg[totnbins];
        for(int i=0;i<totnbins;i++) rpavg[i] = ZERO;
#endif

#pragma omp for  schedule(dynamic) nowait
#endif
        /*---Loop-over-lattice1--------------------*/
        for(int64_t index1=0;index1<totncells;index1++) {

#ifndef SILENT
#if defined(USE_OMP) && defined(_OPENMP)
            if (omp_get_thread_num() == 0)
#endif
                my_progressbar(numdone,&interrupted);


#if defined(USE_OMP) && defined(_OPENMP)
#pragma omp atomic
#endif
            numdone++;
#endif //SILENT


            /* Calculate over all ngb cells */
            const cellarray_index *first  = &(lattice1[index1]);
            if(first->nelements == 0) {
                continue;
            }
            DOUBLE *x1 = X1 + first->start;
            DOUBLE *y1 = Y1 + first->start;
            DOUBLE *z1 = Z1 + first->start;
            const int64_t N1 = first->nelements;
            if(autocorr == 1) {
                int same_cell = 1;
                countpairs_rp_pi_driver(x1, y1, z1, N1,
                                        x1, y1, z1, N1,
                                        same_cell
                                        
#ifdef PERIODIC
                                        ,ZERO, ZERO, ZERO
#endif                                  
                                        ,sqr_rpmax, sqr_rpmin, nrpbin, npibin, rupp_sqr, pimax
#ifdef OUTPUT_RPAVG
                                        ,rpavg
#endif
                                        ,npairs);

            }
            for(int64_t ngb=0;ngb<first->num_ngb;ngb++){
                const cellarray_index *second = first->ngb_cells[ngb];
                if(second->nelements == 0) {
                    continue;
                }
                const int same_cell = 0;
                DOUBLE *x2 = X2 + second->start;
                DOUBLE *y2 = Y2 + second->start;
                DOUBLE *z2 = Z2 + second->start;
#ifdef PERIODIC
                const DOUBLE off_xwrap = first->xwrap[ngb];
                const DOUBLE off_ywrap = first->ywrap[ngb];
                const DOUBLE off_zwrap = first->zwrap[ngb];
#endif
                
                const int64_t N2 = second->nelements;

                countpairs_rp_pi_driver(x1, y1, z1, N1,
                                        x2, y2, z2, N2,
                                        same_cell
#ifdef PERIODIC
                                        ,off_xwrap, off_ywrap, off_zwrap
#endif                                  
                                        ,sqr_rpmax, sqr_rpmin, nrpbin, npibin, rupp_sqr, pimax
#ifdef OUTPUT_RPAVG
                                        ,rpavg
#endif
                                        ,npairs);

            }//loop over ngb cells
        }//index1 loop over totncells
#if defined(USE_OMP) && defined(_OPENMP)
        for(int i=0;i<totnbins;i++) {
            all_npairs[tid][i] = npairs[i];
#ifdef OUTPUT_RPAVG
            all_rpavg[tid][i] = rpavg[i];
#endif
        }
    }//close the omp parallel region
#endif

#ifndef SILENT    
    finish_myprogressbar(&interrupted);
#endif


#if defined(USE_OMP) && defined(_OPENMP)
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


    //The code does not double count for autocorrelations
    //which means the npairs and rpavg values need to be doubled;
    if(autocorr == 1) {
        const uint64_t int_fac = 2;
#ifdef OUTPUT_RPAVG
        const DOUBLE dbl_fac = (DOUBLE) 2.0;
#endif
        for(int i=0;i<totnbins;i++) {
            npairs[i] *= int_fac;
#ifdef OUTPUT_RPAVG
            rpavg[i] *= dbl_fac;
#endif
        }
    }

    
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
            assert(index < totnbins && "index must be within range");
            results->npairs[index] = npairs[index];
#ifdef OUTPUT_RPAVG
            results->rpavg[index] = rpavg[index];
#else
            results->rpavg[index] = 0.0;
#endif
        }
    }

    free(rupp);
    const int free_wraps = periodic == 1 ? 1:0;
    free_cellarray_index(lattice1,totncells, free_wraps);
    if(autocorr == 0) {
        free(lattice2);
    }
    
#if defined(USE_OMP) && defined(_OPENMP)
    matrix_free((void **) all_npairs, numthreads);
#ifdef OUTPUT_RPAVG
    matrix_free((void **) all_rpavg, numthreads);
#endif
#endif

    return results;
}
