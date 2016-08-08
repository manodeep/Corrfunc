/* File: countpairs_wp.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

#include "defs.h"
#include "function_precision.h"
#include "cellarray.h" //definition of struct cellarray_index
#include "utils.h" //all of the utilities
#include "gridlink.h"//function proto-type for gridlink
#include "countpairs_wp.h" //function proto-type

#include "wp_driver.h" //function dispatch for the kernels
#include "progressbar.h" //for the progressbar

#if defined(USE_OMP) && defined(_OPENMP)
#include <omp.h>
#endif


void free_results_wp(results_countpairs_wp *results)
{
    if(results == NULL)
        return;

    free(results->npairs);
    free(results->rupp);
    free(results->wp);
    free(results->rpavg);
}



int countpairs_wp(const int64_t ND, DOUBLE * restrict X, DOUBLE * restrict Y, DOUBLE * restrict Z,
                  const double boxsize,
                  const int numthreads,
                  const char *binfile,
                  const double pimax,
                  results_countpairs_wp *results,
                  const struct config_options *options)
{
    if( ! (options->float_type == sizeof(float) || options->float_type == sizeof(double))){
        fprintf(stderr,"ERROR: In %s> Can only handle doubles or floats. Got an array of size = %zu\n",
                __FUNCTION__, options->float_type);
        return EXIT_FAILURE;
    }
    
    if(options->float_type != sizeof(DOUBLE)) {
        fprintf(stderr,"ERROR: In %s> Can only handle arrays of size=%zu. Got an array of size = %zu\n",
                __FUNCTION__, sizeof(DOUBLE), options->float_type);
        return EXIT_FAILURE;
    }

    int bin_refine_factor=2,zbin_refine_factor=1;
    int nmesh_x, nmesh_y, nmesh_z;
    
    /***********************
     *initializing the  bins
     ************************/
    double *rupp;
    double rpmin,rpmax;
    int nrpbins;
    setup_bins(binfile,&rpmin,&rpmax,&nrpbins,&rupp);
    assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
    assert(nrpbins > 0 && "Number of rp bins is valid");

    const DOUBLE xmin = 0.0, xmax=boxsize;
    const DOUBLE ymin = 0.0, ymax=boxsize;
    const DOUBLE zmin = 0.0, zmax=boxsize;
    DOUBLE rupp_sqr[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        rupp_sqr[i] = rupp[i]*rupp[i];
    }

    if(rpmax < 0.05*boxsize) bin_refine_factor = 1;
    
    const DOUBLE sqr_rpmin = rupp_sqr[0];
    const DOUBLE sqr_rpmax = rupp_sqr[nrpbins-1];

    //set up the 3-d grid structure. Each element of the structure contains a
    //pointer to the cellarray structure that itself contains all the points
    cellarray_index *lattice = gridlink_index(ND, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
        fprintf(stderr,"countpairs_wp> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n");
        bin_refine_factor *=2;
        zbin_refine_factor *=2;
        free(lattice);
        lattice = gridlink_index(ND, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    }
    const int64_t totncells = nmesh_x*nmesh_y*(int64_t) nmesh_z;
    const int periodic = 1;
    const int autocorr = 1;
    assign_ngb_cells(lattice, lattice, totncells, bin_refine_factor, bin_refine_factor, zbin_refine_factor, nmesh_x, nmesh_y, nmesh_z, boxsize, boxsize, boxsize, autocorr, periodic);

#if defined(USE_OMP) && defined(_OPENMP)
    //openmp specific constructs
    omp_set_num_threads(numthreads);
    uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, nrpbins);
    DOUBLE **all_rpavg;
    if(options->need_avg_sep) {
        all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE), numthreads, nrpbins);
    }

#else//sequential mode follows
    uint64_t npair[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        npair[i]=0;
    }
    DOUBLE rpavg[nrpbins];
    if(options->need_avg_sep) {
        for(int i=0;i<nrpbins;i++) {
            rpavg[i]=0;
        }
    }
#endif// USE_OMP


    int interrupted=0;
    int64_t numdone=0;
    if(options->verbose) {
        init_my_progressbar(totncells,&interrupted);
    }

    
#if defined(USE_OMP) && defined(_OPENMP)
#pragma omp parallel shared(numdone)
    {
        const int tid = omp_get_thread_num();
        uint64_t npair[nrpbins];
        for(int i=0;i<nrpbins;i++) {
            npair[i]=0;
        }
        DOUBLE rpavg[nrpbins];
        if(options->need_avg_sep) {
            for(int i=0;i<nrpbins;i++) {
                rpavg[i]=0.0;
            }
        }


#pragma omp for schedule(dynamic) nowait
#endif//USE_OMP
        for(int index1=0;index1<totncells;index1++) {


            if(options->verbose) {
#if defined(USE_OMP) && defined(_OPENMP)
                if (omp_get_thread_num() == 0)
#endif
                    my_progressbar(numdone,&interrupted);
                
                
#if defined(USE_OMP) && defined(_OPENMP)
#pragma omp atomic
#endif
                numdone++;
            }

            
            /* First do the same-cell calculations */
            const cellarray_index *first  = &(lattice[index1]);
            if(first->nelements == 0) {
                continue;
            }

            int same_cell = 1;
            DOUBLE *x1 = X + first->start;
            DOUBLE *y1 = Y + first->start;
            DOUBLE *z1 = Z + first->start;
            const int64_t N1 = first->nelements;
            int status;
            if(options->need_avg_sep) {
                status = wp_driver(x1, y1, z1, N1,
                                   x1, y1, z1, N1, same_cell,
                                   sqr_rpmax, sqr_rpmin, nrpbins, rupp_sqr, pimax,
                                   ZERO, ZERO, ZERO
                                   ,rpavg
                                   ,options
                                   ,npair);
            } else {
                status = wp_driver(x1, y1, z1, N1,
                                   x1, y1, z1, N1, same_cell,
                                   sqr_rpmax, sqr_rpmin, nrpbins, rupp_sqr, pimax,
                                   ZERO, ZERO, ZERO
                                   ,NULL
                                   ,options
                                   ,npair);
            }
            if(status != EXIT_SUCCESS) {
                exit(status);
            }
                
            for(int64_t ngb=0;ngb<first->num_ngb;ngb++){
                cellarray_index *second = first->ngb_cells[ngb];
                DOUBLE *x2 = X + second->start;
                DOUBLE *y2 = Y + second->start;
                DOUBLE *z2 = Z + second->start;
                const int64_t N2 = second->nelements;
                const DOUBLE off_xwrap = first->xwrap[ngb];
                const DOUBLE off_ywrap = first->ywrap[ngb];
                const DOUBLE off_zwrap = first->zwrap[ngb];
                same_cell = 0;
                if(options->need_avg_sep) {
                    status = wp_driver(x1, y1, z1, N1,
                                       x2, y2, z2, N2, same_cell,
                                       sqr_rpmax, sqr_rpmin, nrpbins, rupp_sqr, pimax,
                                       off_xwrap, off_ywrap, off_zwrap
                                       ,rpavg
                                       ,options
                                       ,npair);
                } else {
                    status = wp_driver(x1, y1, z1, N1,
                                       x2, y2, z2, N2, same_cell,
                                       sqr_rpmax, sqr_rpmin, nrpbins, rupp_sqr, pimax,
                                       off_xwrap, off_ywrap, off_zwrap
                                       ,NULL
                                       ,options
                                       ,npair);
                }
                if(status != EXIT_SUCCESS) {
                    exit(status);
                }
            }//ngb loop
        }//index1 loop
#if defined(USE_OMP) && defined(_OPENMP)
        for(int j=0;j<nrpbins;j++) {
            all_npairs[tid][j] = npair[j];
            if(options->need_avg_sep) {
                all_rpavg[tid][j] = rpavg[j];
            }
        }
    }//omp parallel
#endif

    if(options->verbose) {
        finish_myprogressbar(&interrupted);
    }


#if defined(USE_OMP) && defined(_OPENMP)
    uint64_t npair[nrpbins];
    DOUBLE rpavg[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        npair[i] = 0;
        if(options->need_avg_sep) {
            rpavg[i] = ZERO;
        }
    }

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<nrpbins;j++) {
            npair[j] += all_npairs[i][j];
            if(options->need_avg_sep) {
                rpavg[j] += all_rpavg[i][j];
            }
        }
    }
    matrix_free((void **) all_npairs,numthreads);
    if(options->need_avg_sep) {
        matrix_free((void **) all_rpavg, numthreads);
    }
#endif//USE_OMP

    if(options->need_avg_sep) {
        for(int i=0;i<nrpbins;i++) {
            if(npair[i] > 0) {
                rpavg[i] /= (DOUBLE) npair[i];
            }
        }
    }

    const int free_wraps = 1;
    free_cellarray_index(lattice, totncells, free_wraps);

    //Pack in the results
    results->nbin  = nrpbins;
    results->pimax = pimax;
    results->npairs = my_malloc(sizeof(uint64_t), nrpbins);
    results->wp = my_malloc(sizeof(DOUBLE), nrpbins);
    results->rupp   = my_malloc(sizeof(DOUBLE), nrpbins);
    results->rpavg  = my_malloc(sizeof(DOUBLE), nrpbins);

    const DOUBLE avgweight2 = 1.0, avgweight1 = 1.0;
    const DOUBLE density=0.5*avgweight2*ND/(boxsize*boxsize*boxsize);//pairs are not double-counted
    DOUBLE rlow=0.0 ;
    DOUBLE prefac_density_DD=avgweight1*ND*density;
    DOUBLE twice_pimax = 2.0*pimax;

    for(int i=0;i<nrpbins;i++) {
        results->npairs[i] = npair[i];
        results->rupp[i] = rupp[i];
        if(options->need_avg_sep) {
            results->rpavg[i] = rpavg[i];
        } else {
            results->rpavg[i] = 0.0;
        }
        const DOUBLE weight0 = (DOUBLE) results->npairs[i];
        /* compute xi, dividing summed weight by that expected for a random set */
        const DOUBLE vol=M_PI*(results->rupp[i]*results->rupp[i]-rlow*rlow)*twice_pimax;
        const DOUBLE weightrandom = prefac_density_DD*vol;
        results->wp[i] = (weight0/weightrandom-1)*twice_pimax;
        rlow=results->rupp[i];
    }
    free(rupp);
    return EXIT_SUCCESS;
}
