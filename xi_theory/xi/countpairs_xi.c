/* File: countpairs_xi.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "gridlink.h"//function proto-type for gridlink
#include "countpairs_xi.h" //function proto-type
#include "cellarray.h" //definition of struct cellarray
#include "utils.h" //all of the utilities

#include "xi_driver.h"//driver header for xi

#ifndef SILENT
#include "progressbar.h" //for the progressbar
#endif

#if defined(USE_OMP) && defined(_OPENMP)
#include <omp.h>
#endif

#ifndef PERIODIC
#warning "xi is only valid for PERIODIC boundary conditions. Ignoring the Makefile (non)-definition of PERIODIC"
#endif


void free_results_xi(results_countpairs_xi *results)
{
    if(results == NULL)
        return;

    free(results->rupp);
    free(results->xi);
    free(results->npairs);
    free(results->rpavg);
}


results_countpairs_xi countpairs_xi(const int64_t ND, DOUBLE * restrict X, DOUBLE * restrict Y, DOUBLE * restrict Z,
                                    const double boxsize,
#if defined(USE_OMP) && defined(_OPENMP)
                                    const int numthreads,
#endif
                                    const char *binfile)
{
    //How many bins to subdivide rmax into -> affects runtime on O(20-30%) levels.
    //Check with your typical use-case and set appropriately. Values of 1,2 and 3 are
    //all you might need to check.
    int bin_refine_factor=1;

    /***********************
     *initializing the  bins
     ************************/
    double *rupp;
    int nbins;
    double rpmin,rpmax;
    setup_bins(binfile,&rpmin,&rpmax,&nbins,&rupp);
    assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
    assert(nbins > 0 && "Number of rp bins must be > 0");

    /*---Create 3-D lattice--------------------------------------*/
    int nmesh_x=0,nmesh_y=0,nmesh_z=0;
    const DOUBLE xmin = 0.0, xmax=boxsize;
    const DOUBLE ymin = 0.0, ymax=boxsize;
    const DOUBLE zmin = 0.0, zmax=boxsize;

    cellarray_index_particles *lattice = gridlink_index_particles(ND, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
        fprintf(stderr,"%s> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n",__FUNCTION__);
        bin_refine_factor *=2;
        free(lattice);
        lattice = gridlink_index_particles(ND, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, rpmax, bin_refine_factor, bin_refine_factor, bin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    }
    const int64_t totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;
    const int autocorr = 1;
    const int periodic = 1;
    assign_ngb_cells_index_particles(lattice, lattice, totncells, bin_refine_factor, bin_refine_factor, bin_refine_factor, nmesh_x, nmesh_y, nmesh_z, boxsize, boxsize, boxsize, autocorr, periodic);

#if !(defined(USE_OMP) && defined(_OPENMP))
    uint64_t npairs[nbins];
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nbins];
#endif
    for(int i=0; i < nbins;i++) {
        npairs[i] = 0;
#ifdef OUTPUT_RPAVG
        rpavg[i] = 0.0;
#endif
    }

#else
    omp_set_num_threads(numthreads);
    uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, nbins);
#ifdef OUTPUT_RPAVG
    DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE),numthreads,nbins);
#endif
#endif

    DOUBLE rupp_sqr[nbins];
    for(int i=0; i < nbins;i++) {
        rupp_sqr[i] = rupp[i]*rupp[i];
    }


    const DOUBLE pimax = rpmax;
    const DOUBLE sqr_rpmax=rupp_sqr[nbins-1];
    const DOUBLE sqr_rpmin=rupp_sqr[0];


#ifndef SILENT
    int interrupted=0;
    int64_t numdone=0;
    init_my_progressbar(totncells,&interrupted);
#endif    

    /*---Loop-over-Data1-particles--------------------*/
#if defined(USE_OMP) && defined(_OPENMP)
#ifndef SILENT
#pragma omp parallel shared(numdone)
#else
#pragma omp parallel    
#endif//SILENT    
    {
        const int tid = omp_get_thread_num();
        uint64_t npairs[nbins];
#ifdef OUTPUT_RPAVG
        DOUBLE rpavg[nbins];
#endif
        for(int i=0;i<nbins;i++) {
            npairs[i] = 0;
#ifdef OUTPUT_RPAVG
            rpavg[i] = 0.0;
#endif
        }

#pragma omp for  schedule(dynamic) nowait
#endif
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
#endif//SILENT

            /* First do the same-cell calculations */
            const cellarray_index_particles *first  = &(lattice[index1]);
            if(first->nelements == 0) {
                continue;
            } 

            DOUBLE *x1 = first->x;
            DOUBLE *y1 = first->y;
            DOUBLE *z1 = first->z;
            const int64_t N1 = first->nelements;
            int same_cell = 1;
            xi_driver(x1, y1, z1, N1,
                      x1, y1, z1, N1, same_cell, 
                      sqr_rpmax, sqr_rpmin, nbins, rupp_sqr, pimax,
                      ZERO, ZERO, ZERO
#ifdef OUTPUT_RPAVG
                      ,rpavg
#endif
                      ,npairs);
            
            for(int64_t ngb=0;ngb<first->num_ngb;ngb++){
                const cellarray_index_particles *second = first->ngb_cells[ngb];
                if(second->nelements == 0) {
                    continue;
                }
                DOUBLE *x2 = second->x;
                DOUBLE *y2 = second->y;
                DOUBLE *z2 = second->z;
                const int64_t N2 = second->nelements;
                const DOUBLE off_xwrap = first->xwrap[ngb];
                const DOUBLE off_ywrap = first->ywrap[ngb];
                const DOUBLE off_zwrap = first->zwrap[ngb];
                same_cell = 0;
                xi_driver(x1, y1, z1, N1,
                          x2, y2, z2, N2, same_cell,
                          sqr_rpmax, sqr_rpmin, nbins, rupp_sqr, pimax,
                          off_xwrap, off_ywrap, off_zwrap
#ifdef OUTPUT_RPAVG
                          ,rpavg
#endif
                          ,npairs);
            }//ngb loop
        }//index1 loop
#if defined(USE_OMP) && defined(_OPENMP)
        for(int j=0;j<nbins;j++) {
            all_npairs[tid][j] = npairs[j];
        }
#ifdef OUTPUT_RPAVG
        for(int j=0;j<nbins;j++) {
            all_rpavg[tid][j] = rpavg[j];
        }
#endif

    }//close the omp parallel region
#endif

#ifndef SILENT    
    finish_myprogressbar(&interrupted);
#endif

#if defined(USE_OMP) && defined(_OPENMP)
    uint64_t npairs[nbins];
    for(int i=0;i<nbins;i++) npairs[i] = 0;

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<nbins;j++) {
            npairs[j] += all_npairs[i][j];
        }
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nbins];
    for(int i=0;i<nbins;i++) rpavg[i] = 0.0;

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<nbins;j++) {
            rpavg[j] += all_rpavg[i][j];
        }
    }
    matrix_free((void **) all_rpavg, numthreads);
#endif//OUTPUT_RPAVG
    matrix_free((void **) all_npairs, numthreads);
#endif//USE_OMP


    //So the npairs array contains the number of pairs
    //and the rpavg array contain the *SUM* of separations
    //Let's divide out rpavg by npairs to actually get
    //the mean rpavg

#ifdef OUTPUT_RPAVG
    for(int i=0;i<nbins;i++) {
        if(npairs[i] > 0) {
            rpavg[i] /= (DOUBLE) npairs[i] ;
        }
    }
#endif

    //Pack in the results
    results_countpairs_xi results;
    results.nbin = nbins;
    results.npairs = my_malloc(sizeof(uint64_t), nbins);
    results.xi     = my_malloc(sizeof(DOUBLE)  , nbins);
    results.rupp   = my_malloc(sizeof(DOUBLE)  , nbins);
    results.rpavg  = my_malloc(sizeof(DOUBLE)  , nbins);

    const DOUBLE avgweight2 = 1.0, avgweight1 = 1.0;
    const DOUBLE density=0.5*avgweight2*ND/(boxsize*boxsize*boxsize);//0.5 because pairs are not double-counted
    const DOUBLE prefac_density=avgweight1*ND*density;

    DOUBLE rlow=0.0 ;
    //The first bin contains junk
    for(int i=0;i<nbins;i++) {
        results.npairs[i] = npairs[i];
        results.rupp[i]   = rupp[i];
#ifdef OUTPUT_RPAVG
        results.rpavg[i] = rpavg[i];
#else
        results.rpavg[i] = 0.0;
#endif

        const DOUBLE weight0 = (DOUBLE) results.npairs[i];
        const DOUBLE vol=4.0/3.0*M_PI*(rupp[i]*rupp[i]*rupp[i]-rlow*rlow*rlow);
        /* compute xi, dividing summed weight by that expected for a random set */
        const DOUBLE weightrandom = prefac_density*vol;
        assert(weightrandom > 0 && "Random weight is <= 0.0 - that is impossible");
        results.xi[i] = (weight0/weightrandom-1.0);
        rlow=results.rupp[i];
    }

    free_cellarray_index_particles(lattice, totncells);
    free(rupp);


    return results;

}
