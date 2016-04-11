/* File: wp_kernels.c */
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
#include "wp_kernels.h"//the function prototypes

//directly include the kernel files with
//actual implementation. 
#include "same_cell_wp_kernels.c"
#include "diff_cells_wp_kernels.c"

void same_cell_wp_driver(const cellarray_index * first,
                         DOUBLE * restrict const X,DOUBLE * restrict const Y, DOUBLE * restrict const Z, 
                         const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax
#ifdef OUTPUT_RPAVG
                         ,DOUBLE *src_rpavg
#endif                         
                         ,uint64_t *src_npairs)
{
    DOUBLE *x0 = X + first->start;
    DOUBLE *y0 = Y + first->start;
    DOUBLE *z0 = Z + first->start;
    const int64_t N0 = first->nelements;

#if defined(USE_AVX) && defined(__AVX__)
    /*----------------- AVX --------------------*/
    AVX_FLOATS m_rupp_sqr[nbin];
    for(int i=0;i<nbin;i++) {
        m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
    }
#ifdef OUTPUT_RPAVG
    AVX_FLOATS m_kbin[nbin];
    for(int i=0;i<nbin;i++) {
        m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
    }
#endif//RPAVG + AVX

    //AVX is available and wanted -> call the AVX function
    same_cell_avx_intrinsics(x0, y0, z0, N0,
                             sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax
                             ,m_rupp_sqr
#ifdef OUTPUT_RPAVG
                             ,m_kbin
                             ,src_rpavg
#endif                         
                             ,src_npairs);
    
    //Already dispatched to appropriate AVX function -> return;
    return;
    /*----------------- END OF AVX --------------------*/
#else //AVX

#if defined (__SSE4_2__)

    /*------------------ SSE -------------------*/
    SSE_FLOATS m_rupp_sqr[nbin];
    for(int i=0;i<nbin;i++) {
        m_rupp_sqr[i] = SSE_SET_FLOAT(rupp_sqr[i]);
    }
#ifdef OUTPUT_RPAVG
    SSE_FLOATS m_kbin[nbin];
    for(int i=0;i<nbin;i++) {
        m_kbin[i] = SSE_SET_FLOAT((DOUBLE) i);
    }
#endif//RPAVG + SSE
    same_cell_sse_intrinsics(x0, y0, z0, N0,
                             sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax
                             ,m_rupp_sqr
#ifdef OUTPUT_RPAVG
                             ,m_kbin
                             ,src_rpavg
#endif                         
                             ,src_npairs);
    //appropriate SSE function called -> return;
    return;
    /*----------------- END OF SSE --------------------*/
#else

    /*----------------- FALLBACK CODE --------------------*/
    uint64_t npairs[nbin];
    for(int i=0;i<nbin;i++) {
        npairs[i]=0;
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nbin];
    for(int i=0;i<nbin;i++) {
        rpavg[i]=0;
    }
#endif//OUTPUT_RPAVG

    /* naive implementation that is guaranteed to compile */
    for(int64_t i=0;i<N0;i++) {
        const DOUBLE xpos = *x0++;
        const DOUBLE ypos = *y0++;
        const DOUBLE zpos = *z0++;

        DOUBLE *x1 = x0;
        DOUBLE *y1 = y0;
        DOUBLE *z1 = z0;
        for(int64_t j=i+1;j<N0;j++) {
            const DOUBLE dx = *x1++ - xpos;
            const DOUBLE dy = *y1++ - ypos;
            const DOUBLE dz = *z1++ - zpos;
            if(dz >=pimax) break;

            const DOUBLE r2 = dx*dx + dy*dy;
            if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
#ifdef OUTPUT_RPAVG
            const DOUBLE r = SQRT(r2);
#endif
            for(int kbin=nbin-1;kbin>=1;kbin--){
                if(r2 >= rupp_sqr[kbin-1]) {
                    npairs[kbin]++;
#ifdef OUTPUT_RPAVG
                    rpavg[kbin] += r;
#endif                    
                    break;
                }
            }//searching for kbin loop                                                               
        }
    }

    for(int i=0;i<nbin;i++) {
        src_npairs[i] += npairs[i];
#ifdef OUTPUT_RPAVG
        src_rpavg[i] += rpavg[i];
#endif        
    }
    /*----------------- FALLBACK CODE --------------------*/

#endif//SSE
#endif//AVX

}

void diff_cells_wp_driver(const cellarray_index * first, const cellarray_index *second, 
                          DOUBLE * restrict const X,DOUBLE * restrict const Y, DOUBLE * restrict const Z, 
                          const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                          const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#ifdef OUTPUT_RPAVG
                          ,DOUBLE *src_rpavg
#endif                         
                          ,uint64_t *src_npairs)
{

    DOUBLE *x0 = X + first->start;
    DOUBLE *y0 = Y + first->start;
    DOUBLE *z0 = Z + first->start;
    const int64_t N0 = first->nelements;
    
    DOUBLE *x1 = X + second->start;
    DOUBLE *y1 = Y + second->start;
    DOUBLE *z1 = Z + second->start;
    const int64_t N1 = second->nelements;
    
#if defined(USE_AVX) && defined(__AVX__)
    /*----------------- AVX --------------------*/
    AVX_FLOATS m_rupp_sqr[nbin];
    for(int i=0;i<nbin;i++) {
        m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
    }
#ifdef OUTPUT_RPAVG
    AVX_FLOATS m_kbin[nbin];
    for(int i=0;i<nbin;i++) {
        m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
    }
#endif//RPAVG + AVX

    //AVX is available and wanted -> call the AVX function
    diff_cells_avx_intrinsics(x0, y0, z0, N0,
                              x1, y1, z1, N1, 
                              sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax,
                              off_xwrap, off_ywrap, off_zwrap
                              ,m_rupp_sqr
#ifdef OUTPUT_RPAVG
                              ,m_kbin
                              ,src_rpavg
#endif                         
                              ,src_npairs);
    
    //Already dispatched to appropriate AVX function -> return;
    return;
    /*----------------- END OF AVX --------------------*/
#else //AVX

#if defined (__SSE4_2__)

    /*------------------ SSE -------------------*/
    SSE_FLOATS m_rupp_sqr[nbin];
    for(int i=0;i<nbin;i++) {
        m_rupp_sqr[i] = SSE_SET_FLOAT(rupp_sqr[i]);
    }
#ifdef OUTPUT_RPAVG
    SSE_FLOATS m_kbin[nbin];
    for(int i=0;i<nbin;i++) {
        m_kbin[i] = SSE_SET_FLOAT((DOUBLE) i);
    }
#endif//RPAVG + SSE
    diff_cells_sse_intrinsics(x0, y0, z0, N0,
                              x1, y1, z1, N1,
                              sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax,
                              off_xwrap, off_ywrap, off_zwrap
                             ,m_rupp_sqr
#ifdef OUTPUT_RPAVG
                             ,m_kbin
                             ,src_rpavg
#endif                         
                             ,src_npairs);
    //appropriate SSE function called -> return;
    return;
    /*----------------- END OF SSE --------------------*/
#else

    /*----------------- FALLBACK CODE --------------------*/
    uint64_t npairs[nbin];
    for(int i=0;i<nbin;i++) {
        npairs[i]=0;
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nbin];
    for(int i=0;i<nbin;i++) {
        rpavg[i]=0;
    }
#endif//OUTPUT_RPAVG

    /* naive implementation that is guaranteed to compile */
    for(int64_t i=0;i<N0;i++) {
        const DOUBLE xpos = *x0++ + off_xwrap;
        const DOUBLE ypos = *y0++ + off_ywrap;
        const DOUBLE zpos = *z0++ + off_zwrap;

        DOUBLE *localx1 = (DOUBLE *) x1; 
        DOUBLE *localy1 = (DOUBLE *) y1;
        DOUBLE *localz1 = (DOUBLE *) z1; 
        
        int64_t j=0;
        while(j < N1){
            const DOUBLE dz = *localz1++ - zpos;
            if(dz > -pimax) break;
            j++;
        }
        localz1 -= 1;
        localx1 += j;
        localy1 += j;
        
        for(;j<N1;j++) {
            const DOUBLE dx = *localx1++ - xpos;
            const DOUBLE dy = *localy1++ - ypos;
            const DOUBLE dz = *localz1++ - zpos;
            if(dz >=pimax) break;

            const DOUBLE r2 = dx*dx + dy*dy;
            if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
#ifdef OUTPUT_RPAVG
            const DOUBLE r = SQRT(r2);
#endif
            for(int kbin=nbin-1;kbin>=1;kbin--){
                if(r2 >= rupp_sqr[kbin-1]) {
                    npairs[kbin]++;
#ifdef OUTPUT_RPAVG
                    rpavg[kbin] += r;
#endif                    
                    break;
                }
            }//searching for kbin loop                                                               
        }
    }

    for(int i=0;i<nbin;i++) {
        src_npairs[i] += npairs[i];
#ifdef OUTPUT_RPAVG
        src_rpavg[i] += rpavg[i];
#endif        
    }
    /*----------------- FALLBACK CODE --------------------*/

#endif//SSE
#endif//AVX

}




