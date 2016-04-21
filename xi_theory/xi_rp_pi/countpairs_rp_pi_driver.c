/* File: countpairs_rp_pi_driver.c */
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

#include "defs.h"//minimal header macros that might be useful generally
#include "countpairs_rp_pi_driver.h"//the function prototypes + function_precision.h

//directly include the kernel file with
//actual implementations. The appropriate
//AVX or SSE header files are included in the
//kernel files

#include "countpairs_rp_pi_kernels.c"

void countpairs_rp_pi_driver(DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                             DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1,
                             const int same_cell
#ifdef PERIODIC
                             ,const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#endif                        
                             ,const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const int npibin, const DOUBLE *rupp_sqr, const DOUBLE pimax
                             
#ifdef OUTPUT_RPAVG
                             ,DOUBLE *src_rpavg
#endif
                             ,uint64_t *src_npairs)
{

#if defined(USE_AVX) && defined(__AVX__)
    /*----------------- AVX --------------------*/
    AVX_FLOATS m_rupp_sqr[nbin];
    for(int i=0;i<nbin;i++) {
        m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
    }
    AVX_FLOATS m_kbin[nbin];
    for(int i=0;i<nbin;i++) {
        m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
    }

    //AVX is available and wanted -> call the AVX function
    countpairs_rp_pi_avx_intrinsics(x0, y0, z0, N0,
                                    x1, y1, z1, N1, same_cell, 
                                    sqr_rpmax, sqr_rpmin, nbin, npibin, rupp_sqr, pimax
#ifdef PERIODIC                              
                                    ,off_xwrap, off_ywrap, off_zwrap
#endif                              
                                    ,m_rupp_sqr
                                    ,m_kbin
#ifdef OUTPUT_RPAVG
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
    SSE_FLOATS m_kbin[nbin];
    for(int i=0;i<nbin;i++) {
        m_kbin[i] = SSE_SET_FLOAT((DOUBLE) i);
    }

    countpairs_rp_pi_sse_intrinsics(x0, y0, z0, N0,
                                    x1, y1, z1, N1, same_cell,
                                    sqr_rpmax, sqr_rpmin,  nbin, npibin, rupp_sqr, pimax
#ifdef PERIODIC
                                    ,off_xwrap, off_ywrap, off_zwrap
#endif                              
                                    ,m_rupp_sqr
                                    ,m_kbin
#ifdef OUTPUT_RPAVG
                                    ,src_rpavg
#endif                         
                                    ,src_npairs);
    //appropriate SSE function called -> return;
    return;
    /*----------------- END OF SSE --------------------*/
#else

    /*----------------- FALLBACK CODE --------------------*/
    const int64_t totnbins = (npibin+1)*(nbin+1);
    uint64_t npairs[totnbins];
    for(int i=0;i<totnbins;i++) {
        npairs[i]=0;
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[totnbins];
    for(int i=0;i<totnbins;i++) {
        rpavg[i]=0;
    }
#endif//OUTPUT_RPAVG

#warning USING NAIVE CODE

    const DOUBLE dpi = pimax/npibin;
    const DOUBLE inv_dpi = 1.0/dpi;
    
    /* naive implementation that is guaranteed to compile */
    for(int64_t i=0;i<N0;i++) {
#ifdef PERIODIC        
        const DOUBLE xpos = *x0++ + off_xwrap;
        const DOUBLE ypos = *y0++ + off_ywrap;
        const DOUBLE zpos = *z0++ + off_zwrap;
#else
        const DOUBLE xpos = *x0++;
        const DOUBLE ypos = *y0++;
        const DOUBLE zpos = *z0++;
#endif        

        DOUBLE *localz1 = z1;
        int64_t j= 0;
        if(same_cell == 1) {
            j = i+1;
            localz1 += j;
        } else {
            while(j < N1){
                const DOUBLE dz = *localz1++ - zpos;
                if(dz > -pimax) break;
                j++;
            }
            localz1--;
        }
        DOUBLE *localx1 = (DOUBLE *) x1 + j;
        DOUBLE *localy1 = (DOUBLE *) y1 + j;

        for(;j<N1;j++) {
            const DOUBLE dx = *localx1++ - xpos;
            const DOUBLE dy = *localy1++ - ypos;
            const DOUBLE dz = FABS((*localz1++ - zpos));
            if(dz >= pimax) continue;
            
            const DOUBLE r2 = dx*dx + dy*dy ;
            if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
#ifdef OUTPUT_RPAVG
            const DOUBLE r = SQRT(r2);
#endif            
            
            int pibin = (int) (dz*inv_dpi);
            pibin = pibin > npibin ? npibin:pibin;
            for(int kbin=nbin-1;kbin>=1;kbin--) {
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
    }
    for(int i=0;i<totnbins;i++) {
        src_npairs[i] += npairs[i];
#ifdef OUTPUT_RPAVG
        src_rpavg[i] += rpavg[i];
#endif        
    }
    /*----------------- FALLBACK CODE --------------------*/
    return;

#endif//SSE
#endif//AVX

}    
