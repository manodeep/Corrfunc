/* File: wp_driver.c */
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
#include "wp_driver.h"//the function prototypes

//directly include the kernel files with
//actual implementation. 
#include "wp_kernels.c"

void wp_driver(DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
               DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell, 
               const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
               const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#ifdef OUTPUT_RPAVG
               ,DOUBLE *src_rpavg
#endif                         
               ,uint64_t *src_npairs)
{
    const int unroll_fac=4;
    if(N1 < unroll_fac*NVEC) {
        wp_fallback(x0, y0, z0, N0,
                    x1, y1, z1, N1, same_cell,
                    sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax,
                    off_xwrap, off_ywrap, off_zwrap
#ifdef OUTPUT_RPAVG
                    ,src_rpavg
#endif                         
                    ,src_npairs);
        return;
    }


	//Seriously this is the declaration for the function pointers...here be dragons.
	void (*allfunctions[]) ( DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                             DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell,
                             const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                             const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#ifdef OUTPUT_RPAVG
                             ,DOUBLE *src_rpavg
#endif                         
                             ,uint64_t *src_npairs)
		= {
#if defined(__AVX__) && defined(USE_AVX)		   
        wp_avx_intrinsics,
#endif			 
#ifdef __SSE4_2__
        wp_sse_intrinsics,
#endif
        wp_fallback
	};
    //the fastest available code will always be at index 0.
    const int func_dispatch = 0;
    
/* #if (!(defined(__AVX__) && defined(USE_AVX))) || !defined(__SSE4_2__) */
/*     //Old hardware probably -> no AVX or SSE -> use fallback */
/*     const int func_dispatch = 0; // */
/* #else     */
/*     //At least AVX or SSE is available. */
/* #if defined(__AVX__) && defined(USE_AVX) */
/*     //call AVX for > 2*NVEC -> unroll AVX loop by 2 */
/*     const int func_dispatch = 0; */
/* #else */
/* #if defined(__SSE4_2__) */
/*     const int func_dispatch = 0; */
/* #endif//SSE4.2     */
/* #endif//__AVX__ */
/* #endif//__AVX__ || SSE */

    //Now call the appropriate function -> look, runtime dispatch !
    (allfunctions[func_dispatch])(x0, y0, z0, N0,
                                  x1, y1, z1, N1, same_cell,
                                  sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax,
                                  off_xwrap, off_ywrap, off_zwrap
#ifdef OUTPUT_RPAVG
                                  ,src_rpavg
#endif                         
                                  ,src_npairs);
    
    
}

