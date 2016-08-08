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

int wp_driver(DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
              DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell, 
              const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
              const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
              ,DOUBLE *src_rpavg
              ,const struct config_options *options
              ,uint64_t *src_npairs)
{
	//Seriously this is the declaration for the function pointers...here be dragons.
	int (*allfunctions[]) ( DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                            DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell,
                            const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                            const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
                            ,DOUBLE *src_rpavg
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
    const int num_functions = sizeof(allfunctions)/sizeof(void *);
    
    //the fastest available code will always be at index 0.
    int function_dispatch=0;
    if(options->instruction_set != 0) {
        switch(options->instruction_set) {
        case(AVX):function_dispatch=0;break;
        case(SSE):function_dispatch=1;break;
        default:function_dispatch=2;break;
        }
    }
    if(function_dispatch >= num_functions) {
        fprintf(stderr,"ERROR: Could not resolve the correct function.\n Function index = %d must lie between [0, %d)\n",
                function_dispatch, num_functions);
    }
    
    //Now call the appropriate function -> look, runtime dispatch !
    int status = (allfunctions[function_dispatch])(x0, y0, z0, N0,
                                               x1, y1, z1, N1, same_cell,
                                               sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax,
                                               off_xwrap, off_ywrap, off_zwrap
                                               ,src_rpavg
                                               ,src_npairs);
    return status;
}

