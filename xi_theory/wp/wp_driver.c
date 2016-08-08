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
#include <stdint.h>
#include <inttypes.h>

#include "defs.h"
#include "function_precision.h"

#include "cpu_features.h"//for instrset_detect
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
    static int initialized = 0;
    static int (*function)( DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                            DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell,
                            const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                            const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
                            ,DOUBLE *src_rpavg
                            ,uint64_t *src_npairs) = NULL;
        
    //the fastest available code will always be at index 0.
    if(initialized == 0) {
        //Seriously this is the declaration for the function pointers...here be dragons.
        int (*allfunctions[]) ( DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                                DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell,
                                const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                                const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
                                ,DOUBLE *src_rpavg
                                ,uint64_t *src_npairs)
            = {
#ifdef __AVX__
            wp_avx_intrinsics,
#endif			 
#ifdef __SSE4_2__
            wp_sse_intrinsics,
#endif
            wp_fallback
        };
        const int num_functions = sizeof(allfunctions)/sizeof(void *);
        const int fallback_offset = num_functions - 1;
        const int highest_isa = instrset_detect();

        int curr_offset = 0;
        
        /* Now check if AVX is supported by the CPU */
        int avx_offset = fallback_offset;
#ifdef __AVX__
        avx_offset = highest_isa >= 7 ? curr_offset:fallback_offset;
        curr_offset++;
#endif        

        /* Is the SSE function supported at runtime and enabled at compile-time?*/
        int sse_offset = fallback_offset;
#ifdef __SSE4_2__
        sse_offset = highest_isa >= 6 ? curr_offset:fallback_offset;
        curr_offset++;
#endif
        if( curr_offset != fallback_offset) {
            fprintf(stderr,"ERROR: Bug in code (current offset = %d *should equal* fallback function offset = %d)\n",
                    curr_offset, fallback_offset);
            return EXIT_FAILURE;
        } 
        
        int function_dispatch=0;
        /* Check that cpu supports feature */
        if(options->instruction_set != 0) {
            switch(options->instruction_set) {
            case(AVX):
                function_dispatch=avx_offset;break;
            case(SSE):function_dispatch=sse_offset;break;
            default:function_dispatch=fallback_offset;break;
            }
        }
        if(function_dispatch >= num_functions) {
            fprintf(stderr,"In %s> ERROR: Could not resolve the correct function.\n Function index = %d must lie between [0, %d)\n",
                    __FUNCTION__, function_dispatch, num_functions);
            return EXIT_FAILURE;
        }
        function = allfunctions[function_dispatch];
        initialized = 1;
    }
    if(function == NULL) {
        fprintf(stderr,"In %s> ERROR: function to be dispatched can not be NULL", __FUNCTION__);
        return EXIT_FAILURE;
    } 
    
    //Now call the appropriate function -> look, runtime dispatch !
    int status = function(x0, y0, z0, N0,
                          x1, y1, z1, N1, same_cell,
                          sqr_rpmax, sqr_rpmin,  nbin, rupp_sqr, pimax,
                          off_xwrap, off_ywrap, off_zwrap
                          ,src_rpavg
                          ,src_npairs);
    return status;
}

