/* File: wp_kernels.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h"
#include <inttypes.h>

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif

void same_cell_wp_driver(const cellarray_index * first,
                         DOUBLE * restrict const X,DOUBLE * restrict const Y, DOUBLE * restrict const Z, 
                         const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax
#ifdef OUTPUT_RPAVG
                         ,DOUBLE *src_rpavg
#endif                         
                         ,uint64_t *src_npairs);


    
void diff_cells_wp_driver(const cellarray_index * first, const cellarray_index *second,
                          DOUBLE * restrict const X,DOUBLE * restrict const Y,DOUBLE * restrict const Z,
                          const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                          const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#ifdef OUTPUT_RPAVG
                          ,DOUBLE *src_rpavg
#endif                         
                          ,uint64_t *src_npairs);
    
#ifdef __cplusplus
}
#endif
