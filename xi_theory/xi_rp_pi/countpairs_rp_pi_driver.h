/* File: countpairs_rp_pi_driver.h */
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
                             ,uint64_t *src_npairs);

    
#ifdef __cplusplus
}
#endif
