/* File: wp_driver.h */
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

#include "defs.h"
#include "function_precision.h"
#include <inttypes.h>

    int wp_driver(DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                  DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell,
                  const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE *rupp_sqr, const DOUBLE pimax,
                  const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
                  ,DOUBLE *src_rpavg
                  ,const struct config_options *options
                  ,uint64_t *src_npairs);
    
#ifdef __cplusplus
}
#endif
