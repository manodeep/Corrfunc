/* File: countpairs_wp.h */
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
#include <stdint.h>
    
    //define the results structure
    typedef struct{
        uint64_t *npairs;
        double *wp;
        double *rupp;
        double *rpavg;
        double *weightavg;
        double pimax;
        int nbin;
    } results_countpairs_wp;
    
    extern int countpairs_wp(const int64_t ND1, void * restrict X1, void * restrict Y1, void * restrict Z1,
                             const double boxsize,
                             const int numthreads,
                             const char *binfile,
                             const double pimax,
                             results_countpairs_wp *result,
                             struct config_options *options,
                             struct extra_options *extra) __attribute__((warn_unused_result));

    extern void free_results_wp(results_countpairs_wp *results);

#ifdef __cplusplus
}
#endif
