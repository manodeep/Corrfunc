/* File: countpairs.h */
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

            
#include "defs.h"//for struct config_options 
#include <stdint.h> //for uint64_t

//define the results structure
  typedef struct{
    uint64_t *npairs;
    double *rupp;
    double *rpavg;
    double *weightavg;
    int nbin;
  } results_countpairs;
  
  extern int countpairs(const int64_t ND1, void *X1, void *Y1, void  *Z1,
                        const int64_t ND2, void *X2, void *Y2, void  *Z2,
                        const int numthreads,
                        const int autocorr,
                        const char *binfile,
                        results_countpairs *results,
                        struct config_options *options,
                        struct extra_options *extra) __attribute__((warn_unused_result));
  
  extern void free_results(results_countpairs *results);

#ifdef __cplusplus
}
#endif
