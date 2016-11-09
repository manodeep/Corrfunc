/* File: countpairs_xi.h */
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

#include "defs.h" //for struct config_options
#include <stdint.h> //for uint64_t


//define the results structure
    typedef struct{
        uint64_t *npairs;
        double *xi;
        double *rupp;
        double *ravg;
        double *weightavg;
        int nbin;
    } results_countpairs_xi;

    extern void free_results_xi(results_countpairs_xi *results);
    extern int countpairs_xi(const int64_t ND1, void * restrict X1, void * restrict Y1, void * restrict Z1,
                             const double boxsize,
                             const int numthreads,
                             const char *binfile,
                             results_countpairs_xi *results,
                             struct config_options *options,
                             struct extra_options *extra);
    
#ifdef __cplusplus
}
#endif
