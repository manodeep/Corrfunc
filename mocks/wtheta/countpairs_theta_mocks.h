/* File: countpairs_theta_mocks.h */
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
#include <stdint.h> //for uint64_t

    //define the results structure
    typedef struct{
        uint64_t *npairs;
        double *theta_upp;
        double *theta_avg;
        double *weightavg;
        int nbin;
    } results_countpairs_theta;

    extern int countpairs_theta_mocks(const int64_t ND1, void *phi1, void *theta1,
                                      const int64_t ND2, void *phi2, void *theta2,
                                      const int numthreads,
                                      const int autocorr,
                                      const char *binfile,
                                      results_countpairs_theta *results,
                                      struct config_options *options, struct extra_options *extra);
    
    extern void free_results_countpairs_theta(results_countpairs_theta *results);

#ifdef __cplusplus
}
#endif
