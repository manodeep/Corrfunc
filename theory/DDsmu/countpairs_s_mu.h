/* File: countpairs_s_mu.h */
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
        double *supp;
        double *savg;
        double mu_max;
        double mu_min;//not used -> assumed to be 0.0
        double *weightavg;
        int nsbin;
        int nmu_bins;
    } results_countpairs_s_mu;

    extern int countpairs_s_mu(const int64_t ND1, void *X1, void *Y1, void *Z1,
                               const int64_t ND2, void *X2, void *Y2, void *Z2,
                               const int numthreads,
                               const int autocorr,
                               const char *sbinfile,
                               const double mu_max,
                               const int nmu_bins, 
                               results_countpairs_s_mu *results,
                               struct config_options *options,
                               struct extra_options *extra);
    
    extern void free_results_s_mu(results_countpairs_s_mu *results);

#ifdef __cplusplus
}
#endif
