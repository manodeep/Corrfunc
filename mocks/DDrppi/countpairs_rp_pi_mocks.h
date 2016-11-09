/* File: countpairs_rp_pi_mocks.h */
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
        double *rupp;
        double *rpavg;
        double *weightavg;
        double pimax;
        int nbin;
        int npibin;
    } results_countpairs_mocks;
    
    int countpairs_mocks(const int64_t ND1, void *theta1, void *phi1, void *czD1,
                         const int64_t ND2, void *theta2, void *phi2, void *czD2,
                         const int numthreads,
                         const int autocorr,
                         const char *binfile,
                         const double pimax,
                         const int cosmology,
                         results_countpairs_mocks *results,
                         struct config_options *options, struct extra_options *extra);

    void free_results_mocks(results_countpairs_mocks *results);

#ifdef __cplusplus
}
#endif
