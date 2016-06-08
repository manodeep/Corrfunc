/* File: countpairs_rp_pi.h */
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

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t

    //define the results structure
    typedef struct{
        uint64_t *npairs;
        DOUBLE *rupp;
        DOUBLE *rpavg;
        DOUBLE pimax;
        int nbin;
        int npibin;
    } results_countpairs_rp_pi;

    results_countpairs_rp_pi * countpairs_rp_pi(const int64_t ND1, DOUBLE *X1, DOUBLE *Y1, DOUBLE *Z1,
                                                const int64_t ND2, DOUBLE *X2, DOUBLE *Y2, DOUBLE *Z2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                const int numthreads,
#endif
                                                const int autocorr,
                                                const char *binfile,
                                                const DOUBLE pimax)  __attribute__((warn_unused_result));

    void free_results_rp_pi(results_countpairs_rp_pi **results);

#ifdef __cplusplus
}
#endif
