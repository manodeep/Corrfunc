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

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t


    //define the results structure
    typedef struct{
        uint64_t *npairs;
        DOUBLE *theta_upp;
        DOUBLE *theta_avg;
        int nbin;
    } results_countpairs_theta;


    results_countpairs_theta * countpairs_theta_mocks(const int64_t ND1, DOUBLE *phi1, DOUBLE *theta1,
                                                      const int64_t ND2, DOUBLE *phi2, DOUBLE *theta2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                      const int numthreads,
#endif
                                                      const int autocorr,
                                                      const char *binfile) __attribute__((warn_unused_result));

    results_countpairs_theta * countpairs_theta_brute_force(const int64_t ND1, DOUBLE *phi1, DOUBLE *theta1,
                                                            const int64_t ND2, DOUBLE *phi2, DOUBLE *theta2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                            const int numthreads,
#endif
                                                            const int autocorr,
                                                            const char *binfile) __attribute__((warn_unused_result));


    void free_results_countpairs_theta(results_countpairs_theta **results);
    void check_ra_dec(const int64_t N, DOUBLE *phi, DOUBLE *theta);

#ifdef __cplusplus
}
#endif
