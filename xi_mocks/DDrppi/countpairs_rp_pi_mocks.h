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
    } results_countpairs_mocks;

    results_countpairs_mocks * countpairs_mocks(const int64_t ND1, DOUBLE *theta1, DOUBLE *phi1, DOUBLE *czD1,
                                                const int64_t ND2, DOUBLE *theta2, DOUBLE *phi2, DOUBLE *czD2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                const int numthreads,
#endif
                                                const int autocorr,
                                                const char *binfile,
                                                const DOUBLE pimax,
                                                const int cosmology) __attribute__((warn_unused_result));

    void check_ra_dec_cz(const int64_t N, DOUBLE *phi, DOUBLE *theta, DOUBLE *cz);

    void free_results_mocks(results_countpairs_mocks **results);

#ifdef __cplusplus
}
#endif
