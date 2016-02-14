/* File: countspheres_mocks.h */
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
#include "cellarray_mocks.h"
#include <inttypes.h> //for uint64_t


    //define the results structure
    typedef struct{
        DOUBLE **pN;
        DOUBLE rmax;
        int nbin;
        int nc;
        int num_pN;
    } results_countspheres_mocks;


    int count_neighbors(const DOUBLE xcen,const DOUBLE ycen,const DOUBLE zcen,const DOUBLE smin,const DOUBLE inv_rcube,const DOUBLE rmax,
                        const int ngrid,const cellarray *lattice, const int nthreshold, const int bin_refine_factor) __attribute__((warn_unused_result));

    results_countspheres_mocks * countspheres_mocks(const int64_t Ngal, DOUBLE *xgal, DOUBLE *ygal, DOUBLE *zgal,
                                                    const int64_t Nran, DOUBLE *xran, DOUBLE *yran, DOUBLE *zran,
                                                    const int threshold_neighbors,
                                                    const DOUBLE rmax, const int nbin, const int nc,
                                                    const int num_pN,
                                                    const char *centers_file,
                                                    const int cosmology)  __attribute__((warn_unused_result));


    void free_results_countspheres_mocks(results_countspheres_mocks **results);

#ifdef __cplusplus
}
#endif
