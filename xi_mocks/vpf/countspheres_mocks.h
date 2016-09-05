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

#include "defs.h"//for definition of struct config_options
#include <inttypes.h> //for uint64_t


    //define the results structure
    typedef struct{
        double **pN;
        double rmax;
        int nbin;
        int nc;
        int num_pN;
    } results_countspheres_mocks;

    extern int countspheres_mocks(const int64_t Ngal, void *xgal, void *ygal, void *zgal,
                                  const int64_t Nran, void * restrict xran, void * restrict yran, void * restrict zran,
                                  const int threshold_neighbors,
                                  const double rmax, const int nbin, const int nc,
                                  const int num_pN,
                                  const char *centers_file,
                                  const int cosmology,
                                  results_countspheres_mocks *results,
                                  struct config_options *options, struct extra_options *extra);
    
    extern void free_results_countspheres_mocks(results_countspheres_mocks *results);

#ifdef __cplusplus
}
#endif
