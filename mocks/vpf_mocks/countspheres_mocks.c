/* File: countspheres_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "countspheres_mocks.h"//for definition of struct results_countspheres_mocks

#include "countspheres_mocks_impl_float.h"//for the actual implementations of float/double
#include "countspheres_mocks_impl_double.h"
#include "utils.h"//for matrix_free

void free_results_countspheres_mocks(results_countspheres_mocks *results)
{
    if(results == NULL)
        return;

    matrix_free((void **) results->pN, results->nbin);
    results->pN = NULL;
    results->rmax = 0.0;
    results->nbin = 0;
    results->nc = 0;
    results->num_pN = 0;
}



int countspheres_mocks(const int64_t Ngal, void *xgal, void *ygal, void *zgal,
                       const int64_t Nran, void * restrict xran, void * restrict yran, void * restrict zran,
                       const int threshold_neighbors,
                       const double rmax, const int nbin, const int nc,
                       const int num_pN,
                       const char *centers_file,
                       const int cosmology,
                       results_countspheres_mocks *results,
                       struct config_options *options,
                       struct extra_options *extra)
{

    if( ! (options->float_type == sizeof(float) || options->float_type == sizeof(double))){
        fprintf(stderr,"ERROR: In %s> Can only handle doubles or floats. Got an array of size = %zu\n",
                __FUNCTION__, options->float_type);
        return EXIT_FAILURE;
    }

    if( strncmp(options->version, STR(VERSION), sizeof(options->version)/sizeof(char)-1) != 0) {
        fprintf(stderr,"Error: Do not know this API version = `%s'. Expected version = `%s'\n", options->version, STR(VERSION));
        return EXIT_FAILURE;
    }

    if(options->float_type == sizeof(float)) {
        return countspheres_mocks_float(Ngal, (float *) xgal, (float *)ygal, (float *)zgal,
                                        Nran, (float *) xran, (float *) yran, (float *) zran,
                                        threshold_neighbors,
                                        rmax, nbin, nc,
                                        num_pN,
                                        centers_file,
                                        cosmology,
                                        results,
                                        options,
                                        extra);
    } else {
        return countspheres_mocks_double(Ngal, (double *) xgal, (double *)ygal, (double *)zgal,
                                         Nran, (double *) xran, (double *) yran, (double *) zran,
                                         threshold_neighbors,
                                         rmax, nbin, nc,
                                         num_pN,
                                         centers_file,
                                         cosmology,
                                         results,
                                         options,
                                         extra);
    }
}
