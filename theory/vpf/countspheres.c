/* File: countspheres.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "utils.h" //for matrix_free
#include "countspheres.h" //API Header

//Actual implementations
#include "countspheres_impl_float.h"
#include "countspheres_impl_double.h"


void free_results_countspheres(results_countspheres *results)
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



int countspheres(const int64_t np, void * restrict X, void * restrict Y, void * restrict Z,
                 const double rmax, const int nbin, const int nc,
                 const int num_pN,
                 unsigned long seed,
                 results_countspheres *results,
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
        return countspheres_float(np,  (float *) X,  (float *) Y, (float *) Z,
                                  rmax, nbin, nc,
                                  num_pN,
                                  seed,
                                  results,
                                  options,
                                  extra);
    } else {
        return countspheres_double(np,  (double *) X, (double *) Y, (double *) Z,
                                   rmax, nbin, nc,
                                   num_pN,
                                   seed,
                                   results,
                                   options,
                                   extra);
    }

}
