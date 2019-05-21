/* File: countpairs_wp.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "countpairs_wp.h" //function proto-type for API
#include "countpairs_wp_impl_double.h"//actual implementations for double
#include "countpairs_wp_impl_float.h"//actual implementations for float

void free_results_wp(results_countpairs_wp *results)
{
    if(results == NULL)
        return;

    free(results->rupp);results->rupp = NULL;
    free(results->npairs);results->npairs = NULL;
    free(results->rpavg);results->rpavg = NULL;
    free(results->weightavg);results->weightavg = NULL;
    free(results->wp);results->wp = NULL;

    results->nbin = 0;
    results->pimax = 0.0;
}

int countpairs_wp(const int64_t ND, void * restrict X, void * restrict Y, void * restrict Z,
                  const double boxsize,
                  const int numthreads,
                  const char *binfile,
                  const double pimax,
                  results_countpairs_wp *results,
                  struct config_options *options,
                  struct extra_options *extra)
{
    if( ! (options->float_type == sizeof(float) || options->float_type == sizeof(double))){
        fprintf(stderr,"ERROR: In %s> Can only handle doubles or floats. Got an array of size = %zu\n",
                __FUNCTION__, options->float_type);
        return EXIT_FAILURE;
    }

    if( strncmp(options->version, STR(VERSION), sizeof(options->version)/sizeof(char)-1 ) != 0) {
        fprintf(stderr,"Error: Do not know this API version = `%s'. Expected version = `%s'\n", options->version, STR(VERSION));
        return EXIT_FAILURE;
    }

    if(options->float_type == sizeof(float)) {
      return countpairs_wp_float(ND, (float *) X, (float *) Y, (float *) Z,
                                 boxsize,
                                 numthreads,
                                 binfile,
                                 pimax,
                                 results,
                                 options,
                                 extra);
    } else {
      return countpairs_wp_double(ND, (double *) X, (double *) Y, (double *) Z,
                                  boxsize,
                                  numthreads,
                                  binfile,
                                  pimax,
                                  results,
                                  options,
                                  extra);
    }
}
