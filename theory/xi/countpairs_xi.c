/* File: countpairs_xi.c */
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

#include "countpairs_xi.h" //function proto-type for API
#include "countpairs_xi_impl_double.h"//actual implementations for double
#include "countpairs_xi_impl_float.h"//actual implementations for float

void free_results_xi(results_countpairs_xi *results)
{
    if(results == NULL)
        return;

    free(results->rupp);
    free(results->xi);
    free(results->npairs);
    free(results->ravg);
    free(results->weightavg);
}


int countpairs_xi(const int64_t ND, void * restrict X, void * restrict Y, void * restrict Z,
                  const double boxsize,
                  const int numthreads,
                  const char *binfile,
                  results_countpairs_xi *results,
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
        return countpairs_xi_float(ND, (float * restrict) X, (float * restrict) Y, (float * restrict) Z,
                                   boxsize,
                                   numthreads,
                                   binfile,
                                   results,
                                   options,
                                   extra);
    } else {
        return countpairs_xi_double(ND, (double * restrict) X, (double * restrict) Y, (double * restrict) Z,
                                    boxsize,
                                    numthreads,
                                    binfile,
                                    results,
                                    options,
                                    extra);
    }
}
