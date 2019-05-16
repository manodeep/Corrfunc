/* File: countpairs_theta_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "countpairs_theta_mocks.h" //function proto-type for API
#include "countpairs_theta_mocks_impl_double.h"//actual implementations for double
#include "countpairs_theta_mocks_impl_float.h"//actual implementations for float

void free_results_countpairs_theta(results_countpairs_theta *results)
{
    if(results == NULL)
        return;

    free(results->theta_upp);results->theta_upp = NULL;
    free(results->npairs);results->npairs = NULL;
    free(results->theta_avg);results->theta_avg = NULL;
    free(results->weightavg);results->weightavg = NULL;

    results->nbin = 0;
}



int countpairs_theta_mocks(const int64_t ND1, void *phi1, void *theta1,
                           const int64_t ND2, void *phi2, void *theta2,
                           const int numthreads,
                           const int autocorr,
                           const char *binfile,
                           results_countpairs_theta *results,
                           struct config_options *options, struct extra_options *extra)
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
        return countpairs_theta_mocks_float(ND1, (float *) phi1, (float *) theta1,
                                            ND2, (float *) phi2, (float *) theta2,
                                            numthreads,
                                            autocorr,
                                            binfile,
                                            results,
                                            options,
                                            extra);
    } else {
        return countpairs_theta_mocks_double(ND1, (double *) phi1, (double *) theta1,
                                             ND2, (double *) phi2, (double *) theta2,
                                             numthreads,
                                             autocorr,
                                             binfile,
                                             results,
                                             options,
                                             extra);
    }
}
