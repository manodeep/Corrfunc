/* File: countpairs_s_mu.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "countpairs_s_mu.h" //function proto-type for API
#include "countpairs_s_mu_impl_double.h"//actual implementations for double
#include "countpairs_s_mu_impl_float.h"//actual implementations for float

void free_results_s_mu(results_countpairs_s_mu *results)
{
    if(results==NULL)
        return;

    free(results->npairs);results->npairs = NULL;
    free(results->supp);results->supp = NULL;
    free(results->savg);results->savg = NULL;
    free(results->weightavg);results->weightavg = NULL;

    results->mu_max = 0.0;
    results->mu_min = 0.0;
    results->nsbin = 0;
    results->nmu_bins = 0;
}


int countpairs_s_mu(const int64_t ND1, void *X1, void *Y1, void *Z1,
                    const int64_t ND2, void *X2, void *Y2, void *Z2,
                    const int numthreads,
                    const int autocorr,
                    const char *sbinfile,
                    const double mu_max,
                    const int nmu_bins,
                    results_countpairs_s_mu *results,
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
        return countpairs_s_mu_float(ND1, (float *) X1, (float *) Y1, (float *) Z1,
                                     ND2, (float *) X2, (float *) Y2, (float *) Z2,
                                     numthreads,
                                     autocorr,
                                     sbinfile,
                                     mu_max,
                                     nmu_bins,
                                     results,
                                     options,
                                     extra);
    } else {
        return countpairs_s_mu_double(ND1, (double *) X1, (double *) Y1, (double *) Z1,
                                      ND2, (double *) X2, (double *) Y2, (double *) Z2,
                                      numthreads,
                                      autocorr,
                                      sbinfile,
                                      mu_max,
                                      nmu_bins,
                                      results,
                                      options,
                                      extra);
    }
}
