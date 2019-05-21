/* File: countpairs_s_mu_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "countpairs_s_mu_mocks.h" //function proto-type for API
#include "countpairs_s_mu_mocks_impl_double.h"//actual implementations for double
#include "countpairs_s_mu_mocks_impl_float.h"//actual implementations for float

void free_results_mocks_s_mu(results_countpairs_mocks_s_mu *results)
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


int countpairs_mocks_s_mu(const int64_t ND1, void *phi1, void *theta1, void *czD1,
                          const int64_t ND2, void *phi2, void *theta2, void *czD2,
                          const int numthreads,
                          const int autocorr,
                          const char *sbinfile,
                          const double mu_max,
                          const int nmu_bins,
                          const int cosmology,
                          results_countpairs_mocks_s_mu *results,
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
        return countpairs_mocks_s_mu_float(ND1, (float *) phi1, (float *) theta1, (float *) czD1,
                                           ND2, (float *) phi2, (float *) theta2, (float *) czD2,
                                           numthreads,
                                           autocorr,
                                           sbinfile,
                                           mu_max,
                                           nmu_bins,
                                           cosmology,
                                           results,
                                           options,
                                           extra);
    } else {
        return countpairs_mocks_s_mu_double(ND1, (double *) phi1, (double *) theta1, (double *) czD1,
                                            ND2, (double *) phi2, (double *) theta2, (double *) czD2,
                                            numthreads,
                                            autocorr,
                                            sbinfile,
                                            mu_max,
                                            nmu_bins,
                                            cosmology,
                                            results,
                                            options,
                                            extra);
    }
}
