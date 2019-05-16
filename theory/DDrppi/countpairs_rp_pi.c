/* File: countpairs_rp_pi.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "countpairs_rp_pi.h" //function proto-type for API
#include "countpairs_rp_pi_impl_double.h"//actual implementations for double
#include "countpairs_rp_pi_impl_float.h"//actual implementations for float

void free_results_rp_pi(results_countpairs_rp_pi *results)
{
    if(results==NULL)
        return;

    free(results->rupp);results->rupp = NULL;
    free(results->npairs);results->npairs = NULL;
    free(results->rpavg);results->rpavg = NULL;
    free(results->weightavg);results->weightavg = NULL;

    results->nbin = 0;
    results->npibin = 0;
    results->pimax = 0.0;
}


int countpairs_rp_pi(const int64_t ND1, void *X1, void *Y1, void *Z1,
                     const int64_t ND2, void *X2, void *Y2, void *Z2,
                     const int numthreads,
                     const int autocorr,
                     const char *binfile,
                     const double pimax,
                     results_countpairs_rp_pi *results,
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
      return countpairs_rp_pi_float(ND1, (float *) X1, (float *) Y1, (float *) Z1,
                                    ND2, (float *) X2, (float *) Y2, (float *) Z2,
                                    numthreads,
                                    autocorr,
                                    binfile,
                                    pimax,
                                    results,
                                    options,
                                    extra);
    } else {
        return countpairs_rp_pi_double(ND1, (double *) X1, (double *) Y1, (double *) Z1,
                                       ND2, (double *) X2, (double *) Y2, (double *) Z2,
                                       numthreads,
                                       autocorr,
                                       binfile,
                                       pimax,
                                       results,
                                       options,
                                       extra);
    }
}
