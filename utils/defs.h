/* File: defs.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <math.h>

#define ADD_DIFF_TIME(t0,t1)     fabs((t1.tv_sec - t0.tv_sec) + 1e-6*(t1.tv_usec - t0.tv_usec))
#define ALIGNMENT                32

typedef enum {
    AVX=0,
    SSE,
    FALLBACK
} isa;//name for instruction sets

struct config_options
{
    size_t float_type; /* floating point type -> vectorized supports double/float; fallback can support long double*/
    int verbose; /* Outputs progressbar and times */

    int instruction_set; /* select instruction set to run on */
    int need_avg_sep; /* <rp> or <\theta> is required */
    int autocorr;
    
    /* Options for theory*/
    int periodic; /* count in periodic mode? flag ignored for wp/xi */
    int sort_on_z;/* option to sort particles based on their Z co-ordinate in gridlink*/
    double boxsize;

    /* Options for mocks */
    //cosmology struct. Intentionally left anoynoymous, so I can
    //directly access the fields.
    struct{
        double OMEGA_M;
        double OMEGA_B;
        double OMEGA_L;
        double HUBBLE;
        double LITTLE_H;
        double SIGMA_8;
        double NS;
    };
    
    /* the link_in_* variables control how the 3-D cell structure is created */
    int link_in_cz;/* Not implemented but present to enable brute-force calculation*/
    int link_ra; /* relevant for DDrppi_mocks and DDtheta_mocks.*/
    int link_in_dec;/* */

    /* Replaces the divide in DDrppi_mocks by a reciprocal and a Newton-Raphson step*/
    int fast_divide;
};

