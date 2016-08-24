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

#define STRINGIFY(x)   #x
#define STR(x) STRINGIFY(x)

typedef enum {
    FALLBACK=0, /* No special options */
    SSE=1,  /* 64 bit vectors */
    SSE2=2, /* 128 bit vectors */
    SSE3=3, /* 128 bit vectors */
    SSSE3=4, /* 128 bit vectors */
    SSE4=5,/* 128bit vectors */
    SSE42=6, /* 128bit vectors with blend operations */
    AVX=7, /* 256bit vector width */
    AVX2=8,  /* AVX2 (integer operations)*/
    AVX512F=9,/* AVX 512 Foundation */
    NUM_ISA 
} isa;//name for instruction sets -> corresponds to the return from instrset_detect in cpu_features.c

struct config_options
{
    char version[32];/* fill in the version number */
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
    int link_in_dec;/* relevant for DDthteta_mocks */
    int link_in_ra; /* relevant for DDtheta_mocks.*/

    /* Replaces the divide in DDrppi_mocks in AVX mode by a reciprocal and a Newton-Raphson step. */
    int fast_divide;

    /* Fast arccos for wtheta (effective only when OUTPUT_THETAAVG is enabled) */
    int fast_acos;

    /* Reserving to maintain ABI compatibility for the future */
    int reserved[32];
};

