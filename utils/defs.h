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
#include <stdbool.h>

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


#define OPTIONS_HEADER_SIZE     (1024)
struct config_options
{
    /* The fields should appear here in decreasing order of 
       alignment requirements. Generally speaking, alignment
       is at least the sizeof the variable type. double has
       8 byte alignment, int has 4 bytes, char has 1 byte etc...
       (size_t could be 4 or 8 bytes depending on compilation 
       mode)
     */

    /* Theory option for periodic boundaries */
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

    size_t float_type; /* floating point type -> vectorized supports double/float; fallback can support long double*/
    int instruction_set; /* select instruction set to run on */

    char version[32];/* fill in the version number */

    bool verbose; /* Outputs progressbar and times */
    bool need_avg_sep; /* <rp> or <\theta> is required */
    bool autocorr;
    
    /* Options for theory*/
    bool periodic; /* count in periodic mode? flag ignored for wp/xi */
    bool sort_on_z;/* option to sort particles based on their Z co-ordinate in gridlink*/

    /* For DDrppi_mocks and vpf, flag to indicate cz is already co-moving distance */
    bool is_comoving_dist;
    
    /* the link_in_* variables control how the 3-D cell structure is created */
    bool link_in_dec;/* relevant for DDthteta_mocks */
    bool link_in_ra; /* relevant for DDtheta_mocks.*/

    /* Replaces the divide in DDrppi_mocks in AVX mode by a reciprocal and a Newton-Raphson step. */
    bool fast_divide;

    /* Fast arccos for wtheta (effective only when OUTPUT_THETAAVG is enabled) */
    bool fast_acos;

    /* Reserving to maboolain ABI compatibility for the future */
    /* Note that the math here assumes no padding bytes, that's because of the 
       order in which the fields are declared (largest to smallest alignments)  */
    uint8_t reserved[OPTIONS_HEADER_SIZE - 32*sizeof(char) - sizeof(size_t) - 8*sizeof(double) - sizeof(int) - 11*sizeof(bool)];
};

/* Taken from http://stackoverflow.com/questions/19403233/compile-time-struct-size-check-error-out-if-odd 
   which is in turn taken from the linux kernel */
/* #define BUILD_BUG_OR_ZERO(e) (sizeof(struct{ int:-!!(e);})) //gives me unused value warning */
/* #define ENSURE_STRUCT_SIZE(e, size)  BUILD_BUG_OR_ZERO(sizeof(e) != size) */

#define BUILD_BUG_OR_ZERO(cond) typedef char assertion_on_mystruct[( !!(cond) )*2-1 ] 
#define ENSURE_STRUCT_SIZE(e, size)  BUILD_BUG_OR_ZERO(sizeof(e) == size)
