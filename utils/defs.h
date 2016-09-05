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
#include <string.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ADD_DIFF_TIME(t0,t1)     fabs((t1.tv_sec - t0.tv_sec) + 1e-6*(t1.tv_usec - t0.tv_usec))
#define ALIGNMENT                32

#define STRINGIFY(x)   #x
#define STR(x) STRINGIFY(x)

#define API_VERSION          STR("2.0.0")

typedef enum {
  DEFAULT=-42,/* present simply to make the enum a signed int*/
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

    /* Measures the time spent in the C API while accessed from python.
       Enabled when the flag c_timer is set
     */
    double c_api_time;
    
    size_t float_type; /* floating point type -> vectorized supports double/float; fallback can support long double*/
    int instruction_set; /* select instruction set to run on */

    char version[32];/* fill in the version number */
    uint8_t verbose; /* Outputs progressbar and times */
    uint8_t c_api_timer; /* Measures time spent in the C function */

    /* Options valid for both theory and mocks */
    uint8_t need_avg_sep; /* <rp> or <\theta> is required */
    uint8_t autocorr;/* Only one dataset is required */
    
    /* Options for theory*/
    uint8_t periodic; /* count in periodic mode? flag ignored for wp/xi */
    uint8_t sort_on_z;/* option to sort particles based on their Z co-ordinate in gridlink*/

    /* For DDrppi_mocks and vpf*/
    uint8_t is_comoving_dist;/* flag to indicate cz is already co-moving distance */
    
    /* the link_in_* variables control how the 3-D cell structure is created */
    uint8_t link_in_dec;/* relevant for DDthteta_mocks */
    uint8_t link_in_ra; /* relevant for DDtheta_mocks.*/

    /* Replaces the divide in DDrppi_mocks in AVX mode by a reciprocal and a Newton-Raphson step. */
    uint8_t fast_divide;

    /* Fast arccos for wtheta (effective only when OUTPUT_THETAAVG is enabled) */
    uint8_t fast_acos;

    /* Reserving to maintain ABI compatibility for the future */
    /* Note that the math here assumes no padding bytes, that's because of the 
       order in which the fields are declared (largest to smallest alignments)  */
    uint8_t reserved[OPTIONS_HEADER_SIZE - 33*sizeof(char) - sizeof(size_t) - 9*sizeof(double) - sizeof(int) - 11*sizeof(uint8_t)];
};

/* Taken from http://stackoverflow.com/questions/19403233/compile-time-struct-size-check-error-out-if-odd 
   which is in turn taken from the linux kernel */
/* #define BUILD_BUG_OR_ZERO(e) (sizeof(struct{ int:-!!(e);})) */
/* #define ENSURE_STRUCT_SIZE(e, size)  BUILD_BUG_OR_ZERO(sizeof(e) != size) */
/* However, the previous one gives me an unused-value warning and I do not want 
   to turn that compiler warning off. Therefore, this version, which results in 
   an unused local typedef warning is used. I turn off the corresponding warning 
   in common.mk (-Wno-unused-local-typedefs) via CFLAGS
*/
#define BUILD_BUG_OR_ZERO(cond, msg) typedef volatile char assertion_on_##msg[( !!(cond) )*2-1 ] 
#define ENSURE_STRUCT_SIZE(e, size)                 BUILD_BUG_OR_ZERO(sizeof(e) == size, sizeof_struct_config_options)

static inline struct config_options get_config_options(void)
{
    ENSURE_STRUCT_SIZE(struct config_options, OPTIONS_HEADER_SIZE);//compile-time check for making sure struct is correct size

    if(strncmp(API_VERSION, STR(VERSION), 32) != 0) {
        fprintf(stderr,"Error: Version mismatch between header and Makefile. Header claims version = `%s' while Makefile claims version = `%s'\n"
                "Library header probably needs to be updated\n", API_VERSION, STR(VERSION));
        exit(EXIT_FAILURE);
    }
    struct config_options options;
    memset(&options, 0, OPTIONS_HEADER_SIZE);
    snprintf(options.version, sizeof(options.version)/sizeof(char)-1, "%s", API_VERSION);
#ifdef DOUBLE_PREC    
    options.float_type = sizeof(double);
#else
    options.float_type = sizeof(float);
#endif    
#ifndef SILENT
    options.verbose = 1;
#endif
    
#ifdef OUTPUT_RPAVG    
    options.need_avg_sep = 1;
#endif
#ifdef PERIODIC
    options.periodic = 1;
#endif    

#ifdef __AVX__
    options.instruction_set = AVX;
#elif defined(__SSE4_2__)
    options.instruction_set = SSE42;
#else
    options.instruction_set = FALLBACK;
#endif

    /* Options specific to mocks */

    /* Options for DDrppi_mocks */
#ifdef FAST_DIVIDE
    options.fast_divide=1;
#endif

    
    /* Options for wtheta*/
#ifdef OUTPUT_THETAAVG
    options.need_avg_sep = 1;
#endif

#ifdef LINK_IN_DEC
    options.link_in_dec = 1;
#endif
#ifdef LINK_IN_RA
    options.link_in_ra=1;
    options.link_in_dec=1;
#endif

#ifdef FAST_ACOS
    options.fast_acos=1;
#endif    

#ifdef COMOVING_DIST
    options.is_comoving_dist=1;
#endif
    
    return options;
}


#define EXTRA_OPTIONS_HEADER_SIZE     (1024)
struct extra_options
{
    void **weights;//pointer to an array of pointers to store the weight arrays
    void (*weightfunc)(void);//Treacherous territory, generic weighting function pointer
    uint64_t num_weights;//number of valid weight arrays
    uint64_t weighting_func_type;//way to type-cast the generic weightfunc into the actual
                                //function. 
    
    uint8_t reserved[EXTRA_OPTIONS_HEADER_SIZE - sizeof(void **) - sizeof(void *) - 2*sizeof(uint64_t)];
};

static inline int get_extra_options(struct extra_options *extra)
{    
    ENSURE_STRUCT_SIZE(struct extra_options, EXTRA_OPTIONS_HEADER_SIZE);//compile-time check for making sure struct is correct size
    if(extra == NULL) {
        return EXIT_FAILURE;
    }

    memset(extra, 0, EXTRA_OPTIONS_HEADER_SIZE);
    /*Pre-allocate space for 2 sets of weights array pointers */
    extra->num_weights = 2;
    extra->weights = malloc(sizeof(*(extra->weights)) * extra->num_weights);
    if(extra->weights == NULL) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static inline void free_extra_options(struct extra_options *extra)
{
    for(uint64_t i=0;i<extra->num_weights;i++) {
        free(extra->weights[i]);
    }
    free(extra->weights);
    extra->weights = NULL;
    extra->num_weights = 0;
}    
    

#ifdef __cplusplus
}
#endif
