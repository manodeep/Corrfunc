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

#define MAX_NUM_WEIGHTS 10

typedef struct
{
    void *weights[MAX_NUM_WEIGHTS];  // This will be of shape weights[num_weights][num_particles]
    int64_t num_weights;
} weight_struct;

typedef enum {
  NONE=-42, /* default */
  PAIR_PRODUCT=0,
  NUM_WEIGHT_TYPE 
} weight_method_t; // type of weighting to apply

/* Gives the number of weight arrays required by the given weighting method
 */
static inline int get_num_weights_by_method(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return 1;
        default:
        case NONE:
            return 0;
    }
}

/* Maps a name to weighting method
   `method` will be set on return.
 */
static inline int get_weight_method_by_name(const char *name, weight_method_t *method){
    if(name == NULL || strcmp(name, "") == 0){
        *method = NONE;
        return EXIT_SUCCESS;
    }
    if(strcmp(name, "pair_product") == 0 || strcmp(name, "p") == 0){
        *method = PAIR_PRODUCT;
        return EXIT_SUCCESS;
    }
        
    return EXIT_FAILURE;
}
    
struct extra_options
{
    // Two possible weight_structs (at most we will have two loaded sets of particles)
    weight_struct weights0;
    weight_struct weights1;
    weight_method_t weight_method; // the function that will get called to give the weight of a particle pair
    uint8_t reserved[EXTRA_OPTIONS_HEADER_SIZE - 2*sizeof(weight_struct) - sizeof(weight_method_t)];
};

// Here we want to return an int because malloc may fail (unlike get_config_options)
// weight_method determines the number of various weighting arrays that we allocate
static inline int get_extra_options(struct extra_options *extra, const weight_method_t weight_method)
{    
    ENSURE_STRUCT_SIZE(struct extra_options, EXTRA_OPTIONS_HEADER_SIZE);//compile-time check for making sure struct is correct size
    if(extra == NULL) {
        return EXIT_FAILURE;
    }

    memset(extra, 0, EXTRA_OPTIONS_HEADER_SIZE);
    extra->weight_method = weight_method;
    
    weight_struct *w0 = &(extra->weights0);
    weight_struct *w1 = &(extra->weights1);
    w0->num_weights = get_num_weights_by_method(extra->weight_method);
    w1->num_weights = w0->num_weights;
    
    if(w0->weights == NULL || w1->weights == NULL) {
        free(w0->weights); free(w1->weights);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static inline void free_extra_options(struct extra_options *extra)
{
    weight_struct *w0 = &(extra->weights0);
    weight_struct *w1 = &(extra->weights1);
    
    // Free particle lists
    for(int64_t i=0; i < w0->num_weights; i++) {
        free(w0->weights[i]);
        w0->weights[i] = NULL;  // avoid double free from aliased w0/w1
    }
    w0->num_weights = 0;
  
    for(int64_t i=0; i < w1->num_weights; i++) {
        free(w1->weights[i]);
        w1->weights[i] = NULL;
    }
    w1->num_weights = 0;
}    


#include "macros.h"

    
#ifdef __cplusplus
}
#endif
