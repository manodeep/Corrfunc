/* File: defs.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

#include "macros.h"
#include "cpu_features.h"

#ifdef __cplusplus
extern "C" {
#endif

#define API_VERSION          STR("2.5.1")

/* Macros as mask for the binning_flags */
/* These consititute the 32 bytes for
the ``uint32_t binning_flags`` */

#define BINNING_REF_MASK         0x0000000F //Last 4 bits for how the bin sizes are calculated is done. Also indicates if refines are in place
#define BINNING_ORD_MASK         0x000000F0 //Next 4 bits for how the 3-D-> 1-D index conversion
/* The upper 24 bits are unused currently */

#define BINNING_DFL   0x0
#define BINNING_CUST  0x1

struct api_cell_timings
{
    int64_t N1;/* Number of points in the first cell*/
    int64_t N2;/* Number of points in the second cell */
    int64_t time_in_ns;/* Time taken in the compute kernel, measured in nano-seconds*/
    int first_cellindex;
    int second_cellindex;
    int tid;/* Thread-id, 0 for serial case, wastes 4 bytes, since thread id is 4bytes integer and not 8 bytes */
};


#define MAX_FAST_DIVIDE_NR_STEPS  3
#define OPTIONS_HEADER_SIZE     1024

#define BOXSIZE_NOTGIVEN (-2.)

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
    union {
        double boxsize;
        double boxsize_x;
    };
    double boxsize_y;
    double boxsize_z;

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

    /* Per cell timers. Keeps track of the number of particles per cell pair
       and time spent to compute the pairs. Might slow down code */
    struct api_cell_timings *cell_timings;
    int64_t totncells_timings;


    size_t float_type; /* floating point type -> vectorized supports double/float; fallback can support long double*/
    int32_t instruction_set; /* select instruction set to run on */

    char version[32];/* fill in the version number */
    uint8_t verbose; /* Outputs progressbar and times */
    uint8_t c_api_timer; /* Measures time spent in the C function */
    uint8_t c_cell_timer;/* Measures time spent per cell-pair. Might slow down the code */

    /* Options valid for both theory and mocks */
    uint8_t need_avg_sep; /* <rp> or <\theta> is required */
    uint8_t autocorr;/* Only one dataset is required */

    /* Options for theory*/
    uint8_t periodic; /* count in periodic mode? flag ignored for wp/xi */
    uint8_t sort_on_z;/* option to sort particles based on their Z co-ordinate in gridlink */

    /* For DDrppi_mocks and vpf*/
    uint8_t is_comoving_dist;/* flag to indicate cz is already co-moving distance */

    /* the link_in_* variables control how the 3-D cell structure is created */
    uint8_t link_in_dec;/* relevant for DDthteta_mocks */
    uint8_t link_in_ra; /* relevant for DDtheta_mocks.*/

    /* Replaces the divide in DDrppi_mocks in AVX mode by a reciprocal and a Newton-Raphson step. */
    uint8_t fast_divide_and_NR_steps;/* Used in AVX512/AVX; if set to 0, the standard (slow) divide is used
                                        If > 0, the value is interpreted as the number of NR steps
                                        i.e., fast_divide_and_NR_steps = 2, performs two steps of Newton-Raphson
                                        Anything greater than ~5, probably makes the code slower than the
                                        divide without any improvement in precision
                                      */


    /* Fast arccos for wtheta (effective only when OUTPUT_THETAAVG is enabled) */
    uint8_t fast_acos;

    /* Enabled by default */
    uint8_t enable_min_sep_opt;/* Whether to enable min. separation optimizations introduced in v2.3*/

    int8_t bin_refine_factors[3];/* Array for the custom bin refine factors in each dim
                                    xyz for theory routines and ra/dec/cz for mocks
                                    Must be signed integers since some for loops might use -bin_refine_factor
                                    as the starting point */

    uint16_t max_cells_per_dim;/* max number of cells per dimension. same for both theory and mocks */

    uint8_t copy_particles;/* whether to make a copy of the particle positions */
    uint8_t use_heap_sort;/* to allow using heap-sort instead of quicksort from sglib (relevant when the input particles are mostly sorted
                           and consequently quicksort becomes an O(N^2) process */
    union{
        uint32_t binning_flags;/* flag for all linking features,
                                  Will contain OR'ed flags from enum from `binning_scheme`
                                  Intentionally set as unsigned int, since in the
                                  future we might want to support some bit-wise OR'ed
                                  functionality */
        uint8_t bin_masks[4];
    };

    /* Reserving to maintain ABI compatibility for the future */
    uint8_t reserved[OPTIONS_HEADER_SIZE - 33*sizeof(char) - sizeof(size_t) - 11*sizeof(double) - 3*sizeof(int)
                     - sizeof(uint16_t) - 16*sizeof(uint8_t) - sizeof(struct api_cell_timings *) - sizeof(int64_t)
                     ];
};

static inline void set_bin_refine_scheme(struct config_options *options, const int8_t flag)
{
    //Set the top (nbits-4) to whatever already exists in binning_flag
    //and then set the bottom 4 bits to BIN_DFL
    options->binning_flags = (options->binning_flags & ~BINNING_REF_MASK) | (flag & BINNING_REF_MASK);
}


static inline void reset_bin_refine_scheme(struct config_options *options)
{
    set_bin_refine_scheme(options, BINNING_DFL);
}

static inline int8_t get_bin_refine_scheme(struct config_options *options)
{
    //Return the last 4 bits as 8 bits int
    return (int8_t) (options->binning_flags & BINNING_REF_MASK);
}

static inline void set_bin_refine_factors(struct config_options *options, const int bin_refine_factors[3])
{
    for(int i=0;i<3;i++) {
        int8_t bin_refine = bin_refine_factors[i];
        if(bin_refine_factors[i] > INT8_MAX) {
            fprintf(stderr,"Warning: bin refine factor[%d] can be at most %d. Found %d instead\n", i,
                    INT8_MAX, bin_refine_factors[i]);
            bin_refine = 1;
        }
        options->bin_refine_factors[i] = bin_refine;

    }
    /*
      Note, programmatically setting the refine factors resets the binning flag to "BINNING_DFL"
      BINNING_CUST is only set via function parameters, or explicitly */
    reset_bin_refine_scheme(options);
}

static inline void set_custom_bin_refine_factors(struct config_options *options, const int bin_refine_factors[3])
{
    set_bin_refine_factors(options, bin_refine_factors);
    set_bin_refine_scheme(options, BINNING_CUST);
}

static inline void reset_bin_refine_factors(struct config_options *options)
{
    /* refine factors of 2,2,1 in the xyz dims
       seems to produce the fastest code */
    options->bin_refine_factors[0] = 2;
    options->bin_refine_factors[1] = 2;
    options->bin_refine_factors[2] = 1;
    reset_bin_refine_scheme(options);
}



static inline void set_max_cells(struct config_options *options, const int max)
{
    if(max <=0) {
        fprintf(stderr,"Warning: Max. cells per dimension was requested to be set to "
                "a negative number = %d...returning\n", max);
        return;
    }

    if(max > INT16_MAX) {
        fprintf(stderr,"Warning: Max cells per dimension is a 2-byte integer and can not "
                "hold supplied value of %d. Max. allowed value for max_cells_per_dim is %d\n",
                max, INT16_MAX);
    }

    options->max_cells_per_dim = max;
}

static inline void reset_max_cells(struct config_options *options)
{
    options->max_cells_per_dim = NLATMAX;
}


static inline struct config_options get_config_options(void)
{
    ENSURE_STRUCT_SIZE(struct config_options, OPTIONS_HEADER_SIZE);//compile-time check for making sure struct is correct size

    if(strncmp(API_VERSION, STR(VERSION), 32) != 0) {
        fprintf(stderr,"Error: Version mismatch between header and Makefile. Header claims version = `%s' while Makefile claims version = `%s'\n"
                "Library header probably needs to be updated\n", API_VERSION, STR(VERSION));
        exit(EXIT_FAILURE);
    }
    struct config_options options;
    BUILD_BUG_OR_ZERO(sizeof(options.max_cells_per_dim) == sizeof(int16_t), max_cells_per_dim_must_be_16_bits);
    BUILD_BUG_OR_ZERO(sizeof(options.binning_flags) == sizeof(uint32_t), binning_flags_must_be_32_bits);
    BUILD_BUG_OR_ZERO(sizeof(options.bin_refine_factors[0]) == sizeof(int8_t), bin_refine_factors_must_be_8_bits);

    memset(&options, 0, OPTIONS_HEADER_SIZE);
    snprintf(options.version, sizeof(options.version)/sizeof(char)-1, "%s", API_VERSION);

    // If periodic, BOXSIZE_NOTGIVEN requires the user to set a boxsize.
    // A value of 0 will use automatic detection of the particle extent.
    // -1 makes that dimension non-periodic.
    options.boxsize_x = BOXSIZE_NOTGIVEN;
    options.boxsize_y = BOXSIZE_NOTGIVEN;
    options.boxsize_z = BOXSIZE_NOTGIVEN;

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

#ifdef __AVX512F__
    options.instruction_set = AVX512F;
#elif defined(__AVX2__)
    options.instruction_set = AVX2;
#elif defined(__AVX__)
    options.instruction_set = AVX;
#elif defined(__SSE4_2__)
    options.instruction_set = SSE42;
#else
    options.instruction_set = FALLBACK;
#endif

    /* Options specific to mocks */
    /* Options for DDrppi_mocks (FAST_DIVIDE is also applicable for both DDsmu, and DDsmu_mocks) */
#if defined(FAST_DIVIDE)
#if FAST_DIVIDE > MAX_FAST_DIVIDE_NR_STEPS
    options.fast_divide_and_NR_steps = MAX_FAST_DIVIDE_NR_STEPS;
#else
    options.fast_divide_and_NR_steps = FAST_DIVIDE;
#endif
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

#ifdef ENABLE_MIN_SEP_OPT
    //Introduced in Corrfunc v2.3
    options.enable_min_sep_opt=1;/* optimizations based on min. separation between cell-pairs. Enabled by default */
#endif

#ifdef FAST_ACOS
    options.fast_acos=1;
#endif

#ifdef COMOVING_DIST
    options.is_comoving_dist=1;
#endif

#ifdef COPY_PARTICLES
    /* Config options introduced in Corrfunc v2.3*/
    options.copy_particles = 1;/* make a copy of particles (positions and weights) (by default) */
#else
    // Using the input particles -> positions will have to re-ordered
    // Setting the next option will mean that the particles will be re-ordered
    // into their input order when the calculation completes. Usually relevant when
    // there are other "properties" arrays for the same particle; and changing the
    // positions would
    options.copy_particles = 0;
#endif //Create a copy of particle positions (doubles the memory usage)

    /* For the thread timings */
    options.totncells_timings = 0;
    /* If the API level timers are requested, then
       this pointer will have to be allocated */
    options.cell_timings = NULL;

    /*Setup the binning options */
    reset_max_cells(&options);
    reset_bin_refine_factors(&options);
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
    // These should not be strncmp because we want the implicit length comparison of strcmp.
    // It is still safe because one of the args is a string literal.
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

// weight_method determines the number of various weighting arrays that we allocate
static inline struct extra_options get_extra_options(const weight_method_t weight_method)
{
    struct extra_options extra;
    ENSURE_STRUCT_SIZE(struct extra_options, EXTRA_OPTIONS_HEADER_SIZE);//compile-time check for making sure struct is correct size
    memset(&extra, 0, EXTRA_OPTIONS_HEADER_SIZE);

    extra.weight_method = weight_method;

    weight_struct *w0 = &(extra.weights0);
    weight_struct *w1 = &(extra.weights1);
    w0->num_weights = get_num_weights_by_method(extra.weight_method);
    w1->num_weights = w0->num_weights;

    return extra;
}

static inline void print_cell_timings(struct config_options *options)
{
    fprintf(stderr,"#########################################################################\n");
    fprintf(stderr,"#  Cell_1    Cell_2          N1          N2        Time_ns     ThreadID  \n");
    fprintf(stderr,"#########################################################################\n");
    for(int64_t i=0;i<options->totncells_timings;i++) {
        fprintf(stderr,"%8d %8d %12"PRId64" %12"PRId64" %12"PRId64" %12d\n",
                options->cell_timings[i].first_cellindex,
                options->cell_timings[i].second_cellindex,
                options->cell_timings[i].N1,
                options->cell_timings[i].N2,
                options->cell_timings[i].time_in_ns,
                options->cell_timings[i].tid);
    }

}

static inline void free_cell_timings(struct config_options *options)
{
    if(options->totncells_timings > 0 && options->cell_timings != NULL) {
        free(options->cell_timings);
    }
    options->totncells_timings = 0;

    return;
}

static inline void allocate_cell_timer(struct config_options *options, const int64_t num_cell_pairs)
{
    if(options->totncells_timings >= num_cell_pairs) return;

    free_cell_timings(options);
    options->cell_timings = calloc(num_cell_pairs, sizeof(*(options->cell_timings)));
    if(options->cell_timings == NULL) {
        fprintf(stderr,"Warning: In %s> Could not allocate memory to store the API timings per cell. \n",
                __FUNCTION__);
    } else {
        options->totncells_timings = num_cell_pairs;
    }

    return;
}

static inline void assign_cell_timer(struct api_cell_timings *cell_timings, const int64_t num_cell_pairs, struct config_options *options)
{
    /* Does the existing thread timings pointer have enough memory allocated ?*/
    allocate_cell_timer(options, num_cell_pairs);

    /* This looks like a repeated "if" condition but it is not. Covers the case for the calloc failure above */
    if(options->totncells_timings >= num_cell_pairs) {
        memmove(options->cell_timings, cell_timings, sizeof(struct api_cell_timings) * num_cell_pairs);
    }
}


#include "macros.h"


#ifdef __cplusplus
}
#endif
