/* File: tests_common.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <inttypes.h>

#ifndef MAXLEN
#define MAXLEN 500
#endif

#include "defs.h"
#include "utils.h"

#ifdef INTEGRATION_TESTS

#pragma message("Running (SLOW) integration tests")

/* Define the instruction sets that are supported by the compiler */
const isa valid_instruction_sets[] = {FALLBACK
#ifdef __SSE4_2__
                                      ,SSE42
#endif
#ifdef __AVX__
                                      ,AVX
#endif
#ifdef __AVX2__
                                      ,AVX2
#endif
#ifdef __AVX512F__
                                      ,AVX512F
#endif
};

/* Strings corresponding to the instruction sets in the array `valid_instruction_sets` */
const char isa_name[][20] = {"FALLBACK"
#ifdef __SSE4_2__
                             ,"SSE42"
#endif
#ifdef __AVX__
                             , "AVX"
#endif
#ifdef __AVX2__
                             , "AVX2"
#endif
#ifdef __AVX512F__
                             , "AVX512F"
#endif
};

/* This is a fun C tid-bit. The sizeof(valid_instruction_sets) refers to the total bytes
   required to store the array. As in the typeof valid_instruction_sets is int[5] when
   all 5 instructions sets are supported */
const int num_instructions = sizeof(valid_instruction_sets)/sizeof(valid_instruction_sets[0]);

/* The max. value of bin refine factor to probe. Each of bin refinements factors is set from [1, max_binref]
 (inclusive) */
const int min_bin_ref = 1, max_bin_ref = 3;

/* Macro to setup the loop over instruction sets, various bin factors and then run
 the tests */
#define BEGIN_INTEGRATION_TEST_SECTION                                  \
    do {                                                                \
           int dotest = 1;                                              \
           const isa old_isa = options.instruction_set;                 \
           const int old_min_sep_opt = options.enable_min_sep_opt;      \
           const int old_copy_parts  = options.copy_particles;      \
           struct timespec t0, t1;                                      \
           for(int iset=0;iset<num_instructions;iset++) {               \
               int fastest_bin_ref[] = {1, 1, 1};                       \
               double fastest_time = 1e30;                              \
               options.instruction_set = valid_instruction_sets[iset];  \
               for(int bfx=min_bin_ref;bfx<=max_bin_ref;bfx++) {    \
                   for(int bfy=min_bin_ref;bfy<=max_bin_ref;bfy++) { \
                       for(int bfz=min_bin_ref;bfz<=max_bin_ref;bfz++) { \
                           for(int copy_parts=1;copy_parts >= 0; copy_parts--) { \
                               options.copy_particles = copy_parts;         \
                               for(int enable_min_sep_opt=0;enable_min_sep_opt<=1;enable_min_sep_opt++) { \
                                   if(dotest == 1) {                    \
                                       const int bf[] = {bfx, bfy, bfz}; \
                                       set_custom_bin_refine_factors(&options, bf); \
                                       options.enable_min_sep_opt = enable_min_sep_opt; \
                                       fprintf(stderr,"Bin refs = (%d, %d, %d), copy_particles = %5s, and min. sep. opt %8s ...", \
                                               options.bin_refine_factors[0], \
                                               options.bin_refine_factors[1], \
                                               options.bin_refine_factors[2], \
                                               copy_parts == 1 ? "TRUE":"FALSE", \
                                               enable_min_sep_opt == 0 ? "DISABLED":"ENABLED"); \
                                       current_utc_time(&t0);


/* Clean up the integration tests (close the loops and check for error) */
#define END_INTEGRATION_TEST_SECTION(code_to_free_results_memory)                                  \
                                       current_utc_time(&t1);                                              \
                                       double time_to_run = REALTIME_ELAPSED_NS(t0, t1); \
                                       if(time_to_run < fastest_time) { \
                                           fastest_time = time_to_run;  \
                                           memcpy(&fastest_bin_ref, &bf, sizeof(bf)); \
                                       }                                \
                                       if(ret != EXIT_SUCCESS) {        \
                                           fprintf(stderr, ANSI_COLOR_RED "FAILED"); \
                                           dotest = 0;                  \
                                       } else {                         \
                                           fprintf(stderr,ANSI_COLOR_GREEN "PASSED"); \
                                           code_to_free_results_memory; /* whatever code is required to free the memory in results struct */ \
                                       }                                \
                                       fprintf(stderr, ANSI_COLOR_RESET " (isa = %s). Time = %6.2lf seconds \n", \
                                               isa_name[iset],time_to_run * 1e-9); \
                                   } /* close the enable_min_sep_opt condition*/ \
                               }/* copy particle positions*/                                           \
                           }/* close the dotest if condition*/ \
                       }/*bin ref z*/                   \
                   }/*bin ref y*/                   \
               }/*bin ref x*/                                           \
               if(ret == EXIT_SUCCESS) {        \
                   fprintf(stderr, ANSI_COLOR_MAGENTA "Fastest time = %8.2lf seconds with bin-ref = {%d, %d, %d}" ANSI_COLOR_RESET "\n", \
                           fastest_time*1e-9,                           \
                           fastest_bin_ref[0],                          \
                           fastest_bin_ref[1],                          \
                           fastest_bin_ref[2]);                         \
               }                                                        \
           } /*instruction set */                                       \
           reset_bin_refine_factors(&options);                          \
           options.instruction_set = old_isa;                           \
           options.enable_min_sep_opt = old_min_sep_opt;                \
           options.copy_particles = old_copy_parts;            \
    } while(0)

    //wtheta has 3 implementations (brute-force, link-in-dec and link-in-dec + link-in-ra)
    //For developer testing, multiple bin refine factors are tested as well as the
    //all three of the linking logic.

    // (dec_link, ra_link) == (0, 0) -> brute-force
    // (dec_link, ra_link) == (1, 0) -> dec-linking only
    // (dec_link, ra_link) == (1, 1) -> dec + ra linking

    /* The order of the for loop breaks the convention "RA before DEC"
       -- This is because the binning in RA can only be done if the binning
       in DEC is enabled. Therefore, it makes more sense to loop in RA *only*
       after the DEC binning is decided.
     */

#define BEGIN_DDTHETA_INTEGRATION_TEST_SECTION                               \
    do {                                                                     \
        struct timespec t0, t1;                                              \
        const isa old_isa = options.instruction_set;                         \
        const int old_min_sep_opt = options.enable_min_sep_opt;              \
        const int old_copy_parts  = options.copy_particles;             \
        int dotest = 1;                                                 \
        for(int iset=0;iset<num_instructions;iset++) {                       \
            options.instruction_set = valid_instruction_sets[iset];          \
            for(int dec_link=0;dec_link <= 1;dec_link++) {                   \
                for(int ra_link=0;ra_link <= dec_link; ra_link++) {          \
                    int fastest_bin_ref[] = {1, 1, 1};                       \
                    int fastest_isa = 0;                                     \
                    double fastest_time = 1e30;                              \
                    for(int dec_bin_ref=min_bin_ref;dec_bin_ref<=max_bin_ref;dec_bin_ref++) { \
                        for(int ra_bin_ref=min_bin_ref;ra_bin_ref<=max_bin_ref;ra_bin_ref++) { \
                            for(int copy_parts=1;copy_parts >= 0; copy_parts--) { \
                                options.copy_particles = copy_parts;         \
                                for(int enable_min_sep_opt=0;enable_min_sep_opt<=1;enable_min_sep_opt++) { \
                                    if(dotest == 1) {                        \
                                        if(dec_link == 0 && ra_link == 0 && (dec_bin_ref != min_bin_ref || \
                                                                             ra_bin_ref != min_bin_ref || \
                                                                             copy_parts != 0 || \
                                                                             enable_min_sep_opt != 0)) continue; \
                                        if(dec_link == 1 && ra_link == 0 && ra_bin_ref != min_bin_ref) continue; \
                                        if(dec_link == 1 && ra_link == 0 && ra_bin_ref != min_bin_ref) continue; \
                                        const int bf[] = {ra_bin_ref, dec_bin_ref, -1}; \
                                        set_custom_bin_refine_factors(&options, bf); \
                                        options.link_in_dec=dec_link;             \
                                        options.link_in_ra=ra_link;               \
                                        options.enable_min_sep_opt = enable_min_sep_opt; \
                                        if(dec_link == 1) {                 \
                                            fprintf(stderr,"With (dec, ra)-linking = (%1d, %1d), (dec, ra) bin-ref = (%d, %d), "\
                                                    "copy_particles = %5s, min. sep. opt. %8s ...", \
                                                    dec_link, ra_link,  \
                                                    options.bin_refine_factors[1], \
                                                    options.bin_refine_factors[0], \
                                                    copy_parts == 1 ? "TRUE":"FALSE", \
                                                    enable_min_sep_opt == 0 ? "DISABLED":"ENABLED"); \
                                        }                                   \
                                        current_utc_time(&t0);      \


/* Clean up the integration tests (close the loops and check for error) */
#define END_DDTHETA_INTEGRATION_TEST_SECTION(code_to_free_results_memory)                          \
                                        current_utc_time(&t1);                    \
                                        double time_to_run = REALTIME_ELAPSED_NS(t0, t1); \
                                        if(time_to_run < fastest_time) {          \
                                            fastest_time = time_to_run;           \
                                            fastest_isa = iset;                   \
                                            memcpy(&fastest_bin_ref, &bf, sizeof(bf)); \
                                        }                                         \
                                        if(ret != EXIT_SUCCESS) {                 \
                                            fprintf(stderr, ANSI_COLOR_RED "FAILED"); \
                                            dotest = 0;                           \
                                        } else {                                  \
                                            fprintf(stderr,ANSI_COLOR_GREEN "PASSED"); \
                                            code_to_free_results_memory; /* whatever code is required to free the memory in results struct */ \
                                        }                                         \
                                        fprintf(stderr, ANSI_COLOR_RESET " (isa = %s). Time = %6.2lf seconds \n", \
                                                isa_name[iset],time_to_run * 1e-9); \
                                    } /* closes dotest */                         \
                                } /* cloeses copy_particles */   \
                            } /* closes min-sep-opt*/                         \
                        } /* dec_bin_ref*/                                    \
                    } /* ra_bin_ref*/                                         \
                    if(ret == EXIT_SUCCESS) {                                 \
                        fprintf(stderr, ANSI_COLOR_MAGENTA "Fastest time = %8.2lf seconds with (dec, ra) bin-ref = {%d, %d} and instruction_set = %s"  ANSI_COLOR_RESET "\n", \
                            fastest_time*1e-9,                                \
                            fastest_bin_ref[1],                               \
                            fastest_bin_ref[0],                               \
                            isa_name[fastest_isa]);                           \
                    } /* closes ret==EXIT_SUCCESS */                          \
                } /* ra_link loop */                                          \
            } /*dec_link loop */                                              \
        } /* isa loop */                                                      \
        reset_bin_refine_factors(&options);                                   \
        options.instruction_set = old_isa;                                    \
        options.copy_particles = old_copy_parts;                        \
        options.enable_min_sep_opt = old_min_sep_opt;                   \
    } while(0)

#else
/* Running regular tests -> no need for exhaustive testing */
#define BEGIN_INTEGRATION_TEST_SECTION  do {
#define END_INTEGRATION_TEST_SECTION(code_to_free_results_memory)    } while(0)

#define BEGIN_DDTHETA_INTEGRATION_TEST_SECTION do {
#define END_DDTHETA_INTEGRATION_TEST_SECTION(code_to_free_results_memory)   } while(0)

#endif/*INTEGRATION_TESTS*/


#ifdef _OPENMP
const int nthreads=4;
#else
const int nthreads=1;
#endif

const double maxdiff = 1e-9;
const double maxreldiff = 1e-6;

char binfile[]="../tests/bins";
char angular_binfile[]="../tests/angular_bins";
double pimax=40.0;
double theory_mu_max=0.5;
double mocks_mu_max=1.0;
int nmu_bins=10;
double boxsize=420.0;
