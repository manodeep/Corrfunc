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

#warning "Running (SLOW) integration tests"

/* Define the instruction sets that are supported by the compiler */
const isa valid_instruction_sets[] = {FALLBACK
#ifdef __SSE4_2__
                                      ,SSE42
#endif                                      
#ifdef __AVX__                                      
                                      ,AVX
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
};

/* This is a fun C tid-bit. The sizeof(valid_instruction_sets) refers to the total bytes
   required to store the array. As in the typeof valid_instruction_sets is int[3] when
   all 3 instructions sets are supported */
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
           int fastest_bin_ref[] = {1, 1, 1};                           \
           double fastest_time = 1e30;                                  \
           struct timespec t0, t1;                                      \
           for(int iset=0;iset<num_instructions;iset++) {                     \
               options.instruction_set = valid_instruction_sets[iset];    \
               for(int bfx=min_bin_ref;bfx<=max_bin_ref;bfx++) {         \
                   for(int bfy=min_bin_ref;bfy<=max_bin_ref;bfy++) {     \
                       for(int bfz=min_bin_ref;bfz<=max_bin_ref;bfz++) { \
                           if(dotest == 1) {                            \
                               const int bf[] = {bfx, bfy, bfz};        \
                               set_custom_bin_refine_factors(&options, bf); \
                               fprintf(stderr,"Running with bin refs = (%d, %d, %d) and instruction set = %s...", \
                                       options.bin_refine_factors[0],   \
                                       options.bin_refine_factors[1],   \
                                       options.bin_refine_factors[2],   \
                                       isa_name[iset]);                   \
                               current_utc_time(&t0);


/* Clean up the integration tests (close the loops and check for error) */
#define END_INTEGRATION_TEST_SECTION                                    \
                               current_utc_time(&t1);                   \
                               double time_to_run = REALTIME_ELAPSED_NS(t0, t1);\
                               if(time_to_run < fastest_time) {         \
                                   fastest_time = time_to_run;          \
                                   memcpy(&fastest_bin_ref, &bf, sizeof(bf));\
                               }                                        \
                               if(ret != EXIT_SUCCESS) {                \
                                   fprintf(stderr, ANSI_COLOR_RED "FAILED"); \
                                   dotest = 0;                          \
                               } else {                                 \
                                   fprintf(stderr,ANSI_COLOR_GREEN "PASSED"); \
                               }                                        \
                               fprintf(stderr, ANSI_COLOR_RESET ". Time taken = %8.2lf seconds \n", time_to_run * 1e-9); \
                           }/* close the dotest if condition*/          \
                       }/*bin ref z*/                                   \
                   }/*bin ref y*/                                       \
               }/*bin ref x*/                                           \
               if(ret == EXIT_SUCCESS) {                                \
                   fprintf(stderr, ANSI_COLOR_MAGENTA "Fastest time = %8.2lf seconds with bin-ref = {%d, %d, %d}" ANSI_COLOR_RESET "\n", \
                           fastest_time*1e-9,                           \
                           fastest_bin_ref[0],                          \
                           fastest_bin_ref[1],                          \
                           fastest_bin_ref[2]);                         \
               }                                                        \
           } /*instruction set */                                       \
           reset_bin_refine_factors(&options);                          \
           options.instruction_set = old_isa;                           \
    } while(0)
#else
/* Running regular tests -> no need for exhaustive testing */
#define BEGIN_INTEGRATION_TEST_SECTION  do {                               
#define END_INTEGRATION_TEST_SECTION    } while(0)

#endif


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

