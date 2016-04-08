/* File: utils.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>//defines int64_t datatype -> *exactly* 8 bytes int
#include<inttypes.h>//defines PRId64 for printing int64_t + includes stdint.h
#include<math.h>
#include<string.h>
#include<limits.h>
#include<assert.h>
#include<time.h>
#include<sys/time.h>
#include<stdarg.h>

#ifdef __cplusplus
 extern "C" {
#endif

     //Just to output some colors
    
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                               \
     do { if (!(EXP)) {                                                 \
             fprintf(stderr,"Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
             fprintf(stderr,__VA_ARGS__);                               \
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please email Manodeep Sinha <manodeep@gmail.com>"ANSI_COLOR_RESET"\n"); \
             exit(EXIT_FAILURE);                                        \
         }                                                              \
     } while (0)
#endif

     
     //routines for file i/o
     extern FILE * my_fopen(const char *fname,const char *mode);
     extern FILE * my_fopen_carefully(const char *fname,void (*header)(FILE *));
     extern size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
     extern size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
     extern int my_fseek(FILE *stream, long offset, int whence);

     //general utilities
     extern char *int2bin(int a, char *buffer, int buf_size) ;
     extern int my_snprintf(char *buffer,int len,const char *format, ...) __attribute__((format(printf,3,4)));
     extern void print_time(struct timeval t0,struct timeval t1,const char *s);
     extern int64_t getnumlines(const char *fname,const char comment);

     //memory routines
     extern void* my_realloc(void *x,size_t size,int64_t N,const char *varname);
     extern void* my_realloc_in_function(void **x,size_t size,int64_t N,const char *varname);
     extern void* my_malloc(size_t size,int64_t N);
     extern  void* my_calloc(size_t size,int64_t N);
     extern void my_free(void ** x);
     extern void **matrix_malloc(size_t size,int64_t nx,int64_t ny);
     extern void **matrix_calloc(size_t size,int64_t nx,int64_t ny);
     extern void matrix_free(void **m,int64_t ny);

     void *** volume_malloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe);
     void *** volume_calloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe);
     void volume_free(void ***v,int64_t nrow,int64_t ncol);

     extern void run_system_call(const char *execstring);

     extern void  setup_bins(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp);

     extern int test_all_files_present(const int nfiles, ...);
     //end function declarations

#ifdef __cplusplus
 }
#endif
