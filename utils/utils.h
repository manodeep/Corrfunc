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
#include <time.h>
#include<sys/time.h>
#include<sys/times.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

//routines for file i/o
extern FILE * my_fopen(const char *fname,const char *mode);
extern FILE * my_fopen_carefully(const char *fname,void (*header)(FILE *));
extern size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
extern size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
extern int my_fseek(FILE *stream, long offset, int whence);

//general utilities
extern void get_max_float(const int64_t ND1, const float *cz1, float *czmax);
extern void get_max_double(const int64_t ND1, const double *cz1, double *czmax);
extern char *int2bin(int a, char *buffer, int buf_size) ;
extern int my_snprintf(char *buffer,int len,const char *format, ...) __attribute__((format(printf,3,4)));
extern char * get_time_string(struct timeval t0,struct timeval t1);
extern void print_time(struct timeval t0,struct timeval t1,const char *s);
extern void current_utc_time(struct timespec *ts);
extern int64_t getnumlines(const char *fname,const char comment);
extern int is_big_endian(void);
extern void byte_swap(char * const in, const size_t size, char *out);

//memory routines
extern void* my_realloc(void *x,size_t size,int64_t N,const char *varname);
extern void* my_realloc_in_function(void **x,size_t size,int64_t N,const char *varname);
extern void* my_malloc(size_t size,int64_t N);
extern void* my_calloc(size_t size,int64_t N);
extern void my_free(void ** x);
extern void **matrix_malloc(size_t size,int64_t nx,int64_t ny);
extern void **matrix_calloc(size_t size,int64_t nx,int64_t ny);
extern int matrix_realloc(void **matrix, size_t size, int64_t nrow, int64_t ncol);
extern void matrix_free(void **m,int64_t ny);

void *** volume_malloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe);
void *** volume_calloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe);
void volume_free(void ***v,int64_t nrow,int64_t ncol);

extern int run_system_call(const char *execstring);

extern int setup_bins(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp);
extern int setup_bins_double(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp);
extern int setup_bins_float(const char *fname,float *rmin,float *rmax,int *nbin,float **rupp);

extern int test_all_files_present(const int nfiles, ...);

/* Floating point comparison utilities */
extern int AlmostEqualRelativeAndAbs_float(float A, float B, const float maxDiff, const float maxRelDiff);
extern int AlmostEqualRelativeAndAbs_double(double A, double B, const double maxDiff, const double maxRelDiff);
    
extern void parallel_cumsum(const int64_t *a, const int64_t N, int64_t *cumsum);

//end function declarations

#ifdef __cplusplus
 }
#endif
