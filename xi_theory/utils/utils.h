#ifndef __UTILS_H
#define __UTILS_H
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>//defines int64_t datatype -> *exactly* 8 bytes int 
#include<inttypes.h>//defines PRId64 for printing int64_t
#include<math.h>
#include<string.h>
#include<limits.h>
#include<assert.h>
#include<time.h>
#include<sys/time.h>
#include<stdarg.h>
/* #include<malloc.h> */

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
//short float_almost_equal(float A, float B, int maxUlps);

//memory routines
extern void* my_realloc(void *x,size_t size,int64_t N,const char *varname);
extern void* my_realloc_in_function(void **x,size_t size,int64_t N,const char *varname);
extern void* my_malloc(size_t size,int64_t N);
extern void* my_align_malloc(size_t size,int64_t N,size_t alignment);
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
extern void  setup_bin_lookup(const char *fname,double *rmin,double *rmax,int *nbin,const size_t nbinlookup,double **rupp,int *binlookup);
extern void  setup_squared_bin_lookup(const char *fname,double *rmin,double *rmax,int *nbin,const size_t nbinlookup,double **rupp,int *binlookup);
extern void  setup_bin_lookup_float(const char *fname,float *rmin,float *rmax,int *nbin,const size_t nbinlookup,float **rupp,int *binlookup);
//end function declarations


#endif
