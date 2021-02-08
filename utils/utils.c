/* File: utils.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/*
  A collection of C wrappers I use. Should be
  very obvious. The ones that are not obvious
  have comments before the function itself.

  Bugs:
  Please email me manodeep at gmail dot com

  Ver 1.0: Manodeep Sinha, 2nd April, 2012
  Ver 1.1: Manodeep Sinha, 14th June, 2012 - replaced
  check_string_copy with a "real" wrapper to
  snprintf.
  Ver 1.2: Manodeep Sinha, Jan 8, 2012 - replaced
  print_time with timeval and gettimeofday
*/

#include<inttypes.h>//defines PRId64 for printing int64_t + includes stdint.h
#include<math.h>
#include<string.h>
#include<limits.h>
#include<stdarg.h>
#include<ctype.h>

#include "macros.h"
#include "utils.h"

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
#include <mach/mach_time.h> /* mach_absolute_time -> really fast */
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

void get_max_float(const int64_t ND1, const float *cz1, float *czmax)
{
    float max=*czmax;
    for(int64_t i=0;i<ND1;i++) {
        if(cz1[i] > max) max = cz1[i];
    }
    *czmax = max;
}

void get_max_double(const int64_t ND1, const double *cz1, double *czmax)
{
    double max=*czmax;
    for(int64_t i=0;i<ND1;i++) {
        if(cz1[i] > max) max = cz1[i];
    }
    *czmax = max;
}


int setup_bins(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp)
{
    //set up the bins according to the binned data file
    //the form of the data file should be <rlow  rhigh ....>
    const int MAXBUFSIZE=1000;
    char buf[MAXBUFSIZE];
    FILE *fp=NULL;
    double low,hi;
    const char comment='#';
    const int nitems=2;
    int nread=0;
    *nbin = ((int) getnumlines(fname,comment))+1;
    *rupp = my_calloc(sizeof(double),*nbin+1);

    fp = my_fopen(fname,"r");
    if(fp == NULL) {
        free(*rupp);
        return EXIT_FAILURE;
    }
    int index=1;
    while(1) {
        if(fgets(buf,MAXBUFSIZE,fp)!=NULL) {
            nread=sscanf(buf,"%lf %lf",&low,&hi);
            if(nread==nitems) {

                if(index==1) {
                    *rmin=low;
                    (*rupp)[0]=low;
                }

                (*rupp)[index] = hi;
                index++;
            }
        } else {
            break;
        }
    }
    *rmax = (*rupp)[index-1];
    fclose(fp);

    (*rupp)[*nbin]=*rmax ;
    (*rupp)[*nbin-1]=*rmax ;

    return EXIT_SUCCESS;
}



int setup_bins_double(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp)
{
    //set up the bins according to the binned data file
    //the form of the data file should be <rlow  rhigh ....>
    const int MAXBUFSIZE=1000;
    char buf[MAXBUFSIZE];
    double low,hi;
    const char comment='#';
    const int nitems=2;
    int nread=0;
    *nbin = ((int) getnumlines(fname,comment))+1;
    *rupp = my_calloc(sizeof(double),*nbin+1);

    FILE *fp = my_fopen(fname,"r");
    if(fp == NULL) {
        free(*rupp);
        return EXIT_FAILURE;
    }
    int index=1;
    while(1) {
        if(fgets(buf,MAXBUFSIZE,fp)!=NULL) {
            nread=sscanf(buf,"%lf %lf",&low,&hi);
            if(nread==nitems) {
                if(index==1) {
                    *rmin=low;
                    (*rupp)[0]=low;
                }

                (*rupp)[index] = hi;
                index++;
            }
        } else {
            break;
        }
    }
    *rmax = (*rupp)[index-1];
    fclose(fp);

    (*rupp)[*nbin]=*rmax ;
    (*rupp)[*nbin-1]=*rmax ;

    return EXIT_SUCCESS;
}

int setup_bins_float(const char *fname,float *rmin,float *rmax,int *nbin,float **rupp)
{
    //set up the bins according to the binned data file
    //the form of the data file should be <rlow  rhigh ....>
    const int MAXBUFSIZE=1000;
    char buf[MAXBUFSIZE];
    float low,hi;
    const char comment='#';
    const int nitems=2;
    int nread=0;
    *nbin = ((int) getnumlines(fname,comment))+1;
    *rupp = my_calloc(sizeof(float),*nbin+1);

    FILE *fp = my_fopen(fname,"r");
    if(fp == NULL) {
        free(*rupp);
        return EXIT_FAILURE;
    }
    int index=1;
    while(1) {
        if(fgets(buf,MAXBUFSIZE,fp)!=NULL) {
            nread=sscanf(buf,"%f %f",&low,&hi);
            if(nread==nitems) {

                if(index==1) {
                    *rmin=low;
                    (*rupp)[0]=low;
                }

                (*rupp)[index] = hi;
                index++;
            }
        } else {
            break;
        }
    }
    *rmax = (*rupp)[index-1];
    fclose(fp);

    (*rupp)[*nbin]=*rmax ;
    (*rupp)[*nbin-1]=*rmax ;

    return EXIT_SUCCESS;
}


int run_system_call(const char *execstring)
{
    int status=system(execstring);
    if(status != EXIT_SUCCESS) {
        fprintf(stderr,"ERROR: executing system command: \n`%s'\n\n",execstring);
        perror(NULL);
    }
    return EXIT_FAILURE;
}



FILE * my_fopen(const char *fname,const char *mode)
{
    FILE *fp = fopen(fname,mode);
    if(fp == NULL){
        fprintf(stderr,"Could not open file `%s'\n",fname);
        perror(NULL);
    }
    return fp;//Could be NULL
}

/*
  The following function opens a file (if it already exists)
  in append mode. If the file doesn't exist, then the function
  creates one, calls the *header() function [which presumably
  prints a header to the file] and then returns the file pointer.

  As usual, you need to be careful with the file you are appending
  to -> otherwise you might end up with a ginormous file. Usually,
  I do a system("rm -f filename") before the loop where the file
  might be created/modified and remove the file from previous
  runs.
*/

FILE * my_fopen_carefully(const char *fname,void (*header)(FILE *))
{
    FILE *fp = fopen(fname,"r");//note I am using fopen and not my_fopen.

    if(fp == NULL) {
        /*file does not exist -> open with "w" */
        fp = my_fopen(fname,"w");//using my_fopen here.
        if(fp != NULL) {
            (*header)(fp);/* print the header */
        }
    } else {
        fclose(fp);
        fp = my_fopen(fname,"a+");//open with append mode
    }

    return fp;
}


size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t nwritten;
    nwritten = fwrite(ptr, size, nmemb, stream);
    if(nwritten != nmemb){
        fprintf(stderr,"I/O error (fwrite) has occured.\n");
        fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu \n",nmemb,nwritten);
        perror(NULL);
        return -1;
    }
    return nwritten;
}

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t nread;
    nread = fread(ptr, size, nmemb, stream);
    if(nread != nmemb) {
        fprintf(stderr,"I/O error (fread) has occured.\n");
        fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu\n",nmemb,nread);
        perror(NULL);
        return -1;
    }
    return nread;
}

int my_fseek(FILE *stream, long offset, int whence)
{
    int err=fseek(stream,offset,whence);
    if(err != 0) {
        fprintf(stderr,"ERROR: Could not seek `%ld' bytes into the file..exiting\n",offset);
        perror(NULL);
    }
    return err;
}


// A real wrapper to snprintf that will exit() if the allocated buffer length
// was not sufficient. Usage is the same as snprintf

int my_snprintf(char *buffer,int len,const char *format, ...)
{
    va_list args;
    int nwritten=0;

    va_start(args,format);
    nwritten=vsnprintf(buffer, (size_t) len, format, args );
    va_end(args);
    if (nwritten > len || nwritten < 0) {
        fprintf(stderr,"ERROR: printing to string failed (wrote %d characters while only %d characters were allocated)\n",nwritten,len);
        fprintf(stderr,"Increase `len'=%d in the header file\n",len);
        return -1;
    }
    return nwritten;
}

int is_big_endian(void)
{
    union {
        uint32_t i;
        char c[4];
    } e = { 0x01000000 };

    return e.c[0];
}

void byte_swap(char * const in, const size_t size, char *out)
{
    if(size > 16) {
        fprintf(stderr,"WARNING: In %s> About to byte_swap %zu bytes but no intrinsic C data-type exists with size larger than 16 bytes",
                __FUNCTION__, size);
    }
    //point to the last byte
    char *in_char = (char *) in + (size - 1UL);

    //point to the first byte in output
    char *out_char = out;

    //Start filling in from the front in the output string
    //taking input from the end of the input
    for(size_t i=0;i<size;i++) {
        *out_char = *in_char;
        out_char++;
        in_char--;
    }

}



//Taken from the inter-webs: http://stackoverflow.com/questions/1024389/print-an-int-in-binary-representation-using-c
char * int2bin(int a, char *buffer, int buf_size)
{
    buffer += (buf_size - 1);
    for (int i = 31; i >= 0; i--) {
        *buffer-- = (a & 1) + '0';

        a >>= 1;
    }

    return buffer;
}


/*
Can not remember where I (MS) got this from. Fairly sure
stackoverflow was involved.
Finally taken from http://stackoverflow.com/a/6719178/2237582 */
void current_utc_time(struct timespec *ts)
{

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    static mach_timebase_info_data_t    sTimebaseInfo = {.numer=0, .denom=0};
    uint64_t start = mach_absolute_time();
    if ( sTimebaseInfo.denom == 0 ) {
        mach_timebase_info(&sTimebaseInfo);
    }

    ts->tv_sec = 0;//(start * sTimebaseInfo.numer/sTimebaseInfo.denom) * tv_nsec;
    ts->tv_nsec = start * sTimebaseInfo.numer / sTimebaseInfo.denom;

#if 0
    //Much slower implementation for clock
    //Slows down the code by up to 4x
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
#endif

#else
    clock_gettime(CLOCK_REALTIME, ts);
#endif
}



/*
  I like this particular function. Generic replacement for printing
  (in meaningful units) the actual execution time of a code/code segment.

  The function call should be like this:

  ---------------------------
  struct timeval t_start,t_end;
  gettimeofday(&t_start,NULL);
  do_something();
  gettimeofday(&t_end,NULL);
  print_time(t_start,t_end,"do something");
  ---------------------------

  if the code took 220 mins 30.1 secs
  -> print_time will output `Time taken to execute `do something' = 3 hours 40 mins 30.1 seconds


  (code can be easily extended to include `weeks' as a system of time unit. left to the reader)
*/


char * get_time_string(struct timeval t0,struct timeval t1)
{
  const size_t MAXLINESIZE = 1024;
  char *time_string = my_malloc(sizeof(char), MAXLINESIZE);
  double timediff = t1.tv_sec - t0.tv_sec;
  double ratios[] = {24*3600.0,  3600.0,  60.0,  1};

  if(timediff < ratios[2]) {
      my_snprintf(time_string, MAXLINESIZE,"%6.3lf secs",1e-6*(t1.tv_usec-t0.tv_usec) + timediff);
  }  else {
      double timeleft = timediff;
      size_t curr_index = 0;
      int which = 0;
      while (which < 4) {
          double time_to_print = floor(timeleft/ratios[which]);
          if (time_to_print > 1) {
              timeleft -= (time_to_print*ratios[which]);
              char units[4][10]  = {"days", "hrs" , "mins", "secs"};
              char tmp[MAXLINESIZE];
              my_snprintf(tmp, MAXLINESIZE, "%5d %s",(int)time_to_print,units[which]);
              const size_t len = strlen(tmp);
              const size_t required_len = curr_index + len + 1;
              XRETURN(MAXLINESIZE >= required_len, NULL,
                      "buffer overflow will occur: string has space for %zu bytes while concatenating requires at least %zu bytes\n",
                      MAXLINESIZE, required_len);
              strcpy(time_string + curr_index, tmp);
              curr_index += len;
          }
          which++;
      }
  }

  return time_string;
}

void print_time(struct timeval t0,struct timeval t1,const char *s)
{
    double timediff = t1.tv_sec - t0.tv_sec;
    double ratios[] = {24*3600.0,  3600.0,  60.0,  1};
    fprintf(stderr,"Time taken to execute '%s'  = ",s);
    if(timediff < ratios[2]) {
        fprintf(stderr,"%6.3lf secs",1e-6*(t1.tv_usec-t0.tv_usec) + timediff);
    }  else {
        double timeleft = timediff;
        int which = 0;
        while (which < 4) {
            double time_to_print = floor(timeleft/ratios[which]);
            if (time_to_print > 1) {
                char units[4][10]  = {"days", "hrs" , "mins", "secs"};
                timeleft -= (time_to_print*ratios[which]);
                fprintf(stderr,"%5d %s",(int)time_to_print,units[which]);
            }
            which++;
        }
    }
    fprintf(stderr,"\n");
}


//wrapper for realloc. varname should contain the name of the
//variable being re-allocated -> helps debugging in case of a crash.

void* my_realloc(void *x,size_t size,int64_t N,const char *varname)
{
    void *tmp = realloc(x,N*size);

    if (tmp==NULL) {
        fprintf(stderr,"ERROR: Could not reallocate for %"PRId64" elements with %zu size for variable `%s' ..aborting\n",N,size,varname);
        perror(NULL);
    }

    return tmp;

}

void* my_malloc(size_t size,int64_t N)
{
    void *x = malloc(N*size);
    if (x==NULL){
        fprintf(stderr,"malloc for %"PRId64" elements with %zu bytes failed...\n",N,size);
        perror(NULL);
    }

    return x;
}



void* my_calloc(size_t size,int64_t N)
{
    void *x = calloc((size_t) N, size);
    if (x==NULL)    {
        fprintf(stderr,"malloc for %"PRId64" elements with %zu size failed...\n",N,size);
        perror(NULL);
    }

    return x;
}



//real free. Use only if you are going to check the
//pointer variable afterwards for NULL.
void my_free(void ** x)
{
    /* my_free(void *x) would also free the
       memory but then x would be a local variable
       and the pointer itself in the calling routine
       could not be set to NULL. Hence the pointer
       to pointer business.
    */

    if(*x!=NULL)
        free(*x);//free the memory

    *x=NULL;//sets the pointer in the calling routine to NULL.
}


void **matrix_malloc(size_t size,int64_t nrow,int64_t ncol)
{
    void **m = (void **) my_malloc(sizeof(void *),nrow);
    if(m == NULL) {
        return NULL;
    }

    for(int i=0;i<nrow;i++) {
        m[i] = (void *) my_malloc(size,ncol);
        /* Check if allocation failed */
        if(m[i] == NULL) {
            /* Free up all the memory allocated so far */
            for(int j=i-1;j>=0;j--) {
                free(m[j]);
            }
            free(m);
            return NULL;
        }
    }

    return m;
}

void **matrix_calloc(size_t size,int64_t nrow,int64_t ncol)
{
    void **m = (void **) my_calloc(sizeof(void *),nrow);
    if(m == NULL) {
        return m;
    }
    for(int i=0;i<nrow;i++) {
        m[i] = (void *) my_calloc(size,ncol);
        /* Check if allocation failed */
        if(m[i] == NULL) {
            /* Free up all the memory allocated so far */
            for(int j=i-1;j>=0;j--) {
                free(m[j]);
            }
            free(m);
            return NULL;
        }
    }

    return m;
}


// Resize a matrix.  Returns EXIT_SUCCESS or EXIT_FAILURE.
// Presently only resizing the last dimension is supported, due to
// potential memory leaks when shrinking the first dimension
int matrix_realloc(void **matrix, size_t size, int64_t nrow, int64_t ncol){
    void *tmp;
    for(int i = 0; i < nrow; i++){
        tmp = my_realloc(matrix[i], size, ncol, "matrix_realloc");
        if(tmp == NULL){
            return EXIT_FAILURE;
        }
        matrix[i] = tmp;
    }

    return EXIT_SUCCESS;
}


void matrix_free(void **m,int64_t nrow)
{
    if(m == NULL)
        return;

    for(int i=0;i<nrow;i++)
        free(m[i]);

    free(m);
}




void *** volume_malloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe)
{
    void ***v = (void ***) my_malloc(sizeof(void **),nrow);
    if( v == NULL) {
        return NULL;
    }
    for(int i=0;i<nrow;i++) {
        v[i] = (void *) my_malloc(sizeof(void *),ncol);
        if(v[i] == NULL) {
            /* Free up all the memory allocated so far */
            for(int jj=i-1;jj>=0;jj--) {
                for(int k=0;k<ncol;k++) {
                    free(v[jj][k]);
                }
            }
            free(v);
            return NULL;
        }

        for(int j=0;j<ncol;j++) {
            v[i][j] = my_malloc(size,nframe);
            if(v[i][j] == NULL) {
                /* Free up all the memory allocated so far */
                /* First free up all columns in this row*/
                for(int k=ncol-1;k>=0;k--) {
                    free(v[i][k]);
                }
                /* Now free all previous rows with all ncols */
                for(int jj=i-1;jj>=0;jj--) {
                    for(int k=0;k<ncol;k++) {
                        free(v[jj][k]);
                    }
                }
                free(v);
                return NULL;
            }
        }
    }

    return v;
}

void *** volume_calloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe)
{
    void ***v = (void ***) my_malloc(sizeof(void **),nrow);
    if(v == NULL) {
        return NULL;
    }

    for(int i=0;i<nrow;i++) {
        v[i] = (void *) my_malloc(sizeof(void *),ncol);
        if(v[i] == NULL) {
            /* Free up all the memory allocated so far */
            for(int jj=i-1;jj>=0;jj--) {
                for(int k=0;k<ncol;k++) {
                    free(v[jj][k]);
                }
            }
            free(v);
            return NULL;
        }

        for(int j=0;j<ncol;j++) {
            v[i][j] = my_calloc(size,nframe);
            if(v[i][j] == NULL) {
                /* Free up all the memory allocated so far */
                /* First free up all columns in this row*/
                for(int k=ncol-1;k>=0;k--) {
                    free(v[i][k]);
                }
                /* Now free all previous rows with all ncols */
                for(int jj=i-1;jj>=0;jj--) {
                    for(int k=0;k<ncol;k++) {
                        free(v[j][k]);
                    }
                }
                free(v);
                return NULL;
            }
        }
    }

    return v;

}



void volume_free(void ***v,int64_t nrow,int64_t ncol)
{
    for(int i=0;i<nrow;i++) {
        for(int j=0;j<ncol;j++) {
            free(v[i][j]);
        }

        free(v[i]);
    }

    free(v);
}



int64_t getnumlines(const char *fname,const char comment)
{
    const int MAXLINESIZE = 10000;
    int64_t nlines=0;
    char str_line[MAXLINESIZE];

    FILE *fp = my_fopen(fname,"rt");
    if(fp == NULL) {
        return -1;
    }

    while(1){
        if(fgets(str_line, MAXLINESIZE,fp)!=NULL) {
            /*
              fgets always terminates the string with a '\0'
              on a successful read
             */
            char *c = &str_line[0];
            while(*c != '\0' && isspace(*c)) {
                c++;
            }
            if(*c != '\0' && *c !=comment) {
                 nlines++;
            }
        } else {
            break;
        }
    }
    fclose(fp);
    return nlines;
}


int test_all_files_present(const int nfiles, ...)
{
    /* sets i'th bit for i'th missing file
       return value is 0 *iff* all files are present
       and readable.
    */

    int absent=0;
    va_list filenames;
    va_start(filenames, nfiles);
    XASSERT(nfiles <= 31, "Can only test for 31 files simultaneously. nfiles = %d \n",nfiles);
    for(int i=0;i<nfiles;i++) {
        const char *f = va_arg(filenames, const char *);
        FILE *fp = fopen(f,"r");
        if(fp == NULL) {
            absent |= 1;
        } else {
            fclose(fp);
        }
        absent <<= 1;
    }
    va_end(filenames);

    return absent;
}


/* int float_almost_equal(const float A, const float B, int maxUlps) */
/* { */
/*     /\* MS -- taken from */
/*        http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm */
/*     *\/ */

/*     const int upper_limit_maxulps = 4 * 1024 * 1024; */
/*     /\* Make sure maxUlps is non-negative and small enough that the */
/*        default NAN won't compare as equal to anything.*\/ */
/*     if(maxUlps <= 0 || maxUlps >= upper_limit_maxulps){ */
/*         fprintf(stderr,"Error: Comparison between floats should have smaller number of max. units in last place. Found maxUlps = %d (max allowed = %d)\n", */
/*                 maxUlps, upper_limit_maxulps); */
/*         return EXIT_FAILURE; */
/*     } */
/*     int aInt = *(int*)&A; */

/*     /\* Make aInt lexicographically ordered as a twos-complement int*\/ */
/*     if (aInt < 0) */
/*         aInt = 0x80000000 - aInt; */

/*     /\* Make bInt lexicographically ordered as a twos-complement int*\/ */

/*     int bInt = *(int*)&B; */
/*     if (bInt < 0) */
/*         bInt = 0x80000000 - bInt; */

/*     int intDiff = abs(aInt - bInt); */
/*     if (intDiff <= maxUlps) */
/*         return 1; */

/*     return 0; */
/* } */


/* Directly taken from https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ */
int AlmostEqualRelativeAndAbs_float(float A, float B,
                                    const float maxDiff,
                                    const float maxRelDiff)
{
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    float diff = fabsf(A - B);
    if (diff <= maxDiff)
        return EXIT_SUCCESS;

    A = fabsf(A);
    B = fabsf(B);
    float largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}

int AlmostEqualRelativeAndAbs_double(double A, double B,
                                     const double maxDiff,
                                     const double maxRelDiff)
{
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return EXIT_SUCCESS;

    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
        return EXIT_SUCCESS;

    /* fprintf(stderr,"diff = %e largest * maxRelDiff = %e\n", diff, largest * maxRelDiff); */
    return EXIT_FAILURE;
}

/* #undef __USE_XOPEN2K */

/* A parallel cumulative sum
   Output convention is: cumsum[0] = 0; cumsum[N-1] = sum(a[0:N-1]);
   The algorithm is:
   - Divide the array into `nthreads` chunks
   - cumsum within each chunk
   - compute the "offset" for each chunk by summing the cumsum at the tail of all previous chunks
   - apply the offset
*/
void parallel_cumsum(const int64_t *a, const int64_t N, int64_t *cumsum){
    if (N <= 0){
        return;  // nothing to do
    }
    
    #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    #else
    int nthreads = 1;
    #endif
    
    // We will heuristically limit the number of threads
    // if there isn't enough work for multithreading to be efficient.
    // This is also important for the correctness of the algorithm below,
    // since it enforces nthreads <= N
    int64_t min_N_per_thread = 10000;
    if(N/min_N_per_thread < nthreads){
        nthreads = N/min_N_per_thread;
    }
    if(nthreads < 1){
        nthreads = 1;
    }
    
    #ifdef _OPENMP
    #pragma omp parallel num_threads(nthreads)
    #endif
    {
        #ifdef _OPENMP
        int tid = omp_get_thread_num();
        #else
        int tid = 0;
        #endif
        
        int64_t cstart = N*tid/nthreads;
        int64_t cend = N*(tid+1)/nthreads;
        cumsum[cstart] = cstart > 0 ? a[cstart-1] : 0;
        for(int64_t c = cstart+1; c < cend; c++){
            cumsum[c] = a[c-1] + cumsum[c-1];
        }
        
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        
        int64_t offset = 0;
        for(int t = 0; t < tid; t++){
            offset += cumsum[N*(t+1)/nthreads-1];
        }
        
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        
        if(offset != 0){
            for(int64_t c = cstart; c < cend; c++){
                cumsum[c] += offset;
            }
        }
    }
}
