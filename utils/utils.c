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

#include "macros.h"
#include "utils.h"

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

int setup_bins_float(const char *fname,float *rmin,float *rmax,int *nbin,float **rupp)
{
    //set up the bins according to the binned data file
    //the form of the data file should be <rlow  rhigh ....>
    const int MAXBUFSIZE=1000;
    char buf[MAXBUFSIZE];
    FILE *fp=NULL;
    float low,hi;
    const char comment='#';
    const int nitems=2;
    int nread=0;
    *nbin = ((int) getnumlines(fname,comment))+1;
    *rupp = my_calloc(sizeof(float),*nbin+1);

    fp = my_fopen(fname,"r");
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
  char units[4][10]  = {"days", "hrs" , "mins", "secs"};
  int which = 0;

  double timeleft = timediff;
  double time_to_print;
  
  if(timediff < ratios[2]) {
    my_snprintf(time_string, MAXLINESIZE,"%6.3lf secs",1e-6*(t1.tv_usec-t0.tv_usec) + timediff);
  }  else {
    size_t curr_index = 0;
    while (which < 4) {
      time_to_print = floor(timeleft/ratios[which]);
      if (time_to_print > 1) {
        timeleft -= (time_to_print*ratios[which]);
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
    char units[4][10]  = {"days", "hrs" , "mins", "secs"};
    int which = 0;

    double timeleft = timediff;
    double time_to_print;
    fprintf(stderr,"Time taken to execute '%s'  = ",s);

    if(timediff < ratios[2]) {
        fprintf(stderr,"%6.3lf secs",1e-6*(t1.tv_usec-t0.tv_usec) + timediff);
    }  else {
        while (which < 4) {
            time_to_print = floor(timeleft/ratios[which]);
            if (time_to_print > 1) {
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
    void *x = NULL;
    x = malloc(N*size);
    if (x==NULL){
        fprintf(stderr,"malloc for %"PRId64" elements with %zu bytes failed...\n",N,size);
        perror(NULL);
    }
        
    return x;
}



void* my_calloc(size_t size,int64_t N)
{
    void *x = NULL;
    x = calloc((size_t) N, size);
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
    FILE *fp= NULL;
    const int MAXLINESIZE = 10000;
    int64_t nlines=0;
    char str_line[MAXLINESIZE];

    fp = my_fopen(fname,"rt");
    if(fp == NULL) {
        return -1;
    }

    while(1){
        if(fgets(str_line, MAXLINESIZE,fp)!=NULL) {
            //WARNING: this does not remove white-space. You might
            //want to implement that (was never an issue for me)
            if(str_line[0] !=comment)
                nlines++;
        } else
            break;
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


/* #undef __USE_XOPEN2K */
