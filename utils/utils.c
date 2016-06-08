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

#include "utils.h"
/* #define __USE_XOPEN2K */
/* #define _XOPEN_SOURCE_EXTENDED */
/* #define _GNU_SOURCE  */


void  setup_bins(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp)
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
}


void run_system_call(const char *execstring)
{
    int status;
    status=system(execstring);
    if(status != EXIT_SUCCESS) {
        fprintf(stderr,"ERROR: executing system command: \n`%s'\n...exiting\n",execstring);
        exit(EXIT_FAILURE);
    }

}



FILE * my_fopen(const char *fname,const char *mode)
{
    FILE *fp=NULL;
    fp = fopen(fname,mode);
    if(fp == NULL)
        {
            fprintf(stderr,"Could not open file `%s'\n",fname);
            exit(EXIT_FAILURE);
        }
    return fp;
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
    FILE *fp = NULL;
    fp = fopen(fname,"r");//note I am using fopen and not my_fopen.

    if(fp == NULL)
        {
            /*file does not exist -> open with "w" */
            fp = my_fopen(fname,"w");//using my_fopen here.
            (*header)(fp);/* print the header */
        }
    else
        {
            fclose(fp);
            fp = my_fopen(fname,"a+");//open with append mode
        }

    return fp;
}


size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t nwritten;
    nwritten = fwrite(ptr, size, nmemb, stream);
    if(nwritten != nmemb)
        {
            fprintf(stderr,"I/O error (fwrite) has occured.\n");
            fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n",nmemb,nwritten);
            exit(EXIT_FAILURE);
        }
    return nwritten;
}

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t nread;
    nread = fread(ptr, size, nmemb, stream);
    if(nread != nmemb) {
        fprintf(stderr,"I/O error (fread) has occured.\n");
        fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n",nmemb,nread);
        exit(EXIT_FAILURE);
    }
    return nread;
}

int my_fseek(FILE *stream, long offset, int whence)
{
    int err=fseek(stream,offset,whence);
    if(err != 0) {
        fprintf(stderr,"ERROR: Could not seek `%ld' bytes into the file..exiting\n",offset);
        exit(EXIT_FAILURE);
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
        fprintf(stderr,"Increase `len' in the header file ..exiting\n");
        exit(EXIT_FAILURE);
    }
    return nwritten;
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




void print_time(struct timeval t0,struct timeval t1,const char *s)
{
    double timediff = difftime(t1.tv_sec,t0.tv_sec);
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
        my_free((void **) &x);
        exit(EXIT_FAILURE);
    }

    return tmp;

}

void* my_malloc(size_t size,int64_t N)
{
    void *x = NULL;
    x = malloc(N*size);
    if (x==NULL){
        fprintf(stderr,"malloc for %"PRId64" elements with %zu bytes failed..aborting\n",N,size);
        exit(EXIT_FAILURE);
    }
        
    return x;
}



void* my_calloc(size_t size,int64_t N)
{
    void *x = NULL;
    x = calloc((size_t) N, size);
    if (x==NULL)    {
        fprintf(stderr,"malloc for %"PRId64" elements with %zu size failed..aborting\n",N,size);
        exit(EXIT_FAILURE);
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
    void **m;
    m = (void **) my_malloc(sizeof(void *),nrow);
    for(int i=0;i<nrow;i++)
        m[i] = (void *) my_malloc(size,ncol);

    return m;
}

void **matrix_calloc(size_t size,int64_t nrow,int64_t ncol)
{
    void **m;
    m = (void **) my_calloc(sizeof(void *),nrow);
    for(int i=0;i<nrow;i++)
        m[i] = (void *) my_calloc(size,ncol);

    return m;
}



void matrix_free(void **m,int64_t nrow)
{
    for(int i=0;i<nrow;i++)
        free(m[i]);

    free(m);
}




void *** volume_malloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe)
{
    void ***v;
    v = (void ***) my_malloc(sizeof(void **),nrow);
    for(int i=0;i<nrow;i++) {
        v[i] = (void *) my_malloc(sizeof(void *),ncol);
        for(int j=0;j<ncol;j++) {
            v[i][j] = my_malloc(size,nframe);
        }
    }

    return v;
}

void *** volume_calloc(size_t size,int64_t nrow,int64_t ncol,int64_t nframe)
{
    void ***v;
    v = (void ***) my_malloc(sizeof(void **),nrow);
    for(int i=0;i<nrow;i++) {
        v[i] = (void *) my_malloc(sizeof(void *),ncol);
        for(int j=0;j<ncol;j++) {
            v[i][j] = my_calloc(size,nframe);
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
