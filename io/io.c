/* File: io.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <inttypes.h>
#include <string.h>

#include "io.h"
#include "ftread.h"
#include "utils.h"
#include "macros.h"

#ifndef MEMORY_INCREASE_FAC
#define MEMORY_INCREASE_FAC 1.1
#endif

#ifndef MAXLEN
#define MAXLEN 500
#endif

int64_t read_positions(const char *filename, const char *format, const size_t size, const int num_fields, ...)
{
    XRETURN((sizeof(void *) == sizeof(float *) && sizeof(void *) == sizeof(double *)), -1,
            "Size of void pointer = %zu must be the same as size of float pointer = %zu and sizeof double pointers = %zu\n",
            sizeof(void *), sizeof(float *), sizeof(double *));

    void *columns[num_fields];
    int64_t np = read_columns_into_array(filename, format, size, num_fields, columns);

    va_list ap;
    va_start(ap,num_fields);

    for(int i=0;i<num_fields;i++) {
        void **source = va_arg(ap, void **);
        *source =  columns[i];
    }
    va_end(ap);

    return np;
}


int64_t read_columns_into_array(const char *filename, const char *format, const size_t size, const int num_fields, void **data){
    int64_t np;
    XRETURN(num_fields >= 1, -1, "Number of fields to read-in = %d must be at least 1\n", num_fields);
    XRETURN((size == 4 || size == 8), -1, "Size of fields = %zu must be either 4 or 8\n", size);

    {
        //new scope - just to check if file is gzipped.
        //in that case, use gunzip to unzip the file.
        // *ALL* file open calls in this scope
        // are to fopen and *NOT* my_fopen.

        FILE *fp = fopen(filename,"r");
        if(fp == NULL) {
            /* file does not exist. let's see if the filename.gz file does */
            char buf[MAXLEN] = "";
            my_snprintf(buf, MAXLEN,"%s.gz",filename);
            fp = fopen(buf,"r");
            if (fp == NULL) {
                /* well then, file not found*/
                fprintf(stderr,"ERROR: Could not find file: neither as `%s' nor as `%s'\n",filename,buf);
                return -1;
            } else {
                /* found the gzipped file. Use a system call to uncompress. */
                fclose(fp);

                /*
                  Note, I am using `filename` rather than `buf` both
                  because C standards say that using the same buffer
                  as both source and destination is *undefined behaviour*.

                  Check under "NOTES", towards the end of "man 3 snprintf".
                */
                my_snprintf(buf,MAXLEN,"gunzip %s.gz",filename);
                fprintf(stderr,ANSI_COLOR_YELLOW "Could not locate `%s' but found the gzip file `%s.gz'.\nRunning system command `" ANSI_COLOR_BLUE "%s"ANSI_COLOR_YELLOW"' now to uncompress"ANSI_COLOR_RESET "\n",filename,filename,buf);
                int status = run_system_call(buf);
                if(status != EXIT_SUCCESS) {
                    return -1;
                }
            }
        } else {
            //file exists -> nothing to do.
            fclose(fp);
        }
    }


    if(strncmp(format,"f",1)==0) { /*Read-in fast-food file*/
        //read fast-food file
        int idat[5];
        float fdat[9];
        size_t bytes=sizeof(int) + sizeof(fdat) + sizeof(int);//skip fdat
        bytes += sizeof(int) + sizeof(float) + sizeof(int); //skip znow

        FILE *fp = my_fopen(filename,"r");
        if(fp == NULL) {
            return -1;
        }
        int status = my_ftread(idat,sizeof(idat[0]),5,fp);
        if(status != EXIT_SUCCESS) {
            return -1;
        }
        np = (int64_t) idat[1]; //idat[1] is int.

        for(int i=0;i<num_fields;i++) {
            data[i] = my_malloc(size,np);
            if(data[i] == NULL) {
                for(int j=i-1;j>=0;j--) {
                    free(data[j]);
                }
                return -1;
            }
        }

        my_fseek(fp,bytes,SEEK_CUR);
        //Check that the file was written with the requested precision
        unsigned int dummy;
        my_fread(&dummy,sizeof(dummy), 1, fp);
        //so rewind by 4 bytes  prepare for calls to ftread
        my_fseek(fp, -4, SEEK_CUR);
        dummy /= np;
        XRETURN((dummy == 4 || dummy == 8), -1, "Data-type in file = %u must be either 4 byte (float) or 8 byte(double) precision", dummy);

        if(dummy == size) {
            for(int i=0;i<num_fields;i++) {
                status = my_ftread(data[i],size, np, fp);
                if(status != EXIT_SUCCESS) {
                    return -1;
                }
            }
        } else {
#ifndef SILENT
            fprintf(stderr,ANSI_COLOR_MAGENTA"WARNING: File was written in a different precision than requested (file precision = %u requested precision = %zu)"ANSI_COLOR_RESET"\n",dummy,size);
#endif
            //Okay so the file was written in a different precision.
            //First, print a warning message and then read-in correctly with the
            //requested precision
            if(dummy == 4) {
                XRETURN(size == 8, -1, "size = %zu should have been 8 (doubles were expected)\n", size);
                float *tmp = my_malloc(dummy,np);
                if(tmp == NULL) {
                    return -1;
                }

                //read-in the fields
                for(int i=0;i<num_fields;i++) {
                    status = my_ftread(tmp, dummy, np, fp);
                    if(status != EXIT_SUCCESS) {
                        return -1;
                    }
                    double *tmp_pos = (double *) data[i];
                    for(int64_t j=0;j<np;j++) tmp_pos[j] = tmp[j];
                }

                //free memory
                free(tmp);
            } else {
                XRETURN(size == 4, -1, "size = %zu should have been 4 (floats were expected)\n", size);
                double *tmp = my_malloc(dummy,np);
                if(tmp == NULL) {
                    return -1;
                }

                //read-in the fields
                for(int i=0;i<num_fields;i++) {
                    status = my_ftread(tmp, dummy, np, fp);
                    if(status != EXIT_SUCCESS) {
                        return -1;
                    }
                    float *tmp_pos = (float *) data[i];
                    for(int64_t j=0;j<np;j++) tmp_pos[j] = tmp[j];
                }
                //free memory
                free(tmp);
            }
        }

        fclose(fp);
    } else if(strncmp(format,"a",1)==0 || strncmp(format,"c",1)==0) { /* Read in ascii (white-space/comma) separated file*/
        int64_t i;
        int64_t nmax=300000;
        int nread;
        const int MAXBUFSIZE=10000;
        char buffer[MAXBUFSIZE];
        char delimiters[]=" ,\t";//delimiters are white-space, comma and tab

        for(i=0;i<num_fields;i++) {
            data[i] = my_malloc(size,nmax);
            if(data[i] == NULL) {
                for(int j=i-1;j>=0;j--) {
                    free(data[j]);
                }
                return -1;
            }
        }

        FILE *fp = my_fopen(filename,"r");
        if(fp == NULL) {
            return -1;
        }
        i = 0 ;
        while(fgets(buffer,MAXBUFSIZE,fp) != NULL) {
            double tmp;
            char *saveptr;
            int flag = 1;
            char *copy=buffer;
            for(int j=0;j<num_fields;j++,copy=NULL) {
                char *token = strtok_r(copy,delimiters,&saveptr);
                nread = sscanf(token,"%lf",&tmp);
                if(nread == 1) {
                    if(size==4) {
                        //data is float pointer
                        float *tmp_pos = (float *) data[j];
                        tmp_pos[i] = tmp;
                    } else {
                        //data is double pointer
                        double *tmp_pos = (double *) data[j];
                        tmp_pos[i] = tmp;
                    }
                } else {
                    flag = 0;
                }
            }
            if(flag == 1) {
                i++;
            } else {
                fprintf(stderr,ANSI_COLOR_YELLOW "io> WARNING: Could not parse all requested %d fields in line %"PRId64". Skipping that line" ANSI_COLOR_RESET "\n",num_fields,i);
            }

            if(i==nmax) {
                nmax *= MEMORY_INCREASE_FAC;
                while(nmax<=i)
                    nmax+=5;

                for(int j=0;j<num_fields;j++) {
                    char varname[20];
                    snprintf(varname,20,"data[%d]",j);
                    void *pos=NULL;
                    do {
                        pos = my_realloc(data[j],size,nmax,varname);
                        if(pos == NULL) {
                            nmax--;
                        }
                    } while(nmax > i && pos == NULL);

                    if(nmax == i) {
                        /*realloc failed. free memory and return */
                        fprintf(stderr,"In %s> Reallocation failed,  randomly subsampling the input particle set (currently at %"PRId64" particles) might help\n",
                                __FUNCTION__, nmax);
                        for(int k=0;k<num_fields;k++) {
                            free(data[k]);
                        }
                        return -1;
                    }
                    data[j] = pos;
                }
            }
        }
        np=i ;
        nmax=np;

        //release the extra memory.
        for(int j=0;j<num_fields;j++) {
            char varname[20];
            snprintf(varname,20,"data[%d]",j);
            data[j] = my_realloc(data[j],size,nmax,varname);
        }
        fclose(fp);
    } else {
        fprintf(stderr,"ERROR: In %s> Unknown format `%s'\n",__FUNCTION__,format);
        return -1;
    }

    return np;
}
