/* File: DDrppi.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM DDrppi

   --- DDrppi rpmin rpmax nrpbin data1 data2 [pimax] > wpfile
   --- Measure the cross-correlation function xi(rp,pi) for two different
   data files (or autocorrelation if data1=data2).

   * rpmin   = inner radius of smallest bin (in Mpc/h)
   * rpmax   = outer radius of largest bin
   * nrpbin  = number of bins (logarithmically spaced in r)
   * data1   = name of first data file
   * data2   = name of second data file
   *[pimax]  = maximum value of line-of-sight separation (default=40 Mpc/h)
   > DDfile  = name of output file <logrp log(<rp>) pi pairs>
   ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_rp_pi.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities


void Printhelp(void);

int main(int argc, char *argv[])
{

    /*---Arguments-------------------------*/
    char *file1=NULL,*file2=NULL;
    char *fileformat1=NULL,*fileformat2=NULL;
    char *binfile=NULL;
    DOUBLE pimax ;


    /*---Data-variables--------------------*/
    int64_t ND1=0,ND2=0;

    DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;
    DOUBLE *x2=NULL,*y2=NULL,*z2=NULL;//will point to x1/y1/z1 in case of auto-corr


    /*---Corrfunc-variables----------------*/
#if !(defined(USE_OMP) && defined(_OPENMP))
    const char argnames[][30]={"file1","format1","file2","format2","binfile","pimax"};
#else
    int nthreads=2;
    const char argnames[][30]={"file1","format1","file2","format2","binfile","pimax","Nthreads"};
#endif
    int nargs=sizeof(argnames)/(sizeof(char)*30);

    struct timeval t_end,t_start,t0,t1;
    double read_time=0.0;
    gettimeofday(&t_start,NULL);

    /*---Read-arguments-----------------------------------*/
    if(argc< (nargs+1)) {
        Printhelp() ;
        fprintf(stderr,"\nFound: %d parameters\n ",argc-1);
        int i;
        for(i=1;i<argc;i++) {
            if(i <= nargs)
                fprintf(stderr,"\t\t %s = `%s' \n",argnames[i-1],argv[i]);
            else
                fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
        }
        if(i <= nargs) {
            fprintf(stderr,"\nMissing required parameters \n");
            for(i=argc;i<=nargs;i++)
                fprintf(stderr,"\t\t %s = `?'\n",argnames[i-1]);
        }
        return EXIT_FAILURE;
    }

    file1=argv[1];
    fileformat1=argv[2];
    file2=argv[3];
    fileformat2=argv[4];
    binfile=argv[5];

    pimax=40.0;
#ifdef DOUBLE_PREC
    sscanf(argv[6],"%lf",&pimax) ;
#else
    sscanf(argv[6],"%f",&pimax) ;
#endif


#if defined(USE_OMP) && defined(_OPENMP)
    nthreads=atoi(argv[7]);
    assert(nthreads >= 1 && "Number of threads must be at least 1");
#endif


    int autocorr=0;
    if(strcmp(file1,file2)==0) {
        autocorr=1;
    }

    fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
    fprintf(stderr,"\n\t\t -------------------------------------\n");
    for(int i=1;i<argc;i++) {
        if(i <= nargs) {
            fprintf(stderr,"\t\t %-10s = %s \n",argnames[i-1],argv[i]);
        }  else {
            fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
        }
    }
    fprintf(stderr,"\t\t -------------------------------------\n");


    gettimeofday(&t0,NULL);
    /*---Read-data1-file----------------------------------*/
    ND1=read_positions(file1,fileformat1,sizeof(DOUBLE), 3, &x1, &y1, &z1);
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);
    gettimeofday(&t0,NULL);

    if (autocorr==0) {
        /*---Read-data2-file----------------------------------*/
        ND2=read_positions(file2,fileformat2,sizeof(DOUBLE), 3, &x2, &y2, &z2);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);

    } else {
        //None of these are required. But I prefer to preserve the possibility
        ND2 = ND1;
        x2 = x1;
        y2 = y1;
        z2 = z1;
    }

    /*---Count-pairs--------------------------------------*/
    gettimeofday(&t0,NULL);
    results_countpairs_rp_pi *results = countpairs_rp_pi(ND1,x1,y1,z1,
                                                         ND2,x2,y2,z2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                         nthreads,
#endif
                                                         autocorr,
                                                         binfile,
                                                         pimax);



    gettimeofday(&t1,NULL);
    double pair_time = ADD_DIFF_TIME(t0,t1);
    free(x1);free(y1);free(z1);
    if(autocorr == 0) {
        free(x2);free(y2);free(z2);
    }

    const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
    const int npibin = results->npibin;
    for(int i=1;i<results->nbin;i++) {
        const double logrp = LOG10(results->rupp[i]);
        for(int j=0;j<npibin;j++) {
            int index = i*(npibin+1) + j;
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->rpavg[index],logrp,(j+1)*dpi);
        }
    }

    //free memory in results struct
    free_results_rp_pi(&results);

    gettimeofday(&t_end,NULL);
    fprintf(stderr,"xi_rp_pi> Done -  ND1=%12"PRId64" ND2=%12"PRId64". Time taken = %6.2lf seconds. read-in time = %6.2lf seconds pair-counting time = %6.2lf sec\n",
            ND1,ND2,ADD_DIFF_TIME(t_start,t_end),read_time,pair_time);
    return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{
    fprintf(stderr,"=========================================================================\n") ;
    fprintf(stderr,"   --- DDrppi file1 format1 file2 format2 pimax [Nthreads] > DDfile\n") ;
    fprintf(stderr,"   --- Measure the cross-correlation function xi(rp,pi) for two different\n") ;
    fprintf(stderr,"       data files (or autocorrelation if data1=data2).\n") ;
    fprintf(stderr,"     * data1         = name of first data file\n") ;
    fprintf(stderr,"     * format1       = format of first data file  (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * data2         = name of second data file\n") ;
    fprintf(stderr,"     * format2       = format of second data file (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
    fprintf(stderr,"     * pimax         = maximum line-of-sight-separation\n") ;
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"     * numthreads    = number of threads to use\n");
#endif
    fprintf(stderr,"     > DDfile        = name of output file. Contains <npairs rpavg logrp pi>\n") ;

    fprintf(stderr,"\n\tCompile options: \n");
#ifdef PERIODIC
    fprintf(stderr,"Periodic = True\n");
#else
    fprintf(stderr,"Periodic = False\n");
#endif

#ifdef OUTPUT_RPAVG
    fprintf(stderr,"Output RPAVG = True\n");
#else
    fprintf(stderr,"Output RPAVG = False\n");
#endif

#ifdef DOUBLE_PREC
    fprintf(stderr,"Precision = double\n");
#else
    fprintf(stderr,"Precision = float\n");
#endif

#if defined(USE_AVX) && defined(__AVX__)
    fprintf(stderr,"Use AVX = True\n");
#else
    fprintf(stderr,"Use AVX = False\n");
#endif

#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"Use OMP = True\n");
#else
    fprintf(stderr,"Use OMP = False\n");
#endif

    fprintf(stderr,"=========================================================================\n") ;
}
