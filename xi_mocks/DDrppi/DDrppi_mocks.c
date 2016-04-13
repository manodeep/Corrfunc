/* File: DDrppi_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM DDrppi

   --- DDrppi file1 format1 file2 format2 pimax [Nthreads] > DDfile
   --- Measure the cross-correlation function xi(rp,pi) for two different
   data files (or autocorrelation if data1=data2).
   * data1         = name of first data file
   * format1       = format of first data file  (a=ascii, c=csv, f=fast-food)
   * data2         = name of second data file
   * format2       = format of second data file (a=ascii, c=csv, f=fast-food)
   * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)
   * pimax         = maximum line-of-sight-separation
   * cosmology     = flag to pick-up the cosmology combination to use (set as an array of combinations in ../utils/cosmology_params.c)
   * numthreads    = number of threads to use
   > DDfile        = name of output file. Contains <npairs rpavg logrp pi>

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_rp_pi_mocks.h" //function proto-type for countpairs
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
    int64_t ND1,ND2 ;

    DOUBLE *thetaD1,*phiD1,*czD1;
    DOUBLE *thetaD2,*phiD2,*czD2;

    struct timeval t_end,t_start,t0,t1;
    double read_time=0.0;
    gettimeofday(&t_start,NULL);

    /*---Corrfunc-variables----------------*/
#if !(defined(USE_OMP) && defined(_OPENMP))
    const char argnames[][30]={"file1","format1","file2","format2","binfile","pimax","cosmology flag"};
#else
    int nthreads=2;
    const char argnames[][30]={"file1","format1","file2","format2","binfile","pimax","cosmology flag","Nthreads"};
#endif
    int nargs=sizeof(argnames)/(sizeof(char)*30);
    int cosmology=1;

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
    sscanf(argv[6],"%"DOUBLE_FORMAT,&pimax) ;
    cosmology = atoi(argv[7]);

#if defined(USE_OMP) && defined(_OPENMP)
    nthreads=atoi(argv[8]);
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


    /*---Read-data1-file----------------------------------*/
    gettimeofday(&t0,NULL);
    ND1=read_positions(file1,fileformat1,sizeof(DOUBLE), 3, &phiD1, &thetaD1, &czD1);
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);
    gettimeofday(&t0,NULL);

    if (autocorr==0) {
        /*---Read-data2-file----------------------------------*/
        ND2=read_positions(file2,fileformat2,sizeof(DOUBLE), 3, &phiD2, &thetaD2, &czD2);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
    } else {
        //None of these are required. But I prefer to preserve the possibility
        ND2 = ND1;
        thetaD2 = thetaD1;
        phiD2 = phiD1;
        czD2 = czD1;
    }



    /*---Count-pairs--------------------------------------*/
    results_countpairs_mocks *results  = countpairs_mocks(ND1,phiD1,thetaD1,czD1,
                                                          ND2,phiD2,thetaD2,czD2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                          nthreads,
#endif
                                                          autocorr,
                                                          binfile,
                                                          pimax,
                                                          cosmology);

    free(phiD1);free(thetaD1);free(czD1);
    if(autocorr == 0) {
        free(phiD2);free(thetaD2);free(czD2);
    }


    const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
    const int npibin = results->npibin;
    for(int i=1;i<results->nbin;i++) {
        const double logrp = LOG10(results->rupp[i]);
        for(int j=0;j<npibin;j++) {
            const int index = i*(npibin+1) + j;
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->rpavg[index],logrp,(j+1)*dpi);
        }
    }

    free_results_mocks(&results);
    gettimeofday(&t_end,NULL);
    fprintf(stderr,"DDrppi> Done -  ND1=%"PRId64" ND2=%"PRId64". Time taken = %6.2lf seconds, read-in time = %6.2lf seconds \n",ND1,ND2,ADD_DIFF_TIME(t_start,t_end),read_time);
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
    fprintf(stderr,"     * cosmology     = flag to pick-up the cosmology combination to use (set as an array of combinations in ../utils/cosmology_params.c)\n") ;
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"     * numthreads    = number of threads to use\n");
#endif
    fprintf(stderr,"     > DDfile        = name of output file. Contains <npairs rpavg logrp pi>\n") ;

    fprintf(stderr,"\n\tCompile options: \n");

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
