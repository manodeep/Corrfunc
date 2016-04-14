/* File: DDtheta_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_theta_mocks.h" //function proto-type for countpairs_theta_mocks
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities

void Printhelp(void);

int main(int argc, char **argv)
{
    /*---Arguments-------------------------*/
    char *file1=NULL,*file2=NULL;
    char *fileformat1=NULL,*fileformat2=NULL;
    char *binfile=NULL;

#if !(defined(USE_OMP) && defined(_OPENMP))
    const char argnames[][30]={"file1","format1","file2","format2","binfile"};
#else
    int nthreads;
    const char argnames[][30]={"file1","format1","file2","format2","binfile","Nthreads"};
#endif
    int nargs=sizeof(argnames)/(sizeof(char)*30);
    struct timeval tstart,t0,t1;
    double pair_time=0,read_time=0.0;

    gettimeofday(&tstart,NULL);

    /*---Data-variables--------------------*/
    int64_t ND1,ND2 ;
    DOUBLE *thetaD1,*phiD1 ;
    DOUBLE *thetaD2,*phiD2 ;

    /*---Check-arguments-----------------------------------*/
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
#if defined(USE_OMP) && defined(_OPENMP)
    nthreads=atoi(argv[6]);
    assert(nthreads >= 1 && "Number of threads must be at least 1");
#endif

    int autocorr=0;
    if( strcmp(file1,file2)==0) {
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
    ND1=read_positions(file1,fileformat1, sizeof(DOUBLE), 2, &phiD1,&thetaD1);
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);
    fprintf(stderr,"DDtheta_mocks> Read in %"PRId64" particles from file `%s'\n",ND1,file1);

    /*---Read-data2-file----------------------------------*/
    if(autocorr==0) {
        gettimeofday(&t0,NULL);
        ND2=read_positions(file2,fileformat2, sizeof(DOUBLE), 2, &phiD2, &thetaD2);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
        fprintf(stderr,"DDtheta_mocks> Read in %"PRId64" particles from file `%s'\n",ND2,file2);
    } else {
        ND2=ND1;
        thetaD2 = thetaD1;
        phiD2 = phiD1;
    }

    /*---Count-pairs--------------------------------------*/
    gettimeofday(&t0,NULL);
    results_countpairs_theta *results = countpairs_theta_mocks(ND1,phiD1,thetaD1,
                                                               ND2,phiD2,thetaD2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                               nthreads,
#endif
                                                               autocorr,
                                                               binfile) ;

    gettimeofday(&t1,NULL);
    pair_time = ADD_DIFF_TIME(t0,t1);
    free(thetaD1);free(phiD1);
    if(autocorr==0) {
        free(thetaD2);free(phiD2);
    }

    /*---Output-Pairs-------------------------------------*/
    DOUBLE theta_low = results->theta_upp[0];
    for(int i=1;i<results->nbin;i++) {
        fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->theta_avg[i],theta_low,results->theta_upp[i]);
        theta_low=results->theta_upp[i];
    }

    free_results_countpairs_theta(&results);
    gettimeofday(&t1,NULL);
    fprintf(stderr,"DDtheta> Done -  ND1=%"PRId64" ND2=%"PRId64". Time taken = %6.2lf seconds. read-in time = %6.2lf seconds sec pair-counting time = %6.2lf sec\n",
            ND1,ND2,ADD_DIFF_TIME(tstart,t1),read_time,pair_time);
    return EXIT_SUCCESS;

}



/*---Print-help-information---------------------------*/
void Printhelp(void)
{
    fprintf(stderr,"========================================================================================\n") ;
    fprintf(stderr,"   --- DDtheta file1 format1 file2 format2 binfile [Nthreads] > Thetafile\n") ;
    fprintf(stderr,"   --- Measure the cross-correlation function w(theta) for two different\n") ;
    fprintf(stderr,"       data files (or autocorrelation if data1=data2).\n") ;
    fprintf(stderr,"     * file1       = name of first data file\n") ;
    fprintf(stderr,"     * format1     = format of first data file (f for fast-food, a for ascii)\n") ;
    fprintf(stderr,"     * file2       = name of second data file\n") ;
    fprintf(stderr,"     * format2     = format of second data file\n") ;
    fprintf(stderr,"     * binfile     = name of ascii file containing the theta-bins (thetamin thetamax for each bin)\n") ;
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"     * numthreads  = number of threads to use\n");
#endif

#ifdef OUTPUT_THETAAVG
    fprintf(stderr,"     > Thetafile   = name of output file <npairs thetaavg thetamin thetamax>\n") ;
#else
    fprintf(stderr,"     > Thetafile   = name of output file <npairs    0.0   thetamin thetamax>\n") ;
#endif
    fprintf(stderr,"========================================================================================\n") ;
    fprintf(stderr,"\n\tCompile options: \n");
#ifdef LINK_IN_DEC
    fprintf(stderr,"LINK_IN_DEC = True\n");
#else
    fprintf(stderr,"LINK_IN_DEC = False\n");
#endif

#ifdef LINK_IN_RA
    fprintf(stderr,"LINK_IN_RA = True\n");
#else
    fprintf(stderr,"LINK_IN_RA = False\n");
#endif

#ifdef OUTPUT_THETAAVG
    fprintf(stderr,"Output THETAAVG = True\n");
#else
    fprintf(stderr,"Output THETAAVG = False\n");
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
