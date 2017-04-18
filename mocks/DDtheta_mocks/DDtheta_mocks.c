/* File: DDtheta_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/*
--- DDtheta file1 format1 file2 format2 binfile numthreads [weight_method weights_file1 weights_format1 [weights_file2 weights_format2]] > Thetafile
--- Measure the cross-correlation function w(theta) for two different
   data files (or autocorrelation if data1=data2).
 * file1       = name of first data file
 * format1     = format of first data file (f for fast-food, a for ascii)
 * file2       = name of second data file
 * format2     = format of second data file
 * binfile     = name of ascii file containing the theta-bins (thetamin thetamax for each bin)
 * numthreads  = number of threads to use
--- OPTIONAL ARGS:
 * weight_method = the type of pair weighting to apply.  Options are: 'pair_product', 'none'.  Default: 'none'.
 * weights_file1 = name of file containing the weights corresponding to the first data file
 * weights_format1 = format of file containing the weights corresponding to the first data file
 * weights_file2 = name of file containing the weights corresponding to the second data file
 * weights_format2 = format of file containing the weights corresponding to the second data file
---OUTPUT:
 > Thetafile   = name of output file <thetamin thetamax thetaavg=0.0 npairs weightavg >
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_theta_mocks.h" //function proto-type for countpairs_theta_mocks
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities

void Printhelp(void);

int main(int argc, char **argv)
{
    /*---Arguments-------------------------*/
    char *file1=NULL,*file2=NULL, *weights_file1=NULL,*weights_file2=NULL;
    char *fileformat1=NULL,*fileformat2=NULL, *weights_fileformat1=NULL,*weights_fileformat2=NULL;
    char *binfile=NULL;
    char *weight_method_str=NULL;
    int nthreads=1;
    
    weight_method_t weight_method = NONE;
    int num_weights = 0;
    
#if defined(_OPENMP)
    const char argnames[][30]={"file1","format1","file2","format2","binfile","Nthreads"};
#else
    const char argnames[][30]={"file1","format1","file2","format2","binfile"};
#endif
    const char optargnames[][30]={"weight_method", "weights_file1","weights_format1","weights_file2","weights_format2"};
    
    int nargs=sizeof(argnames)/(sizeof(char)*30);
    int noptargs=sizeof(optargnames)/(sizeof(char)*30);
    
    struct timeval tstart,t0,t1;
    double pair_time=0,read_time=0.0;

    gettimeofday(&tstart,NULL);

    /*---Data-variables--------------------*/
    int64_t ND1,ND2 ;
    DOUBLE *thetaD1,*phiD1, *weights1[MAX_NUM_WEIGHTS]={NULL};
    DOUBLE *thetaD2,*phiD2, *weights2[MAX_NUM_WEIGHTS]={NULL};

    /*---Check-arguments-----------------------------------*/
    if(argc< (nargs+1)) {
        Printhelp() ;
        fprintf(stderr,"\nFound: %d parameters\n ",argc-1);
        int i;
        for(i=1;i<argc;i++) {
            if(i <= nargs)
                fprintf(stderr,"\t\t %s = `%s' \n",argnames[i-1],argv[i]);
            else if(i <= nargs + noptargs)
                fprintf(stderr,"\t\t %s = `%s' \n",optargnames[i-1-nargs],argv[i]);
            else
                fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
        }
        fprintf(stderr,"\nMissing required parameters \n");
        for(i=argc;i<=nargs;i++)
            fprintf(stderr,"\t\t %s = `?'\n",argnames[i-1]);
        return EXIT_FAILURE;
    }
    
    /* Validate optional arguments */
    int noptargs_given = argc - (nargs + 1);
    if(noptargs_given != 0 && noptargs_given != 3 && noptargs_given != 5){
        Printhelp();
        fprintf(stderr,"\nFound: %d optional arguments; must be 0 (no weights), 3 (for one set of weights) or 5 (for two sets)\n ", noptargs_given);
        int i;
        for(i=nargs+1;i<argc;i++) {
            if(i <= nargs + noptargs)
                fprintf(stderr,"\t\t %s = `%s' \n",optargnames[i-nargs-1],argv[i]);
            else
                fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
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
    if(nthreads < 1) {
        fprintf(stderr,"Warning: Nthreads = %d should be >=1. Setting nthreads to 1.\n", nthreads);
        nthreads = 1;
    }
#endif

    if(noptargs_given >= 3){
       weight_method_str = argv[nargs + 1];
       int wstatus = get_weight_method_by_name(weight_method_str, &weight_method);
       if(wstatus != EXIT_SUCCESS){
         fprintf(stderr, "Error: Unknown weight method \"%s\"\n", weight_method_str);
         return EXIT_FAILURE;
       }
       num_weights = get_num_weights_by_method(weight_method);
      
       weights_file1 = argv[nargs + 2];
       weights_fileformat1 = argv[nargs + 3];
    }
    if(noptargs_given >= 5){
       weights_file2 = argv[nargs + 4];
       weights_fileformat2 = argv[nargs + 5];
    }

    int autocorr=0;
    if( strcmp(file1,file2)==0) {
        autocorr=1;
    }

    fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
    fprintf(stderr,"\n\t\t -------------------------------------\n");
    for(int i=1;i<argc;i++) {
        if(i <= nargs) {
            fprintf(stderr,"\t\t %-10s = %s \n",argnames[i-1],argv[i]);
        } else if(i <= nargs + noptargs){
            fprintf(stderr,"\t\t %-10s = %s \n",optargnames[i-nargs-1],argv[i]);
        } else {
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
    
    /* Read weights file 1 */
    if(weights_file1 != NULL){
        gettimeofday(&t0,NULL);
        int64_t wND1 = read_columns_into_array(weights_file1,weights_fileformat1, sizeof(DOUBLE), num_weights, (void **) weights1);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
      
        if(wND1 != ND1){
          fprintf(stderr, "Error: read %"PRId64" lines from %s, but read %"PRId64" from %s\n", wND1, weights_file1, ND1, file1);
          return EXIT_FAILURE;
        }
    }

    /*---Read-data2-file----------------------------------*/
    if(autocorr==0) {
        gettimeofday(&t0,NULL);
        ND2=read_positions(file2,fileformat2, sizeof(DOUBLE), 2, &phiD2, &thetaD2);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
        fprintf(stderr,"DDtheta_mocks> Read in %"PRId64" particles from file `%s'\n",ND2,file2);
        
        if(weights_file2 != NULL){
            gettimeofday(&t0,NULL);
            int64_t wND2 = read_columns_into_array(weights_file2,weights_fileformat2, sizeof(DOUBLE), num_weights, (void **) weights2);
            gettimeofday(&t1,NULL);
            read_time += ADD_DIFF_TIME(t0,t1);

            if(wND2 != ND2){
              fprintf(stderr, "Error: read %"PRId64" lines from %s, but read %"PRId64" from %s\n", wND2, weights_file2, ND2, file2);
              return EXIT_FAILURE;
            }
        }
    } else {
        ND2=ND1;
        thetaD2 = thetaD1;
        phiD2 = phiD1;
        for(int w = 0; w < MAX_NUM_WEIGHTS; w++){
          weights2[w] = weights1[w];
        }
    }

    /*---Count-pairs--------------------------------------*/
    gettimeofday(&t0,NULL);
    results_countpairs_theta results;
    
    /* Pack weights into extra options */
    struct extra_options extra = get_extra_options(weight_method);
    for(int w = 0; w < num_weights; w++){
        extra.weights0.weights[w] = (void *) weights1[w];
        extra.weights1.weights[w] = (void *) weights2[w];
    }
    
    struct config_options options = get_config_options();
    int status = countpairs_theta_mocks(ND1,phiD1,thetaD1,
                                        ND2,phiD2,thetaD2,
                                        nthreads,
                                        autocorr,
                                        binfile,
                                        &results,
                                        &options,
                                        &extra);

    gettimeofday(&t1,NULL);
    pair_time = ADD_DIFF_TIME(t0,t1);
    free(thetaD1);free(phiD1);
    for(int w = 0; w < num_weights; w++){
        free(weights1[w]);
    }
    if(autocorr==0) {
        free(thetaD2);free(phiD2);
        for(int w = 0; w < num_weights; w++){
          free(weights2[w]);
        }
    }
    if(status != EXIT_SUCCESS) {
        return status;
    }
    
    /*---Output-Pairs-------------------------------------*/
    DOUBLE theta_low = results.theta_upp[0];
    for(int i=1;i<results.nbin;i++) {
        fprintf(stdout,"%20.8lf %20.8lf %20.8lf %10"PRIu64" %20.8lf\n",theta_low,results.theta_upp[i], results.theta_avg[i], results.npairs[i], results.weightavg[i]);
        theta_low=results.theta_upp[i];
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
#if defined(_OPENMP)
    fprintf(stderr,"   --- DDtheta file1 format1 file2 format2 binfile numthreads [weight_method weights_file1 weights_format1 [weights_file2 weights_format2]] > Thetafile\n") ;
#else
    fprintf(stderr,"   --- DDtheta file1 format1 file2 format2 binfile [weight_method weights_file1 weights_format1 [weights_file2 weights_format2]] > Thetafile\n") ;
#endif    
    fprintf(stderr,"   --- Measure the cross-correlation function w(theta) for two different\n") ;
    fprintf(stderr,"       data files (or autocorrelation if data1=data2).\n") ;
    fprintf(stderr,"     * file1       = name of first data file\n") ;
    fprintf(stderr,"     * format1     = format of first data file (f for fast-food, a for ascii)\n") ;
    fprintf(stderr,"     * file2       = name of second data file\n") ;
    fprintf(stderr,"     * format2     = format of second data file\n") ;
    fprintf(stderr,"     * binfile     = name of ascii file containing the theta-bins (thetamin thetamax for each bin)\n") ;
#if defined(_OPENMP)
    fprintf(stderr,"     * numthreads  = number of threads to use\n");
#endif
    fprintf(stderr,"   --- OPTIONAL ARGS:\n");
    fprintf(stderr,"     * weight_method = the type of pair weighting to apply.  Options are: 'pair_product', 'none'.  Default: 'none'.\n");
    fprintf(stderr,"     * weights_file1 = name of file containing the weights corresponding to the first data file\n");
    fprintf(stderr,"     * weights_format1 = format of file containing the weights corresponding to the first data file\n");
    fprintf(stderr,"     * weights_file2 = name of file containing the weights corresponding to the second data file\n");
    fprintf(stderr,"     * weights_format2 = format of file containing the weights corresponding to the second data file\n");
    fprintf(stderr,"   ---OUTPUT:\n") ;
#ifdef OUTPUT_THETAAVG
    fprintf(stderr,"     > Thetafile   = name of output file < thetamin thetamax thetaavg npairs weightavg >\n") ;
#else
    fprintf(stderr,"     > Thetafile   = name of output file <thetamin thetamax thetaavg=0.0 npairs weightavg >\n") ;
#endif
    fprintf(stderr,"========================================================================================\n") ;
    fprintf(stderr,"\n\tCompile options: \n");
#ifdef OUTPUT_THETAAVG
    fprintf(stderr,"Output THETAAVG = True\n");
#else
    fprintf(stderr,"Output THETAAVG = False\n");
#endif

/* #OPT += -DOUTPUT_THETAAVG */
/* OPT += -DLINK_IN_RA #link_in_dec must be enabled before link_in_ra */
/* #OPT += -DFAST_ACOS ## replaces acos by an 8th order REMEZ polynomial. Results are approximate (max. absolute error 3.6e-9) ~50% boost, but obviously approximate */

#ifdef LINK_IN_DEC
    fprintf(stderr,"Linking in declination = True\n");
#else
    fprintf(stderr,"Linking in declination = False\n");
#endif    

#ifdef LINK_IN_DEC
    //RA linking only works when dec linking is enabled
#ifdef LINK_IN_RA
    fprintf(stderr,"Linking in right ascension  = True\n");
#else
    fprintf(stderr,"Linking in right ascension  = False\n");
#endif//ra
#endif//dec

#ifdef OUTPUT_THETAAVG
    //Only makes sense when <theta> is requested
#ifdef FAST_ACOS
    fprintf(stderr,"Fast (approx) arc-cosine  = True\n");
#else
    fprintf(stderr,"Fast (approx) arc-cosine  = False\n");
#endif//fast acos    
#endif//thetaavg    
    
#ifdef DOUBLE_PREC
    fprintf(stderr,"Precision = double\n");
#else
    fprintf(stderr,"Precision = float\n");
#endif

#if defined(__AVX__)
    fprintf(stderr,"Use AVX = True\n");
#else
    fprintf(stderr,"Use AVX = False\n");
#endif

#if defined(_OPENMP)
    fprintf(stderr,"Use OMP = True\n");
#else
    fprintf(stderr,"Use OMP = False\n");
#endif

    fprintf(stderr,"=========================================================================\n") ;

}
