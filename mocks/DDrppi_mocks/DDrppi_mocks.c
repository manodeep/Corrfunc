/* File: DDrppi_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM DDrppi

--- DDrppi file1 format1 file2 format2 binfile pimax cosmology numthreads [weight_method weights_file1 weights_format1 [weights_file2 weights_format2]] > DDfile
--- Measure the cross-correlation function xi(rp,pi) for two different
   data files (or autocorrelation if data1=data2).
 * file1         = name of first data file
 * format1       = format of first data file  (a=ascii, c=csv, f=fast-food)
 * file2         = name of second data file
 * format2       = format of second data file (a=ascii, c=csv, f=fast-food)
 * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)
 * pimax         = maximum line-of-sight-separation
 * cosmology     = flag to pick-up the cosmology combination to use (set as an array of combinations in ../utils/cosmology_params.c)
 * numthreads    = number of threads to use
--- OPTIONAL ARGS:
 * weight_method = the type of pair weighting to apply.  Options are: 'pair_product', 'none'.  Default: 'none'.
 * weights_file1 = name of file containing the weights corresponding to the first data file
 * weights_format1 = format of file containing the weights corresponding to the first data file
 * weights_file2 = name of file containing the weights corresponding to the second data file
 * weights_format2 = format of file containing the weights corresponding to the second data file
---OUTPUT:
 > DDfile        = name of output file. Contains <npairs rpavg logrp pi weightavg>

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include "defs.h" //for basic API + all macros
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_rp_pi_mocks.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities


void Printhelp(void);

int main(int argc, char *argv[])
{
    /*---Arguments-------------------------*/
    char *file1=NULL,*file2=NULL, *weights_file1=NULL,*weights_file2=NULL;
    char *fileformat1=NULL,*fileformat2=NULL, *weights_fileformat1=NULL,*weights_fileformat2=NULL;
    char *binfile=NULL;
    char *weight_method_str=NULL;
    DOUBLE pimax;
    
    weight_method_t weight_method = NONE;
    int num_weights = 0;

    /*---Data-variables--------------------*/
    int64_t ND1,ND2 ;

    DOUBLE *thetaD1,*phiD1,*czD1, *weights1[MAX_NUM_WEIGHTS]={NULL};
    DOUBLE *thetaD2,*phiD2,*czD2, *weights2[MAX_NUM_WEIGHTS]={NULL};

    struct timeval t_end,t_start,t0,t1;
    double read_time=0.0;
    gettimeofday(&t_start,NULL);
    int nthreads=1;
    
    /*---Corrfunc-variables----------------*/
#if defined(_OPENMP)
    const char argnames[][30]={"file1","format1","file2","format2","binfile","pimax","cosmology flag","numthreads"};
#else
    const char argnames[][30]={"file1","format1","file2","format2","binfile","pimax","cosmology flag"};
#endif
    const char optargnames[][30]={"weight_method", "weights_file1","weights_format1","weights_file2","weights_format2"};
    
    int nargs=sizeof(argnames)/(sizeof(char)*30);
    int noptargs=sizeof(optargnames)/(sizeof(char)*30);
    
    int cosmology=1;

    /*---Read-arguments-----------------------------------*/
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

    pimax=40.0;
    sscanf(argv[6],"%"REAL_FORMAT,&pimax) ;
    cosmology = atoi(argv[7]);

#if defined(_OPENMP)
    nthreads=atoi(argv[8]);
    assert(nthreads >= 1 && "Number of threads must be at least 1");
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
    if(strcmp(file1,file2)==0) {
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
    ND1=read_positions(file1,fileformat1,sizeof(DOUBLE), 3, &phiD1, &thetaD1, &czD1);
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);
    gettimeofday(&t0,NULL);
    
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

    if (autocorr==0) {
        /*---Read-data2-file----------------------------------*/
        ND2=read_positions(file2,fileformat2,sizeof(DOUBLE), 3, &phiD2, &thetaD2, &czD2);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
        
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
        //None of these are required. But I prefer to preserve the possibility
        ND2 = ND1;
        thetaD2 = thetaD1;
        phiD2 = phiD1;
        czD2 = czD1;
        for(int w = 0; w < MAX_NUM_WEIGHTS; w++){
          weights2[w] = weights1[w];
        }
    }



    /*---Count-pairs--------------------------------------*/
    results_countpairs_mocks results;
    struct config_options options = get_config_options();

    /* Pack weights into extra options */
    struct extra_options extra = get_extra_options(weight_method);
    for(int w = 0; w < num_weights; w++){
        extra.weights0.weights[w] = (void *) weights1[w];
        extra.weights1.weights[w] = (void *) weights2[w];
    }
    
    int status = countpairs_mocks(ND1,phiD1,thetaD1,czD1,
                                  ND2,phiD2,thetaD2,czD2,
                                  nthreads,
                                  autocorr,
                                  binfile,
                                  pimax,
                                  cosmology,
                                  &results,
                                  &options,
                                  &extra);

    free(phiD1);free(thetaD1);free(czD1);
    for(int w = 0; w < num_weights; w++){
        free(weights1[w]);
    }
    if(autocorr == 0) {
        free(phiD2);free(thetaD2);free(czD2);
        for(int w = 0; w < num_weights; w++){
          free(weights2[w]);
        }
    }

    if(status != EXIT_SUCCESS) {
        return status;
    }

    const DOUBLE dpi = pimax/(DOUBLE)results.npibin ;
    const int npibin = results.npibin;
    for(int i=1;i<results.nbin;i++) {
        const double logrp = LOG10(results.rupp[i]);
        for(int j=0;j<npibin;j++) {
            const int index = i*(npibin+1) + j;
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf %20.8lf \n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi, results.weightavg[index]);
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
#if defined(_OPENMP)
    fprintf(stderr,"   --- DDrppi file1 format1 file2 format2 binfile pimax cosmology numthreads [weight_method weights_file1 weights_format1 [weights_file2 weights_format2]] > DDfile\n") ;
#else
    fprintf(stderr,"   --- DDrppi file1 format1 file2 format2 binfile pimax cosmology [weight_method weights_file1 weights_format1 [weights_file2 weights_format2]] > DDfile\n") ;
#endif
    fprintf(stderr,"   --- Measure the cross-correlation function xi(rp,pi) for two different\n") ;
    fprintf(stderr,"       data files (or autocorrelation if data1=data2).\n") ;
    fprintf(stderr,"     * data1         = name of first data file\n") ;
    fprintf(stderr,"     * format1       = format of first data file  (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * data2         = name of second data file\n") ;
    fprintf(stderr,"     * format2       = format of second data file (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
    fprintf(stderr,"     * pimax         = maximum line-of-sight-separation\n") ;
    fprintf(stderr,"     * cosmology     = flag to pick-up the cosmology combination to use (set as an array of combinations in ../utils/cosmology_params.c)\n") ;
#if defined(_OPENMP)
    fprintf(stderr,"     * numthreads    = number of threads to use\n");
#endif
    fprintf(stderr,"   --- OPTIONAL ARGS:\n");
    fprintf(stderr,"     * weight_method = the type of pair weighting to apply.  Options are: 'pair_product', 'none'.  Default: 'none'.\n");
    fprintf(stderr,"     * weights_file1 = name of file containing the weights corresponding to the first data file\n");
    fprintf(stderr,"     * weights_format1 = format of file containing the weights corresponding to the first data file\n");
    fprintf(stderr,"     * weights_file2 = name of file containing the weights corresponding to the second data file\n");
    fprintf(stderr,"     * weights_format2 = format of file containing the weights corresponding to the second data file\n");
    fprintf(stderr,"   ---OUTPUT:\n") ;
    fprintf(stderr,"     > DDfile        = name of output file. Contains <npairs rpavg logrp pi weightavg>\n") ;

    fprintf(stderr,"\n\tCompile options: \n");

#ifdef OUTPUT_RPAVG
    fprintf(stderr,"Output RPAVG = True\n");
#else
    fprintf(stderr,"Output RPAVG = False\n");
#endif

#ifdef FAST_DIVIDE
    fprintf(stderr,"Fast (approx) divide = True\n");
#else
    fprintf(stderr,"Fast (approx) divide = False\n");
#endif

#ifdef COMOVING_DIST
    fprintf(stderr,"CZ column contains co-moving distance = True\n");
#else
    fprintf(stderr,"CZ column contains co-moving distance = False\n");
#endif    
    
#ifdef DOUBLE_PREC
    fprintf(stderr,"Precision = double\n");
#else
    fprintf(stderr,"Precision = float\n");
#endif

#if defined(_OPENMP)
    fprintf(stderr,"Use OMP = True\n");
#else
    fprintf(stderr,"Use OMP = False\n");
#endif

    fprintf(stderr,"=========================================================================\n") ;
}
