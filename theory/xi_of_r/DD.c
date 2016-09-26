/* File: DD.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities

void Printhelp(void);

int main(int argc, char *argv[])
{

    /*---Data-variables--------------------*/
    int64_t ND1,ND2 ;

    DOUBLE *x1=NULL,*y1=NULL,*z1=NULL,*weight1=NULL;
    DOUBLE *x2=NULL,*y2=NULL,*z2=NULL,*weight2=NULL;
    int num_fields = 3; // minimum

    char *file1=NULL,*file2=NULL;
    char *fileformat1=NULL,*fileformat2=NULL;
    char *binfile=NULL;
    weight_method_t weight_type = NONE;

    /*---Corrfunc-variables----------------*/
#if !(defined(USE_OMP) && defined(_OPENMP))
    const int nthreads=1;
    const char argnames[][30]={"file1","format1","file2","format2","binfile","[weight_method]"};
#else
    int nthreads=2;
    const char argnames[][30]={"file1","format1","file2","format2","binfile","Nthreads","[weight_method]"};
#endif
    int nargs=sizeof(argnames)/(sizeof(char)*30);
    int n_mandatory_args=nargs-1;

    struct timeval t_end,t_start,t0,t1;
    double read_time=0.0;

    gettimeofday(&t_start,NULL);

    /*---Read-arguments-----------------------------------*/
    if(argc < (n_mandatory_args+1)) {
        Printhelp() ;
        fprintf(stderr,"\nFound: %d parameters\n ",argc-1);
        int i;
        for(i=1;i<argc;i++) {
            if(i <= n_mandatory_args)
                fprintf(stderr,"\t\t %s = `%s' \n",argnames[i-1],argv[i]);
            else
                fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
        }
        if(i <= n_mandatory_args) {
            fprintf(stderr,"\nMissing required parameters \n");
            for(i=argc;i<=n_mandatory_args;i++)
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

    int n_weights = 0;
    if(argc > (n_mandatory_args+1)){
        int status = get_weight_method_by_name(argv[n_mandatory_args+2], &weight_type, &n_weights);
        if(status != EXIT_SUCCESS){
            return status;
        }
        num_fields += n_weights;
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
    // No good way to read a variable number of weights fields without modifying read_positions
    ND1=read_positions(file1,fileformat1, sizeof(DOUBLE), num_fields, &x1, &y1, &z1, &weight1);
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);

    int autocorr=0;
    if( strcmp(file1,file2)==0) {
        autocorr=1;
    }

    gettimeofday(&t0,NULL);
    if (autocorr==0) {
        /*---Read-data2-file----------------------------------*/
        ND2=read_positions(file2,fileformat2, sizeof(DOUBLE), num_fields, &x2, &y2, &z2, &weight2);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
    } else {
        ND2 = ND1;
        x2 = x1;
        y2 = y1;
        z2 = z1;
        weight2 = weight1;
    }

    /*---Count-pairs--------------------------------------*/
    gettimeofday(&t0,NULL);
    
    struct config_options options = get_config_options();
    options.float_type = sizeof(DOUBLE);
    
    struct extra_options extra;
    int status = get_extra_options(&extra, n_weights);
    if(status != EXIT_SUCCESS){
        return status;
    }
    extra.weighting_func_type = weight_method;
    extra.weights[0] = weight1;
    extra.weights[1] = weight2;
    
    results_countpairs results;
    status = countpairs(ND1,x1,y1,z1,
                            ND2,x2,y2,z2,
                            nthreads,
                            autocorr,
                            binfile,
                            &results,
                            &options,
                            &extra);

    
    free(x1);free(y1);free(z1);free(weight1);
    if(autocorr == 0) {
        free(x2);free(y2);free(z2);free(weight2);
    }
    if(status != EXIT_SUCCESS) {
      return status;
    }
    
    gettimeofday(&t1,NULL);
    double pair_time = ADD_DIFF_TIME(t0,t1);

    DOUBLE rlow=results.rupp[0];
    for(int i=1;i<results.nbin;i++) {
        fprintf(stdout,"%e\t%e\t%e\t%12"PRIu64"\n",rlow, results.rupp[i], results.rpavg[i],results.npairs[i]);
        rlow=results.rupp[i];
    }

    //free the memory in the results structx
    free_results(&results);
    gettimeofday(&t_end,NULL);
    fprintf(stderr,"xi_of_r> Done -  ND1=%"PRId64" ND2=%"PRId64". Time taken = %6.2lf seconds. read-in time = %6.2lf seconds sec pair-counting time = %6.2lf sec\n",
            ND1,ND2,ADD_DIFF_TIME(t_start,t_end),read_time,pair_time);
    return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{
    fprintf(stderr,"=========================================================================\n") ;
    fprintf(stderr,"   --- DD file1 format1 file2 format2 binfile [weight_type] > DDfile\n") ;
    fprintf(stderr,"   --- Measure the cross-correlation function DD(r) for two different\n") ;
    fprintf(stderr,"       data files (or autocorrelation if file1=file2).\n") ;
    fprintf(stderr,"     * file1         = name of first data file\n") ;
    fprintf(stderr,"     * format1       = format of first data file  (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * file2         = name of second data file\n") ;
    fprintf(stderr,"     * format2       = format of second data file (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"     * numthreads    = number of threads to use\n");
#endif
    fprintf(stderr,"     * [weight_type] = (optional) type of weighting to apply (n=none, p=pair-product). Default: n. \n");
    fprintf(stderr,"                       pair-product requires one weight column to be present in each data file \n");

#ifdef OUTPUT_RPAVG
    fprintf(stderr,"     > DDfile        = name of output file <npairs rpavg rmin rmax>\n") ;
#else
    fprintf(stderr,"     > DDfile        = name of output file <npairs  0.0  rmin rmax>\n") ;
#endif
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
