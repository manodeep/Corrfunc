/* File: xi.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM xi
--- xi boxsize file format binfile Nthreads [weight_method weights_file weights_format] > xifile
--- Measure the auto-correlation function xi(r) for a single periodic box
 * boxsize      = BoxSize (in same units as X/Y/Z of the data)
 * file         = name of data file
 * format       = format of data file  (a=ascii, c=csv, f=fast-food)
 * binfile      = name of ascii file containing the r-bins (rmin rmax for each bin)
 * Nthreads    = number of threads to use
--- OPTIONAL ARGS:
 * weight_method = the type of pair weighting to apply.  Options are: 'pair_product', 'none'.  Default: 'none'.
 * weights_file = name of file containing the weights corresponding to the data file
 * weights_format = format of file containing the weights corresponding to the data file
---OUTPUT:
 > xifile        = name of output file. Contains <rmin rmax rpavg=0.0 npairs xi weightavg>
   ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_xi.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities



void Printhelp(void);

int main(int argc, char *argv[])
{

    /*---Arguments-------------------------*/
    double boxsize;
    char *file=NULL,*fileformat=NULL,*weights_file=NULL,*weights_fileformat=NULL;;
    char *binfile=NULL;
    char *weight_method_str=NULL;
    
    weight_method_t weight_method = NONE;
    int num_weights = 0;

    /*---Data-variables--------------------*/
    int64_t ND1=0;
    DOUBLE *x1=NULL,*y1=NULL,*z1=NULL,*weights1[MAX_NUM_WEIGHTS]={NULL};;

    int nthreads=1;

    /*---Corrfunc-variables----------------*/
#if defined(_OPENMP)
    const char argnames[][30]={"boxsize","file","format","binfile","Nthreads"};
#else
    const char argnames[][30]={"boxsize","file","format","binfile"};
#endif
    const char optargnames[][30]={"weight_method", "weights_file","weights_format"};

    int nargs=sizeof(argnames)/(sizeof(char)*30);
    int noptargs=sizeof(optargnames)/(sizeof(char)*30);

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
            else if(i <= nargs + noptargs)
                fprintf(stderr,"\t\t %s = `%s' \n",optargnames[i-1-nargs],argv[i]);
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
    
    /* Validate optional arguments */
    int noptargs_given = argc - (nargs + 1);
    if(noptargs_given != 0 && noptargs_given != 3){
        Printhelp();
        fprintf(stderr,"\nFound: %d optional arguments; must be 0 (no weights), 3 (for one set of weights)\n ", noptargs_given);
        int i;
        for(i=nargs+1;i<argc;i++) {
            if(i <= nargs + noptargs)
                fprintf(stderr,"\t\t %s = `%s' \n",optargnames[i-nargs-1],argv[i]);
            else
                fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
        }
        return EXIT_FAILURE;
    }

    boxsize=atof(argv[1]);
    file=argv[2];
    fileformat=argv[3];
    binfile=argv[4];

#if defined(USE_OMP) && defined(_OPENMP)
    nthreads=atoi(argv[5]);
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
      
       weights_file = argv[nargs + 2];
       weights_fileformat = argv[nargs + 3];
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


    gettimeofday(&t0,NULL);
    /*---Read-data1-file----------------------------------*/
    ND1=read_positions(file,fileformat, sizeof(DOUBLE), 3, &x1, &y1, &z1);
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);

    //check that theee positions are within limits
    for(int i=0;i<ND1;i++) {
        assert(x1[i] >= 0.0 && x1[i] <= boxsize && "xpos is within limits [0, boxsize]");
        assert(y1[i] >= 0.0 && y1[i] <= boxsize && "ypos is within limits [0, boxsize]");
        assert(z1[i] >= 0.0 && z1[i] <= boxsize && "zpos is within limits [0, boxsize]");
    }
    
    /* Read weights file */
    if(weights_file != NULL){
        gettimeofday(&t0,NULL);
        int64_t wND1 = read_columns_into_array(weights_file, weights_fileformat, sizeof(DOUBLE), num_weights, (void **) weights1);
        gettimeofday(&t1,NULL);
        read_time += ADD_DIFF_TIME(t0,t1);
      
        if(wND1 != ND1){
          fprintf(stderr, "Error: read %"PRId64" lines from %s, but read %"PRId64" from %s\n", wND1, weights_file, ND1, file);
          return EXIT_FAILURE;
        }
    }


    /*---Count-pairs--------------------------------------*/
    gettimeofday(&t0,NULL);
    results_countpairs_xi results;
    struct config_options options = get_config_options();
    
    /* Pack weights into extra options */
    struct extra_options extra = get_extra_options(weight_method);
    for(int w = 0; w < num_weights; w++){
        extra.weights0.weights[w] = (void *) weights1[w];
    }

    /* If you want to change the bin refine factors */
    /* const int bf[] = {2, 2, 1}; */
    /* set_bin_refine_factors(&options, bf); */
    int status = countpairs_xi(ND1, x1, y1, z1,
                               boxsize,
                               nthreads,
                               binfile,
                               &results,
                               &options,
                               &extra);
    free(x1);free(y1);free(z1);
    for(int w = 0; w < num_weights; w++){
        free(weights1[w]);
    }
    if(status != EXIT_SUCCESS) {
        return status;
    }

    gettimeofday(&t1,NULL);
    double pair_time = ADD_DIFF_TIME(t0,t1);

    //Output the results
    /* Note: we discard the first bin, to mimic the fact that close pairs
     * are disregarded in SDSS data.
     */
    DOUBLE rlow=results.rupp[0];
    for(int i=1;i<results.nbin;++i) {
        fprintf(stdout,"%e\t%e\t%e\t%12"PRIu64"\t%e\t%e\n",rlow,results.rupp[i],results.ravg[i],results.npairs[i],results.xi[i], results.weightavg[i]);
        rlow=results.rupp[i];
    }

    //free the memory in the results struct
    free_results_xi(&results);

    gettimeofday(&t_end,NULL);
    fprintf(stderr,"xi> Done -  ND1=%12"PRId64". Time taken = %6.2lf seconds. read-in time = %6.2lf seconds pair-counting time = %6.2lf sec\n",
            ND1,ADD_DIFF_TIME(t_start,t_end),read_time,pair_time);
    return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{
    fprintf(stderr,"=========================================================================\n") ;
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"   --- xi boxsize file format binfile Nthreads [weight_method weights_file weights_format] > xifile\n");
#else
    fprintf(stderr,"   --- xi boxsize file format binfile [weight_method weights_file weights_format] > xifile\n");
#endif
    fprintf(stderr,"   --- Measure the auto-correlation function xi(r) for a single periodic box\n") ;
    fprintf(stderr,"     * boxsize      = BoxSize (in same units as X/Y/Z of the data)\n") ;
    fprintf(stderr,"     * file         = name of data file\n") ;
    fprintf(stderr,"     * format       = format of data file  (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * binfile      = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"     * Nthreads    = number of threads to use\n");
#endif

    fprintf(stderr,"   --- OPTIONAL ARGS:\n");
    fprintf(stderr,"     * weight_method = the type of pair weighting to apply.  Options are: 'pair_product', 'none'.  Default: 'none'.\n");
    fprintf(stderr,"     * weights_file = name of file containing the weights corresponding to the data file\n");
    fprintf(stderr,"     * weights_format = format of file containing the weights corresponding to the data file\n");
    fprintf(stderr,"   ---OUTPUT:\n") ;
#ifdef OUTPUT_RPAVG
    fprintf(stderr,"     > xifile        = name of output file. Contains <rmin rmax rpavg npairs xi weightavg>\n") ;
#else
    fprintf(stderr,"     > xifile        = name of output file. Contains <rmin rmax rpavg=0.0 npairs xi weightavg>\n") ;
#endif

    fprintf(stderr,"\n\tCompile options: \n");
#ifdef PERIODIC
    fprintf(stderr,"\tPeriodic = True\n");
#else
    fprintf(stderr,"\tPeriodic = False\n");
#endif

#ifdef OUTPUT_RPAVG
    fprintf(stderr,"\tOutput RPAVG = True\n");
#else
    fprintf(stderr,"\tOutput RPAVG = False\n");
#endif

#ifdef DOUBLE_PREC
    fprintf(stderr,"\tPrecision = double\n");
#else
    fprintf(stderr,"\tPrecision = float\n");
#endif

#if defined(__AVX__)
    fprintf(stderr,"\tUse AVX = True\n");
#else
    fprintf(stderr,"\tUse AVX = False\n");
#endif

#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"\tUse OMP = True\n");
#else
    fprintf(stderr,"\tUse OMP = False\n");
#endif

    fprintf(stderr,"=========================================================================\n") ;


}
