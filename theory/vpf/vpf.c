/* File: vpf.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM VPF

--- vpf rmax nbins ncentres num_pN file format seed boxsize > VPFfile
--- Measures the counts-in-spheres in a simulation box
 * rmax         = size of the biggest sphere
 * nbins        = number of bins to use for the counts-in-spheres
 * ncentres     = number of random spheres to place on the cube
 * num_pN        = number of counts-in-spheres to output. [num_pN=1-> P0, num_pN=2->P0,P1, num_pN=3->P0,P1,P2...
 * file         = name of data file
 * format       = format of data file  (a=ascii, c=csv, f=fast-food)
 * seed         = seed for random number generator
 * boxsize      = if periodic, the boxsize to use for the periodic wrap (0 means detect the particle extent)
 > VPFfile      = name of output file <r P0 P1 P2 ...>
   */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countspheres.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities
#include "macros.h"

void Printhelp(void);


int main(int argc, char *argv[])
{
    double rmax;
    int nbin,nc,num_pN;
    unsigned long seed=-1;
    char *file=NULL,*fileformat=NULL;

    double boxsize = -1.;

    /*---VPF-variables----------------*/
    const char argnames[][30]={"rmax","nbins","ncentres","num_pN","file","format","seed","boxsize"};
    int nargs=sizeof(argnames)/(sizeof(char)*30);

    int64_t np;
    DOUBLE *x,*y,*z;

    struct timeval tstart,t0,t1;
    gettimeofday(&tstart,NULL);

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

    rmax = atof(argv[1]);
    nbin = atoi(argv[2]);
    nc   = atoi(argv[3]);
    num_pN = atoi(argv[4]);
    file =  argv[5];
    fileformat = argv[6];
    seed = atol(argv[7]);
    boxsize = atof(argv[8]);

    assert(nbin >=1 && "Number of bins has to be at least 1");
    assert(nc >=1   && "Number of spheres has to be at least 1");
    assert(num_pN >= 1 && "Number of pN's requested must be at least 1");

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

    /*---Read-particle-file-----------------------------------------------------*/
    gettimeofday(&t0,NULL);
    np = read_positions(file,fileformat, sizeof(DOUBLE), 3, &x, &y, &z);
    gettimeofday(&t1,NULL);

    results_countspheres results;
    struct config_options options = get_config_options();
    options.boxsize = boxsize;

    /* If you want to change the bin refine factors */
    /* const int bf[] = {2, 2, 1}; */
    /* set_bin_refine_factors(&options, bf); */

    int status = countspheres(np, x, y, z,
                              rmax, nbin, nc,
                              num_pN,
                              seed,
                              &results,
                              &options,
                              NULL);
    free(x);free(y);free(z);
    if(status != EXIT_SUCCESS) {
        return status;
    }

    //Output the results
    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    for(int ibin=0;ibin<results.nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        fprintf(stdout,"%lf ", r);
        for(int i=0;i<num_pN;i++) {
            fprintf(stdout," %10.4e", (results.pN)[ibin][i]);
        }
        fprintf(stdout,"\n");
    }

    free_results_countspheres(&results);
    gettimeofday(&t1,NULL);
    fprintf(stderr,"vpf> Done. Ngal = %"PRId64". Time taken = %0.2lf seconds \n",np,ADD_DIFF_TIME(tstart,t1));

    return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{

    fprintf(stderr,"=========================================================================\n") ;
    fprintf(stderr,"   --- vpf rmax nbins nspheres numpN file format seed > VPFfile\n") ;
    fprintf(stderr,"   --- Measures the counts-in-spheres in a simulation box\n") ;
    fprintf(stderr,"     * rmax         = size of the biggest sphere\n");
    fprintf(stderr,"     * nbins        = number of bins to use for the counts-in-spheres\n");
    fprintf(stderr,"     * nspheres     = number of random spheres to place on the cube\n");
    fprintf(stderr,"     * numpN        = number of counts-in-spheres to output. [numpN=1-> P0, numpN=2->P0,P1, numpN=3->P0,P1,P2...\n");
    fprintf(stderr,"     * file         = name of data file\n") ;
    fprintf(stderr,"     * format       = format of data file  (a=ascii, c=csv, f=fast-food)\n") ;
    fprintf(stderr,"     * seed         = seed for random number generator\n");
    fprintf(stderr,"     * boxsize      = if periodic, the boxsize to use for the periodic wrap (0 means detect the particle extent)\n");
    fprintf(stderr,"     > VPFfile      = name of output file <r P0 P1 P2 ...>\n") ;
    fprintf(stderr,"\n\tCompile options: \n");
#ifdef PERIODIC
    fprintf(stderr,"Periodic = True\n");
#else
    fprintf(stderr,"Periodic = False\n");
#endif

#ifdef DOUBLE_PREC
    fprintf(stderr,"Precision = double\n");
#else
    fprintf(stderr,"Precision = float\n");
#endif

    fprintf(stderr,"=========================================================================\n") ;
}
