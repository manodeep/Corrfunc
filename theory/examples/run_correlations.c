/* File: run_correlations.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/*
  Example code to show how to use the correlation function libraries
  Author: Manodeep Sinha <manodeep@gmail.com>
  Date: At some point in early 2015.



*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <inttypes.h>

#include "function_precision.h"
#include "io.h"
#include "defs.h"
#include "utils.h"
#include "macros.h"

/* Library proto-types + struct definitions in the ../../include directory
   These are the pair-counters
*/
#include "countpairs.h"
#include "countpairs_rp_pi.h"
#include "countpairs_s_mu.h"
#include "countpairs_wp.h"
#include "countpairs_xi.h"

//for vpf.
#include "countspheres.h"

#ifndef MAXLEN
#define MAXLEN 500
#endif

void Printhelp(void);

void Printhelp(void)
{
    fprintf(stderr,ANSI_COLOR_RED "=========================================================================\n") ;
    fprintf(stderr,"   --- run_correlations file format binfile boxsize [Nthreads] \n") ;
    fprintf(stderr,"   --- Measure the auto-correlation functions DD(r), DD(rp, pi) and wp(rp) for a single file\n");
    fprintf(stderr,"     * file         = name of data file\n") ;
    fprintf(stderr,"     * format       = format of data file  (a=ascii, f=fast-food)\n") ;
    fprintf(stderr,"     * binfile      = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
    fprintf(stderr,"     * boxsize      = BoxSize (in same units as X/Y/Z of the data)\n");
    fprintf(stderr,"     * pimax        = pimax   (in same units as X/Y/Z of the data)\n");
    fprintf(stderr,"     * mu_max       = Max. value of the cosine of the angle to the LOS (must be within [0.0, 1.0])\n");
    fprintf(stderr,"     * nmu_bins     = Number of linear bins to create (the bins themselves range from [0.0, mu_max]\n");
#if defined(_OPENMP)
    fprintf(stderr,"     * numthreads   = number of threads to use\n");
#endif
    fprintf(stderr,"=========================================================================" ANSI_COLOR_RESET "\n") ;
}

int main(int argc, char **argv)
{
    char file[MAXLEN];
    char fileformat[MAXLEN];
    char binfile[MAXLEN];
    DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;
    double boxsize;
    struct timeval t0,t1;
    DOUBLE pimax;
    DOUBLE mu_max=1.0;
    int nmu_bins=10;
    int nthreads=4;//default to 4 threads

    struct config_options options = get_config_options();
    options.verbose = 1;
    options.periodic = 1;
    options.float_type = sizeof(DOUBLE);

#if defined(_OPENMP)
    const char argnames[][30]={"file","format","binfile","boxsize","pimax","mu_max","nmu_bins","Nthreads"};
#else
    const char argnames[][30]={"file","format","binfile","boxsize","pimax","mu_max","nmu_bins"};
#endif
    int nargs=sizeof(argnames)/(sizeof(char)*30);

    if(argc > 1) {
        //Command-line options were supplied - check that they are correct
        if(argc < (nargs + 1) ) {
            //Not enough options were supplied
            Printhelp();
            return EXIT_FAILURE;
        } else {
            //Correct number of options - let's parse them.
            my_snprintf(file,MAXLEN, "%s",argv[1]);
            my_snprintf(fileformat,MAXLEN, "%s",argv[2]);
            my_snprintf(binfile,MAXLEN,"%s",argv[3]);
            boxsize=atof(argv[4]);
            pimax=atof(argv[5]);
            mu_max=atof(argv[6]);
            nmu_bins=atoi(argv[7]);
#if defined(_OPENMP)
            nthreads = atoi(argv[8]);
#endif
        }
    } else {
        my_snprintf(file, MAXLEN, "%s", "../tests/data/gals_Mr19.ff");
        my_snprintf(fileformat, MAXLEN, "%s","f");
        my_snprintf(binfile, MAXLEN,"%s","../tests/bins");
        boxsize=420.0;
        pimax=40.0;
        mu_max=1.0;
        nmu_bins=10;
    }

    fprintf(stderr,ANSI_COLOR_BLUE  "Running `%s' with the parameters \n",argv[0]);
    fprintf(stderr,"\n\t\t -------------------------------------\n");
    fprintf(stderr,"\t\t %-10s = %s \n",argnames[0],file);
    fprintf(stderr,"\t\t %-10s = %s \n",argnames[1],fileformat);
    fprintf(stderr,"\t\t %-10s = %s \n",argnames[2],binfile);
    fprintf(stderr,"\t\t %-10s = %10.4lf\n",argnames[3],boxsize);
    fprintf(stderr,"\t\t %-10s = %10.4lf\n",argnames[4],pimax);
    fprintf(stderr,"\t\t %-10s = %10.4lf\n",argnames[5],mu_max);
    fprintf(stderr,"\t\t %-10s = %d\n",argnames[6],nmu_bins);
#if defined(_OPENMP)
    fprintf(stderr,"\t\t %-10s = %d\n",argnames[7],nthreads);
#endif
    fprintf(stderr,"\t\t -------------------------------------" ANSI_COLOR_RESET "\n");

    //Read-in the data
    const int64_t ND1 = read_positions(file,fileformat,sizeof(*x1),3, &x1, &y1, &z1);

    int autocorr=1;
    DOUBLE *x2 = x1;
    DOUBLE *y2 = y1;
    DOUBLE *z2 = z1;
    int64_t ND2 = ND1;

    //Do the straight-up DD counts
    {
        gettimeofday(&t0,NULL);
#if defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(r) calculation would be:\n `%s %s %s %s %s %s %d'" ANSI_COLOR_RESET "\n",
                "../DD/DD",file,fileformat,file,fileformat,binfile,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(r) calculation would be:\n `%s %s %s %s %s %s'" ANSI_COLOR_RESET "\n",
                "../DD/DD",file,fileformat,file,fileformat,binfile);
#endif

        results_countpairs results;
        int status = countpairs(ND1,x1,y1,z1,
                                ND2,x2,y2,z2,
                                nthreads,
                                autocorr,
                                binfile,
                                &results,
                                &options, NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        gettimeofday(&t1,NULL);
        double pair_time = ADD_DIFF_TIME(t0,t1);
#if 0
        DOUBLE rlow=results.rupp[0];
        for(int i=1;i<results.nbin;i++) {
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results.npairs[i],results.rpavg[i],rlow,results.rupp[i]);
            rlow=results.rupp[i];
        }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done 3-d auto-correlation. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", ND1, pair_time);
        //The results structure contains the pair-counts


        //free the result structure
        free_results(&results);
    }

    //Do the DD(rp, pi) counts
    {
        gettimeofday(&t0,NULL);
#if defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(rp,pi) calculation would be:\n `%s %s %s %s %s %s %lf %d'" ANSI_COLOR_RESET "\n",
                "../DDrppi/DDrppi",file,fileformat,file,fileformat,binfile,pimax,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(rp,pi) calculation would be:\n `%s %s %s %s %s %s %lf'" ANSI_COLOR_RESET "\n",
                "../DDrppi/DDrppi",file,fileformat,file,fileformat,binfile,pimax);
#endif

        results_countpairs_rp_pi results;
        int status = countpairs_rp_pi(ND1,x1,y1,z1,
                                      ND2,x2,y2,z2,
                                      nthreads,
                                      autocorr,
                                      binfile,
                                      pimax,
                                      &results,
                                      &options, NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        gettimeofday(&t1,NULL);
        double pair_time = ADD_DIFF_TIME(t0,t1);
#if 0
        const int npibin = results.npibin;
        const DOUBLE dpi = pimax/(DOUBLE)results.npibin ;
        for(int i=1;i<results.nbin;i++) {
            const double logrp = LOG10(results.rupp[i]);
            for(int j=0;j<npibin;j++) {
                int index = i*(npibin+1) + j;
                fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi);
            }
        }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done DD(rp,pi) auto-correlation. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", ND1, pair_time);


        //free the result structure
        free_results_rp_pi(&results);
    }


    //Do the DD(s, mu) counts
    {
        gettimeofday(&t0,NULL);
#if defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(s,mu) calculation would be:\n `%s %s %s %s %s %s %lf %d %d'"ANSI_COLOR_RESET"\n",
                "../DDsmu/DDsmu",file,fileformat,file,fileformat,binfile,mu_max,nmu_bins,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(s,mu) calculation would be:\n `%s %s %s %s %s %s %lf %d'"ANSI_COLOR_RESET"\n",
                "../DDsmu/DDsmu",file,fileformat,file,fileformat,binfile,mu_max,nmu_bins);
#endif

        results_countpairs_s_mu results;
        int status = countpairs_s_mu(ND1,x1,y1,z1,
                                     ND2,x2,y2,z2,
                                     nthreads,
                                     autocorr,
                                     binfile,
                                     mu_max,
                                     nmu_bins,
                                     &results,
                                     &options, NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        gettimeofday(&t1,NULL);
        double pair_time = ADD_DIFF_TIME(t0,t1);
#if 0
    double smin = results.supp[0];
    const double dmu = mu_max/(double) nmu_bins;
    for(int i=1;i<results.nsbin;i++) {
        const double smax = results.supp[i];
        for(int j=0;j<nmu_bins;j++) {
            int index = i*(nmu_bins+1) + j;
            fprintf(stdout,"%e\t%e\t%e\t%12"PRIu64"\t%e\n", smin, smax, (j+1)*dmu, results.npairs[index], results.weightavg[index]);
        }
        smin = smax;
    }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done DD(s,mu) auto-correlation. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds "ANSI_COLOR_RESET"\n", ND1, pair_time);

        //free the result structure
        free_results_s_mu(&results);
    }



    //Do the wp counts
    {
        gettimeofday(&t0,NULL);
#if defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent wp calculation would be:\n `%s %lf %s %s %s %lf %d'" ANSI_COLOR_RESET "\n",
                "../wp/wp",boxsize,file,fileformat,binfile,pimax,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent wp calculation would be:\n `%s %lf %s %s %s %lf'" ANSI_COLOR_RESET "\n",
                "../wp/wp",boxsize,file,fileformat,binfile,pimax);
#endif
        results_countpairs_wp results;
        options.need_avg_sep = 1;
        options.verbose = 1;
        options.float_type = sizeof(*x1);
        int status = countpairs_wp(ND1,x1,y1,z1,
                                   boxsize,
                                   nthreads,
                                   binfile,
                                   pimax,
                                   &results,
                                   &options, NULL);
        gettimeofday(&t1,NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        double pair_time = ADD_DIFF_TIME(t0,t1);
#if 0
        double rlow=results.rupp[0];
        for(int i=1;i<results.nbin;++i) {
            fprintf(stdout,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",results.wp[i],results.rpavg[i],rlow,results.rupp[i],results.npairs[i]);
            rlow=results.rupp[i];
        }
#endif

        fprintf(stderr,ANSI_COLOR_GREEN "Done wp. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds" ANSI_COLOR_RESET "\n", ND1, pair_time);

        //free the result structure
        free_results_wp(&results);
    }


    //Do xi on the periodic cube
    {
        gettimeofday(&t0,NULL);
#if defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent xi calculation would be:\n `%s %lf %s %s %s %d'" ANSI_COLOR_RESET "\n",
                "../xi/xi",boxsize,file,fileformat,binfile,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent xi calculation would be:\n `%s %lf %s %s %s'" ANSI_COLOR_RESET "\n",
                "../xi/xi",boxsize,file,fileformat,binfile);
#endif
        results_countpairs_xi results;
        int status = countpairs_xi(ND1,x1,y1,z1,
                                   boxsize,
                                   nthreads,
                                   binfile,
                                   &results,
                                   &options, NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }
        gettimeofday(&t1,NULL);
        double pair_time = ADD_DIFF_TIME(t0,t1);
#if 0
        DOUBLE rlow=results.rupp[0];
        for(int i=1;i<results.nbin;++i) {
            fprintf(stdout,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",results.xi[i],results.rpavg[i],rlow,results.rupp[i],results.npairs[i]);
            rlow=results.rupp[i];
        }
#endif

        fprintf(stderr,ANSI_COLOR_GREEN "Done xi. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds" ANSI_COLOR_RESET "\n", ND1, pair_time);

        //free the result structure
        free_results_xi(&results);
    }


    //Now run the vpf -- not OpenMP'ized -- so will only use 1 core.
    {
        gettimeofday(&t0,NULL);
        const double rmax=10.0;
        const int nbin=10;
        const int nc=10000;
        const int num_pN=6;
        unsigned long seed=-1;

        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent vpf calculation would be:\n `%s %lf %d %d %d %s %s %ld'" ANSI_COLOR_RESET "\n",
                "../vpf/vpf",rmax,nbin,nc,num_pN,file,fileformat,seed);

        results_countspheres results;
        int status = countspheres(ND1, x1, y1, z1,
                                  rmax, nbin, nc,
                                  num_pN,
                                  seed,
                                  &results,
                                  &options, NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        gettimeofday(&t1,NULL);
        double sphere_time = ADD_DIFF_TIME(t0,t1);

#if 0
        //Output the results
        const DOUBLE rstep = rmax/(DOUBLE)nbin ;
        for(int ibin=0;ibin<results.nbin;ibin++) {
            const double r=(ibin+1)*rstep;
            fprintf(stdout,"%"REAL_FORMAT" ", r);
            for(int i=0;i<num_pN;i++) {
                fprintf(stdout," %10.4e", (results.pN)[ibin][i]);
            }
            fprintf(stdout,"\n");
        }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done VPF. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds" ANSI_COLOR_RESET "\n", ND1, sphere_time);
        free_results_countspheres(&results);
    }

    free(x1);free(y1);free(z1);
    return EXIT_SUCCESS;
}
