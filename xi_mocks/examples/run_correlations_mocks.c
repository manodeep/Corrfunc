/* File: run_correlations_mocks.c */
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

#include "function_precision.h"
#include "io.h"
#include "defs.h"
#include "utils.h"
#include "cosmology_params.h"

/* Library proto-types + struct definitions in the ../..//include directory */
#include "countpairs_rp_pi_mocks.h"
#include "countpairs_theta_mocks.h"
#include "countspheres_mocks.h"

#ifndef MAXLEN
#define MAXLEN 500
#endif

void Printhelp(void);

void Printhelp(void)
{
    fprintf(stderr,ANSI_COLOR_RED "=========================================================================\n") ;
    fprintf(stderr,"   --- run_correlations_mocks file format binfile boxsize [Nthreads] \n") ;
    fprintf(stderr,"   --- Measure the auto-correlation functions DD(r), DD(rp, pi) and wp(rp) for a single file\n");
    fprintf(stderr,"     * file         = name of data file\n") ;
    fprintf(stderr,"     * format       = format of data file  (a=ascii, f=fast-food)\n") ;
    fprintf(stderr,"     * binfile      = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
    fprintf(stderr,"     * pimax        = pimax   (in same units as X/Y/Z of the data)\n");
    fprintf(stderr,"     * cosmology    = flag to pick-up the cosmology combination to use (set as an array of combinations in ../utils/cosmology_params.c)\n");
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"     * numthreads   = number of threads to use\n");
#endif
    fprintf(stderr,"=========================================================================" ANSI_COLOR_RESET "\n") ;
}

int main(int argc, char **argv)
{
    char file[MAXLEN];
    char fileformat[MAXLEN];
    char binfile[MAXLEN];
    DOUBLE *ra1=NULL,*dec1=NULL,*cz1=NULL;
    struct timeval t0,t1;
    DOUBLE pimax;
    int cosmology=1;

#if !(defined(USE_OMP) && defined(_OPENMP))
    const char argnames[][30]={"file","format","binfile","pimax","cosmology"};
#else
    int nthreads=4;//default to 4 threads
    const char argnames[][30]={"file","format","binfile","pimax","cosmology","Nthreads"};
#endif
    int nargs=sizeof(argnames)/(sizeof(char)*30);

    if(argc > 1) {
        //Command-line options were supplied - check that they are correct
        if(argc < (nargs + 1) ) {
            //Not enough options were supplied
            Printhelp();
            exit(EXIT_FAILURE);
        } else {
            //Correct number of options - let's parse them.
            my_snprintf(file,MAXLEN, "%s",argv[1]);
            my_snprintf(fileformat,MAXLEN, "%s",argv[2]);
            my_snprintf(binfile,MAXLEN,"%s",argv[3]);
            pimax=atof(argv[4]);
            cosmology=atoi(argv[5]);
#if defined(USE_OMP) && defined(_OPENMP)
            nthreads = atoi(argv[6]);
#endif
        }
    } else {
        my_snprintf(file, MAXLEN, "%s", "../tests/data/Mr19_mock_northonly.rdcz.dat");
        my_snprintf(fileformat, MAXLEN, "%s","a");
        my_snprintf(binfile, MAXLEN,"%s","../tests/bins");
        pimax=40.0;
        cosmology=1;
    }

    fprintf(stderr,ANSI_COLOR_BLUE  "Running `%s' with the parameters \n",argv[0]);
    fprintf(stderr,"\n\t\t -------------------------------------\n");
    fprintf(stderr,"\t\t %-10s = %s \n",argnames[0],file);
    fprintf(stderr,"\t\t %-10s = %s \n",argnames[1],fileformat);
    fprintf(stderr,"\t\t %-10s = %s \n",argnames[2],binfile);
    fprintf(stderr,"\t\t %-10s = %10.4lf\n",argnames[3],pimax);
    fprintf(stderr,"\t\t %-10s = %d\n",argnames[4],cosmology);
#if defined(USE_OMP) && defined(_OPENMP)
    fprintf(stderr,"\t\t %-10s = %d\n",argnames[5],nthreads);
#endif
    fprintf(stderr,"\t\t -------------------------------------" ANSI_COLOR_RESET "\n");

    init_cosmology(cosmology);


    //Read-in the data
    const int64_t ND1 = read_positions(file,fileformat,sizeof(DOUBLE),3, &ra1, &dec1, &cz1);

    int autocorr=1;
    DOUBLE *ra2 = ra1;
    DOUBLE *dec2 = dec1;
    DOUBLE *cz2 = cz1;
    int64_t ND2 = ND1;

    //Do the DD(rp, pi) counts
    {
        gettimeofday(&t0,NULL);
#if defined(USE_OMP) && defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(rp,pi) calculation would be:\n `%s %s %s %s %s %s %lf %d %d'" ANSI_COLOR_RESET "\n",
                "../DDrppi/DDrppi_mocks",file,fileformat,file,fileformat,binfile,pimax,cosmology,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent DD(rp,pi) calculation would be:\n `%s %s %s %s %s %s %lf %d'" ANSI_COLOR_RESET "\n",
                "../DDrppi/DDrppi_mocks",file,fileformat,file,fileformat,binfile,pimax,cosmology);
#endif

        results_countpairs_mocks *results  = countpairs_mocks(ND1,ra1,dec1,cz1,
                                                              ND2,ra2,dec2,cz2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                              nthreads,
#endif
                                                              autocorr,
                                                              binfile,
                                                              pimax,
                                                              cosmology);

        gettimeofday(&t1,NULL);
        double pair_time = ADD_DIFF_TIME(t0,t1);
#if 0
        const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
        const int npibin = results->npibin;
        for(int i=1;i<results->nbin;i++) {
            const double logrp = LOG10(results->rupp[i]);
            for(int j=0;j<npibin;j++) {
                int index = i*(npibin+1) + j;
                fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->rpavg[index],logrp,(j+1)*dpi);
            }
        }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done DD(rp,pi) auto-correlation. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", ND1, pair_time);


        //free the result structure
        free_results_mocks(&results);
    }



    //Do the w(theta) counts
    {
        gettimeofday(&t0,NULL);
#if defined(USE_OMP) && defined(_OPENMP)
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent w(theta) calculation would be:\n `%s %s %s %s %s %s %lf %d'" ANSI_COLOR_RESET "\n",
                "../wtheta/DDtheta_mocks",file,fileformat,file,fileformat,binfile,pimax,nthreads);
#else
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent w(theta) calculation would be:\n `%s %s %s %s %s %s %lf'" ANSI_COLOR_RESET "\n",
                "../wtheta/DDtheta_mocks",file,fileformat,file,fileformat,binfile,pimax);
#endif

        results_countpairs_theta *results = countpairs_theta_mocks(ND1,ra1,dec1,
                                                                   ND2,ra2,dec2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                                   nthreads,
#endif
                                                                   autocorr,
                                                                   binfile) ;

        gettimeofday(&t1,NULL);
        DOUBLE pair_time = ADD_DIFF_TIME(t0,t1);

#if 0
        /*---Output-Pairs-------------------------------------*/
        DOUBLE theta_low = results->theta_upp[0];
        for(int i=1;i<results->nbin;i++) {
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->theta_avg[i],theta_low,results->theta_upp[i]);
            theta_low=results->theta_upp[i];
        }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done wtheta. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds" ANSI_COLOR_RESET "\n", ND1, pair_time);

        //free the result structure
        free_results_countpairs_theta(&results);
    }


    //Do the VPF
    {
        gettimeofday(&t0,NULL);
        const double rmax=10.0;
        const int nbin=10;
        const int nc=10000;
        const int num_pN=6;
        const int64_t Nran=nc;//Need to set it to nc so that the loop runs
        DOUBLE *xran=NULL,*yran=NULL,*zran=NULL;
        const int threshold_neighbors=1;
        const char centers_file[]="../tests/data/Mr19_centers_xyz_forVPF_rmax_10Mpc.txt";
        fprintf(stderr,ANSI_COLOR_MAGENTA "Command-line for running equivalent w(theta) calculation would be:\n `%s %lf %d %d %d %lf %s %s %s %s %s %d'" ANSI_COLOR_RESET "\n",
                "../vpf/vpf_mocks",rmax,nbin,nc,num_pN,0.0,file,fileformat,"junk","junkformat",centers_file,cosmology);

        results_countspheres_mocks *results = countspheres_mocks(ND1, ra1, dec1, cz1,
                                                                 Nran, xran, yran, zran,
                                                                 threshold_neighbors,
                                                                 rmax, nbin, nc,
                                                                 num_pN,
                                                                 centers_file,
                                                                 cosmology);



        gettimeofday(&t1,NULL);
        double sphere_time = ADD_DIFF_TIME(t0,t1);

#if 0
        //Output the results
        const DOUBLE rstep = rmax/(DOUBLE)nbin ;
        for(int ibin=0;ibin<results->nbin;ibin++) {
            const double r=(ibin+1)*rstep;
            fprintf(stdout,"%10.2"DOUBLE_FORMAT" ", r);
            for(int i=0;i<num_pN;i++) {
                fprintf(stdout," %10.4e", (results->pN)[ibin][i]);
            }
            fprintf(stdout,"\n");
        }
#endif
        fprintf(stderr,ANSI_COLOR_GREEN "Done VPF. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds" ANSI_COLOR_RESET "\n", ND1, sphere_time);
        free_results_countspheres_mocks(&results);
    }



    free(ra1);free(dec1);free(cz1);
    return EXIT_SUCCESS;
}
