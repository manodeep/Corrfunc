/* File: vpf_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM vpf_mocks

   --- vpf_mocks rmax nbin nc num_pN volume mocksfile mocksformat RANDfile randformat centers-file cosmology > output
   --- compute the void probability function for MOCK galaxies

   * rmax = maximum radius (in h^-1 Mpc)
   * nbin = number of radii (evenly spaced in r)
   * nc = number of centers to place (does not count rejected centers)
   * numpN              = number of counts-in-spheres to output. [numpN=1-> P0, numpN=2->P0,P1, numpN=3->P0,P1,P2...\n");
   * volume             = volume of sample (in Mpc^3/h^3)
   * galaxy file        = contains: ra,dec,cz)
   * galaxy file format = (ascii, fast-food)
   * random file        = (contains: ra,dec,cz)
   * random file format = (ascii, fast-food)
   * centers file       = file containing sphere centers (XYZ). If there are not enough centers, file will be truncated and re-written with centers
   * cosmology          = flag to pick-up the cosmology combination to use (set as an array of combinations in ../utils/cosmology_params.c)
   > output: <R P0 P1 P2 ...>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "defs.h"
#include "function_precision.h"
#include "utils.h"
#include "io.h"
#include "countspheres_mocks.h"


void Printhelp(void) ;

int main(int argc, char *argv[])
{
    /*---Arguments-------------------------*/
    int nbin,nc,num_pN;
    DOUBLE volume,rmax ;

    /*---Particle-distribution-variables---*/
    int64_t Ngal,Nran=0;

    DOUBLE *ra=NULL,*dec=NULL,*cz=NULL;
    DOUBLE *xran=NULL,*yran=NULL,*zran=NULL;
    char *galaxy_file,*galaxy_file_format,*random_file,*random_file_format,*centers_file;
    int cosmology=1;

    struct timeval tstart,t0,t1;

    const char argnames[][100]={"rmax","nbin","ncenters","num_pN","volume","galaxy file","galaxy file-format","randoms file","randoms file-format","centers file","cosmology flag"};
    int nargs=sizeof(argnames)/(sizeof(char)*100);

    gettimeofday(&tstart,NULL);
    /*---Read-arguments-----------------------------------*/
    if(argc < (nargs+1)) {
        Printhelp() ;
        fprintf(stderr,"\nFound: %d parameters\n ",argc-1);
        int i;
        for(i=1;i<argc;i++) {
            if(i <= nargs)
                fprintf(stderr,"\t\t %25s = `%s' \n",argnames[i-1],argv[i]);
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
    sscanf(argv[1],"%"REAL_FORMAT"",&rmax);
    sscanf(argv[2],"%d",&nbin);
    sscanf(argv[3],"%d",&nc);
    sscanf(argv[4],"%d",&num_pN);
    sscanf(argv[5],"%"REAL_FORMAT,&volume);
    galaxy_file=argv[6];
    galaxy_file_format=argv[7];
    random_file=argv[8];
    random_file_format=argv[9];
    centers_file=argv[10];
    cosmology = atoi(argv[11]);

    fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
    fprintf(stderr,"\n\t\t -------------------------------------\n");
    for(int i=1;i<argc;i++) {
        if(i <= nargs) {
            fprintf(stderr,"\t\t %-25s = %s \n",argnames[i-1],argv[i]);
        }  else {
            fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
        }
    }
    fprintf(stderr,"\t\t -------------------------------------\n");
    XRETURN(rmax > 0, EXIT_FAILURE, "rmax=%lf must be > 0\n",rmax);
    XRETURN(nbin > 0, EXIT_FAILURE,"Number of bins=%d must be > 0\n",nbin);
    XRETURN(nc > 0,EXIT_FAILURE,"Number of spheres=%d must be > 0\n",nc);


    int need_randoms = 1;
    int64_t num_centers_in_file=0;
    FILE *fpcen = fopen(centers_file,"r");
    if(fpcen != NULL) {
        double rr = 0.0;
        int num_read = fscanf(fpcen,"%*f %*f %*f %lf",&rr);
        assert(num_read == 1 && "Could not read max. sphere radius from the centers file");
        num_centers_in_file = getnumlines(centers_file,'#');
        if( rr >= rmax && num_centers_in_file >= nc) {
            need_randoms = 0;
            fclose(fpcen);
        } else {
            fclose(fpcen);
            num_centers_in_file = 0;
            need_randoms = 1;
        }
    } else {
        num_centers_in_file = 0;
        need_randoms = 1;
    }
    fprintf(stderr,"vpf_mocks> found %"PRId64" centers (need %d centers) - need randoms = %d\n",num_centers_in_file,nc,need_randoms);

    gettimeofday(&t0,NULL);
    /*---Read-galaxy-data1-file----------------------------------*/
    Ngal=read_positions(galaxy_file,galaxy_file_format, sizeof(DOUBLE), 3, &ra, &dec, &cz);
    gettimeofday(&t1,NULL);
    fprintf(stderr,"vpf_mocks> Ngal = %"PRId64". Time to read-in galaxies=%6.2lf sec\n",Ngal,ADD_DIFF_TIME(t0,t1)) ;

    /*---Read-random-file---------------------------------*/
    if(need_randoms == 1) {
        gettimeofday(&t0,NULL);
        Nran = read_positions(random_file,random_file_format, sizeof(DOUBLE), 3, &xran, &yran, &zran);
        gettimeofday(&t1,NULL);
        fprintf(stderr,"vpf_mocks> Nrandoms = %"PRId64". Time to read-in randoms = %6.2lf sec\n",Nran,ADD_DIFF_TIME(t0,t1)) ;
    }

    /*---Expected-number-of-randoms-in-sphere-------------*/
    int threshold_neighbors;
    if(need_randoms == 1) {
        assert(volume > 0 && "Mock volume must be > 0");
        double Nexpected = (double)Nran*(4.0*M_PI*(rmax*rmax*rmax)/3.)/volume ;
        fprintf(stderr,"vpf_mocks> Expected number of randoms in sphere = %lf\n",Nexpected) ;
        threshold_neighbors = (int) (Nexpected - sqrt(Nexpected));//allow 1-sigma deviation. Conservative
    } else {
        threshold_neighbors = 1;//dummy value -> just to ensure that the check does not compare with uninitialized values
        Nran = nc;// HACK: set Nran to number of spheres requested. Code will not execute loop otherwise
    }

    results_countspheres_mocks results;
    struct config_options options = get_config_options();

    int status = countspheres_mocks(Ngal, ra, dec, cz,
                                    Nran, xran, yran, zran,
                                    threshold_neighbors,
                                    rmax, nbin, nc,
                                    num_pN,
                                    centers_file,
                                    cosmology,
                                    &results,
                                    &options, NULL);

    free(ra);free(dec);free(cz);
    if(need_randoms == 1) {
        free(xran);free(yran);free(zran);
    }
    if(status != EXIT_SUCCESS){
        return status;
    }

    //Output the results
    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    for(int ibin=0;ibin<results.nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        fprintf(stdout,"%10.2"REAL_FORMAT" ", r);
        for(int i=0;i<num_pN;i++) {
            fprintf(stdout," %10.4e", (results.pN)[ibin][i]);
        }
        fprintf(stdout,"\n");
    }
    gettimeofday(&t1,NULL);
    fprintf(stderr,"vpf_mocks> Done. Ngal = %"PRId64". Time taken = %6.2lf sec\n",Ngal,ADD_DIFF_TIME(tstart,t1));

    free_results_countspheres_mocks(&results);
    return EXIT_SUCCESS ;
}

/*---Print-help-information---------------------------*/

void Printhelp(void)
{
    fprintf(stderr,"=========================================================================\n") ;
    fprintf(stderr,"   --- vpf_mocks rmax nbin nc num_pN volume mocksfile mocksformat RANDfile randformat centers-file cosmology > output\n");
    fprintf(stderr,"   --- compute the void probability function for SDSS galaxies\n") ;
    fprintf(stderr,"      * rmax = maximum radius (in h^-1 Mpc)\n") ;
    fprintf(stderr,"      * nbin = number of radii (evenly spaced in r)\n") ;
    fprintf(stderr,"      * nc = number of centers to place (does not count rejected centers)\n") ;
    fprintf(stderr,"      * numpN        = number of counts-in-spheres to output. [numpN=1-> P0, numpN=2->P0,P1, numpN=3->P0,P1,P2...\n");
    fprintf(stderr,"      * volume = volume of sample (in Mpc^3/h^3)\n") ;
    fprintf(stderr,"      * galaxy file, (contains: ra,dec,cz)\n") ;
    fprintf(stderr,"      * galaxy file format (a -> ascii, f->fast-food)\n") ;
    fprintf(stderr,"      * random file, (contains: ra,dec,cz)\n") ;
    fprintf(stderr,"      * random file format (a -> ascii, f-> fast-food)\n");
    fprintf(stderr,"      * file with sphere centers (centers will be read-in if enough centers exist, otherwise centers will be output into this file)\n");
    fprintf(stderr,"      > output: <R P0 P1 P2 ...>\n") ;

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

    fprintf(stderr,"=========================================================================\n") ;
}
