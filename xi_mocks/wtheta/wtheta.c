/* File: wtheta.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM wtheta

   --- wtheta nthetabin ND1 ND2 NR1 NR2 D1D2 D1R2 D2R1 R1R2 > wfile
   --- Measure the angular cross-correlation function w(theta) for two different
   data and random files (or autocorrelation if D1=D2 and R1=R2).
   Uses the Szalay estimator: w(theta) = (D1D2 - D1R2 - D2R1 + R1R2)/RR2

   * nthetabin = number of bins (logarithmically spaced in theta)
   * ND1    = number of data points in data1
   * ND2    = number of data points in data2
   * NR1    = number of random points in random1
   * NR2    = number of random points in random2
   * D1D2   = data1-data2 pairs (output of corrfunc.c)
   * D1R2   = data1-random2 pairs (output of corrfunc.c)
   * D2R1   = data2-random1 pairs (output of corrfunc.c)
   * R1R2   = random1-random2 pairs (output of corrfunc.c)
   > wfile = name of output file <thetamin thetamax thetaavg w(theta) D1D2>
   ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>

#include "utils.h"

void Printhelp(void);

int main(int argc, char **argv)
{
    int i;
    /*---Arguments-------------------------*/
    int nthetabin,ND1,ND2,NR1,NR2 ;
    FILE *fpD1D2,*fpD1R2,*fpD2R1,*fpR1R2 ;
    /*---Data-variables--------------------*/
    double *thetaavg,*D1D2,*D1R2,*D2R1,*R1R2 ;
    /*---Corrfunc-variables----------------*/
    double fN1,fN2,*wtheta ;
    double *rmin,*rmax;

    /*---Read-arguments-----------------------------------*/
    if(argc<10)  {
        Printhelp() ;
        return EXIT_FAILURE;
    }

    sscanf(argv[1],"%d",&nthetabin) ;
    sscanf(argv[2],"%d",&ND1) ;
    sscanf(argv[3],"%d",&ND2) ;
    sscanf(argv[4],"%d",&NR1) ;
    sscanf(argv[5],"%d",&NR2) ;
    fpD1D2=my_fopen(argv[6],"r") ;
    fpD1R2=my_fopen(argv[7],"r") ;
    fpD2R1=my_fopen(argv[8],"r") ;
    fpR1R2=my_fopen(argv[9],"r") ;
    if(fpD1D2 == NULL || fpD1R2 == NULL || fpD2R1 == NULL || fpR1R2 == NULL) {
        return EXIT_FAILURE;
    }


    /*---Read-D1D2-file-----------------------------------*/
    rmin     = my_malloc(sizeof(*rmin), nthetabin);
    rmax     = my_malloc(sizeof(*rmax), nthetabin);
    thetaavg = my_malloc(sizeof(*thetaavg),nthetabin);
    D1D2     = my_malloc(sizeof(*D1D2),    nthetabin);

    i = 0 ;
    while(fscanf(fpD1D2,"%lf %lf %lf %lf%*[^\n]",&rmin[i],&rmax[i],&thetaavg[i],&D1D2[i])!=EOF) {
        i++;
        assert(i<= nthetabin && "Output contains too many bins");
    }
    if(i!=nthetabin) {
        fprintf(stderr,"wtheta> Warning: nthetabin= %d, but D1D2 file has %d lines\n",nthetabin,i) ;
        return EXIT_FAILURE;
    }

    /*---Read-D1R2-file-----------------------------------*/
    D1R2=my_calloc(sizeof(*D1R2),nthetabin);
    i = 0 ;
    while(fscanf(fpD1R2,"%*f %*f %*f %lf%*[^\n]",&D1R2[i])!=EOF) {
        i++ ;
        assert(i<= nthetabin && "Output contains too many bins");
    }

    if(i!=nthetabin) {
        fprintf(stderr,"wtheta> Warning: nthetabin= %d, but D1R2 file has %d lines\n",nthetabin,i) ;
        return EXIT_FAILURE;
    }

    /*---Read-D2R1-file-----------------------------------*/
    D2R1=my_calloc(sizeof(*D2R1),nthetabin);
    i = 0 ;
    while(fscanf(fpD2R1,"%*f %*f %*f %lf%*[^\n]",&D2R1[i])!=EOF) {
        i++ ;
        assert(i<= nthetabin && "Output contains too many bins");
    }
    if(i!=nthetabin) {
        fprintf(stderr,"wtheta> Warning: nthetabin= %d, but D2R1 file has %d lines\n",nthetabin,i) ;
        return EXIT_FAILURE;
    }

    /*---Read-R1R2-file-----------------------------------*/
    R1R2=my_calloc(sizeof(*R1R2),nthetabin);
    i = 0 ;
    while(fscanf(fpR1R2,"%*f %*f %*f %lf%*[^\n]",&R1R2[i])!=EOF) {
        i++ ;
        assert(i<= nthetabin && "Output contains too many bins");
    }
    if(i!=nthetabin) {
        fprintf(stderr,"wtheta> Warning: nthetabin= %d, but R1R2 file has %d lines\n",nthetabin,i) ;
        return EXIT_FAILURE;
    }

    /*---Compute-w(theta)---------------------------------*/
    wtheta=my_calloc(sizeof(*wtheta),nthetabin);

    fN1 = (double)NR1/(double)ND1 ;
    fN2 = (double)NR2/(double)ND2 ;

    for(i=0;i<nthetabin;i++) {
        wtheta[i] = (fN1*fN2*D1D2[i] - fN1*D1R2[i] - fN2*D2R1[i] + R1R2[i])/R1R2[i] ;
    }

    /*---Print-w(theta)-----------------------------------*/
    for(i=0;i<nthetabin;i++) {
        fprintf(stdout,"%9.5lf %9.5lf %9.5lf %10.4lf %12"PRId64"\n",rmin[i],rmax[i],thetaavg[i],wtheta[i],(int64_t) D1D2[i]) ;
    }

    free(wtheta);
    free(R1R2);free(D2R1);free(D1R2);free(D1D2);
    free(thetaavg);free(rmin);free(rmax);
    fclose(fpD1D2);
    fclose(fpD1R2);
    fclose(fpD2R1);
    fclose(fpR1R2);

    return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{
    fprintf(stderr,"=========================================================================\n") ;
    fprintf(stderr,"   --- wtheta nthetabin ND1 ND2 NR1 NR2 D1D2 D1R2 D2R1 R1R2 > wfile\n") ;
    fprintf(stderr,"   --- Measure the angular cross-correlation function w(theta) for two different\n") ;
    fprintf(stderr,"       data and random files (or autocorrelation if D1=D2 and R1=R2).\n") ;
    fprintf(stderr,"       Uses the Szalay estimator: w(theta) = (D1D2 - D1R2 - D2R1 - R1R2)/RR2\n") ;
    fprintf(stderr,"      * nthetabin = number of bins (logarithmically spaced in theta)\n") ;
    fprintf(stderr,"      * ND1    = number of data points in data1\n") ;
    fprintf(stderr,"      * ND2    = number of data points in data2\n") ;
    fprintf(stderr,"      * NR1    = number of random points in random1\n") ;
    fprintf(stderr,"      * NR2    = number of random points in random2\n") ;
    fprintf(stderr,"      * D1D2   = data1-data2 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      * D1R2   = data1-random2 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      * D2R1   = data2-random1 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      * R1R2   = random1-random2 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      > wfile = name of output file <thetamin thetamax thetaavg w(theta) D1D2>\n") ;
    fprintf(stderr,"=========================================================================\n") ;
}
