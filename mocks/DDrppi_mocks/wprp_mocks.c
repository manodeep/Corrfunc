/* File: wprp_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* PROGRAM wprp

--- wprp nrbin ND1 ND2 NR1 NR2 D1D2 D1R2 D2R1 R1R2 [pimax] > wpfile
--- Measure the cross-correlation function wp(rp) for two different
   data and random files (or autocorrelation if D1=D2 and R1=R2).
   Uses the Szalay estimator: xi(rp,pi) = (D1D2 - D1R2 - D2R1 - R1R2)/RR2
  * nrpbin = number of bins (logarithmically spaced in r)
  * ND1    = number of data points in data1
  * ND2    = number of data points in data2
  * NR1    = number of random points in random1
  * NR2    = number of random points in random2
  * D1D2   = data1-data2 pairs (output of corrfunc.c)
  * D1R2   = data1-random2 pairs (output of corrfunc.c)
  * D2R1   = data2-random1 pairs (output of corrfunc.c)
  * R1R2   = random1-random2 pairs (output of corrfunc.c)
  *[pimax] = maximum pi in integral to compute wp(rp) (default = 40 Mpc/h)
  > wpfile = name of output file <rpmin rpmax rpavg logrp wp(rp) DDtot>
   ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#include "macros.h"
#include "utils.h"

void Printhelp(void);

int main(int argc, char *argv[])
{
    int i,j;
    /*---Arguments-------------------------*/
    int nrpbin,ND1,ND2,NR1,NR2 ;
    double pimax ;
    FILE *fpD1D2,*fpD1R2,*fpD2R1,*fpR1R2 ;
    /*---Data-variables--------------------*/
    int ndat ;
    double *logrp1,*ravg1,*pi1,*D1D2,*D1R2,*D2R1,*R1R2 ;
    /*---Corrfunc-variables----------------*/
    int npibin,indx,npimax ;
    double dpi ;
    double fN1,fN2,*xirppi,*wp,*logrp,*rpavg,*DDtot ;



    /*---Read-arguments-----------------------------------*/

    if(argc<10) {
        Printhelp() ;
        return(EXIT_FAILURE);
    }
    sscanf(argv[1],"%d",&nrpbin) ;
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

    pimax=40. ;
    if(argc>10) sscanf(argv[10],"%lf",&pimax) ;

    /*----------------------------------------------------*/

    npibin = 40 ;
    dpi = 1. ;
    npimax = (int)(pimax/dpi) ;
    if(npibin < npimax)
        npibin = npimax;


    /*---Read-D1D2-file-----------------------------------*/
    ndat=nrpbin*npibin ;
    logrp1 = my_calloc(sizeof(*logrp1),ndat) ;
    pi1    = my_calloc(sizeof(*pi1),ndat) ;
    ravg1  = my_calloc(sizeof(*ravg1),ndat);
    D1D2   = my_calloc(sizeof(*D1D2),ndat);

    i = 0 ;
    while(fscanf(fpD1D2,"%lf %lf %lf %lf%*[^\n]",&D1D2[i],&ravg1[i],&logrp1[i],&pi1[i])!=EOF) {
        /* fprintf(stderr,"%lf %lf %lf %lf\n",D1D2[i],ravg1[i],logrp1[i],pi1[i]); */
        i++ ;
        XASSERT(i <= ndat, "Current line parsed = %d must be less than number of data points expected = %d\n"
                "nrpbin = %d npibin = %d", i, ndat, nrpbin, npibin);
    }

    if(i!=ndat) {
        fprintf(stderr,"wprp> Warning: nrpbin*npibin= %d, but D1D2 file has %d lines\n",ndat,i) ;
        return(-1) ;
    }
    fclose(fpD1D2);

    /*---Read-D1R2-file-----------------------------------*/
    D1R2=my_calloc(sizeof(*D1R2),ndat) ;
    i = 0 ;

    while(fscanf(fpD1R2,"%lf%*[^\n]",&D1R2[i])!=EOF) {
        i++ ;
        XASSERT(i <= ndat, "Current line parsed = %d must be less than number of data points expected = %d\n"
                "nrpbin = %d npibin = %d", i, ndat, nrpbin, npibin);
    }
    if(i!=ndat) {
        fprintf(stderr,"wprp> Warning: nrpbin*npibin= %d, but D1R2 file has %d lines\n",ndat,i) ;
        return(-1) ;
    }
    fclose(fpD1R2);

    /*---Read-D2R1-file-----------------------------------*/

    D2R1=my_calloc(sizeof(*D2R1),ndat);
    i = 0 ;
    while(fscanf(fpD2R1,"%lf%*[^\n]",&D2R1[i])!=EOF) {
        i++ ;
        XASSERT(i <= ndat, "Current line parsed = %d must be less than number of data points expected = %d\n"
                "nrpbin = %d npibin = %d", i, ndat, nrpbin, npibin);
    }
    if(i!=ndat) {
        fprintf(stderr,"wprp> Warning: nrpbin*npibin= %d, but D2R1 file has %d lines\n",ndat,i) ;
        return(-1) ;
    }
    fclose(fpD2R1);

    /*---Read-D1R2-file-----------------------------------*/
    R1R2=my_calloc(sizeof(*R1R2),ndat) ;
    i = 0 ;
    while(fscanf(fpR1R2,"%lf%*[^\n]",&R1R2[i])!=EOF) {
        i++ ;
        XASSERT(i <= ndat, "Current line parsed = %d must be less than number of data points expected = %d\n"
                "nrpbin = %d npibin = %d", i, ndat, nrpbin, npibin);
    }
    if(i!=ndat) {
        fprintf(stderr,"wprp> Warning: nrpbin*npibin= %d, but R1R2 file has %d lines\n",ndat,i) ;
        return(-1) ;
    }
    fclose(fpR1R2);

    /*---Compute-xi(rp,pi)--------------------------------*/
    xirppi=my_calloc(sizeof(*xirppi), ndat);
    fN1 = (double)NR1/(double)ND1 ;
    fN2 = (double)NR2/(double)ND2 ;

    for(i=0;i<ndat;i++) {
        xirppi[i] = (fN1*fN2*D1D2[i] - fN1*D1R2[i] - fN2*D2R1[i] + R1R2[i])/(R1R2[i]) ;
    }

    /*---Compute-wp(rp)-----------------------------------*/
    logrp = my_calloc(sizeof(*logrp),nrpbin) ;
    rpavg = my_calloc(sizeof(*rpavg),nrpbin) ;
    DDtot = my_calloc(sizeof(*DDtot),nrpbin) ;
    wp    = my_calloc(sizeof(*wp),nrpbin) ;

    indx=0 ;
    for(i=0;i<nrpbin;i++) {
        logrp[i] = logrp1[indx] ;
        for(j=0;j<npibin;j++) {
            wp[i] += 2.*dpi*xirppi[indx] ;
            rpavg[i] += D1D2[indx]*ravg1[indx] ;
            DDtot[i] += D1D2[indx] ;
            indx++ ;
        }
        if(DDtot[i]>0)
            rpavg[i] /= DDtot[i] ;
        else
            rpavg[i] = 0.0;
    }

    /*---Print-wp(rp)-------------------------------------*/
    double rpmin=0.0;
    {
        double dlogrp = logrp[1]-logrp[0];
        rpmin = pow(10.0,logrp[0]-dlogrp);
    }

    for(i=0;i<nrpbin;i++) {
        double rpmax = pow(10.0,logrp[i]);
        fprintf(stdout,"%12.6lf %12.6lf %12.6lf %12.6lf %12"PRId64"\n",rpmin,rpmax,rpavg[i],wp[i],(int64_t) DDtot[i]) ;//rpavg actually contains the upper limit of the bin
        rpmin = rpmax;
    }

    free(logrp);
    free(rpavg);
    free(DDtot);
    free(wp);
    free(logrp1);free(pi1);free(D1D2);
    free(D1R2);
    free(D2R1);
    free(R1R2);
    free(xirppi);

    return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/

void Printhelp(void)
{
    fprintf(stderr,"=========================================================================\n") ;
    fprintf(stderr,"   --- wprp nrbin ND1 ND2 NR1 NR2 D1D2 D1R2 D2R1 R1R2 [pimax] > wpfile\n") ;
    fprintf(stderr,"   --- Measure the cross-correlation function wp(rp) for two different\n") ;
    fprintf(stderr,"       data and random files (or autocorrelation if D1=D2 and R1=R2).\n") ;
    fprintf(stderr,"       Uses the Szalay estimator: xi(rp,pi) = (D1D2 - D1R2 - D2R1 - R1R2)/RR2\n") ;
    fprintf(stderr,"      * nrpbin = number of bins (logarithmically spaced in r)\n") ;
    fprintf(stderr,"      * ND1    = number of data points in data1\n") ;
    fprintf(stderr,"      * ND2    = number of data points in data2\n") ;
    fprintf(stderr,"      * NR1    = number of random points in random1\n") ;
    fprintf(stderr,"      * NR2    = number of random points in random2\n") ;
    fprintf(stderr,"      * D1D2   = data1-data2 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      * D1R2   = data1-random2 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      * D2R1   = data2-random1 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      * R1R2   = random1-random2 pairs (output of corrfunc.c)\n") ;
    fprintf(stderr,"      *[pimax] = maximum pi in integral to compute wp(rp) (default = 40 Mpc/h)\n") ;
    fprintf(stderr,"      > wpfile = name of output file <rpmin rpmax rpavg logrp wp(rp) DDtot>\n") ;
    fprintf(stderr,"=========================================================================\n") ;
}
