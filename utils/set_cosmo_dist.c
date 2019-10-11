/* File: set_cosmo_dist.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/* function set_cosmo_dist (modelled after cosmodist.c)

   --- Computes the comoving distance as a function of redshift for
   a given cosmological model.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "set_cosmo_dist.h"
#include "cosmology_params.h"

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define epsilon 1e-10

int set_cosmo_dist(const double zmax,const int max_size,double *zc,double *dc,const int lasdamas_cosmology)
{
    /* First, initialize cosmology*/
    int status = init_cosmology(lasdamas_cosmology);
    if(status != EXIT_SUCCESS) {
        return -1;
    }

    int i = 0;
    double Omegak;
    double smallh = 1.0;//Andreas pointed out that I don't need the real value of LITTLE_H
    Omegak = 1.0 - OMEGA_M - OMEGA_L;
    const double Dh = SPEED_OF_LIGHT*0.01/smallh ;// c/(100) -> in units of little h^-1 Mpc

    const double Deltaz = 1.0/max_size;
    const double dz = 1e-2*Deltaz;

    double Eint = 0.0;
    double E2 = 1.0 ;
    double z2 = Deltaz ;
    for(double z=2.0*dz;z<zmax;z+=2.0*dz) {
        double E0 = E2 ;

        double E1 = 1.0/sqrt(OMEGA_M * CUBE(1+z-dz) + Omegak *(1+z-dz) + OMEGA_L);
        E2 = 1.0/sqrt(OMEGA_M * CUBE(1+z) + Omegak *(1+z) + OMEGA_L);
        Eint += dz*(E0 + 4.*E1 + E2)/3. ;

        if(z>(z2-epsilon) && z<(z2+epsilon)) {
            if( i >= max_size )
                break ;

            const double Dc = Eint*Dh ;
            zc[i] = z ;
            dc[i] = Dc ;
            z2 += Deltaz ;
            i++ ;
        }
    }

#ifdef DEBUG
    fprintf(stderr,"%s> (Omega_m, Omega_L, Omega_k, h, zmax) = (%4.2f, %4.2f, %4.2f, %4.2f,%4.2lf)\n", __FUNCTION__, OMEGA_M,OMEGA_L,Omegak,LITTLE_H,zmax) ;
    fprintf(stderr,"%s> tabulated redshift: %g to %g  (distance: %g to %g Mpc)\n", __FUNCTION__, zc[0], zc[i-1], dc[0], dc[i-1]) ;
#endif
    return i ;
}
