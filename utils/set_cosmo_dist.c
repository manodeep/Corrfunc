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


#include "set_cosmo_dist.h"
#include "cosmology_params.h"

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define epsilon 1e-10

int set_cosmo_dist(const double zmax,const int max_size,double *zc,double *dc,const int lasdamas_cosmology)
{
    /* First, initialize cosmology*/
    init_cosmology(lasdamas_cosmology);

    int i = 0;
    double Omegak;
    double dz,z,Deltaz,z2;
    double Eint,E0,E1,E2,Dh,Dc;
    double smallh = 1.0;//Andreas pointed out that I don't need the real value of LITTLE_H
    /* double one_plus_z,one_plus_z_minus_dz,one_plus_z_sqr,one_plus_z_minus_dz_sqr; */
    Omegak = 1.0 - OMEGA_M - OMEGA_L;
    Dh = SPEED_OF_LIGHT*0.01/smallh ;// c/(100) -> in units of little h^-1 Mpc

    dz = 0.000001 ;
    Deltaz = 0.0001 ;

    Eint=0. ;
    E2=1. ;
    z2=Deltaz ;
    for(z=2.*dz;z<zmax;z+=2.*dz) {
        E0 = E2 ;
        /* one_plus_z = 1+z; */
        /* one_plus_z_minus_dz = one_plus_z - dz; */
        /* one_plus_z_sqr = SQR(one_plus_z); */
        /* one_plus_z_minus_dz_sqr = SQR(one_plus_z_minus_dz); */
        /* E1 = 1./sqrt(OMEGA_M*one_plus_z_minus_dz_sqr*one_plus_z_minus_dz + Omegak*one_plus_z_minus_dz_sqr + OMEGA_L) ; */
        /* E2 = 1./sqrt(OMEGA_M*one_plus_z_sqr*one_plus_z + Omegak*one_plus_z_sqr + OMEGA_L) ; */

        E1 = 1.0/sqrt(OMEGA_M * CUBE(1+z-dz) + Omegak *(1+z-dz) + OMEGA_L);
        E2 = 1.0/sqrt(OMEGA_M * CUBE(1+z) + Omegak *(1+z) + OMEGA_L);
        Eint += dz*(E0 + 4.*E1 + E2)/3. ;


        if(z>(z2-epsilon) && z<(z2+epsilon)) {
            if( i >= max_size )
                break ;

            Dc = Eint*Dh ;
            zc[i] = z ;
            dc[i] = Dc ;
            z2 += Deltaz ;
            i++ ;
        }
    }

#ifdef DEBUG
    fprintf(stderr,"cosmodist> (Omega_m, Omega_L, Omega_k, h, zmax) = (%4.2f, %4.2f, %4.2f, %4.2f,%4.2lf)\n",OMEGA_M,OMEGA_L,Omegak,LITTLE_H,zmax) ;
    fprintf(stderr,"cosmodist> tabulated redshift: %g to %g  (distance: %g to %g Mpc)\n", zc[0], zc[i-1], dc[0], dc[i-1]) ;
#endif
    return i ;
}
