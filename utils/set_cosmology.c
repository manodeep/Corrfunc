/* File: set_cosmology.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include "set_cosmology.h"
#include "cosmology_params.h"
#include "set_cosmo_dist.h"
#include "macros.h"

#include<math.h>
#include<gsl/gsl_integration.h>


double get_age(const double z)
{
    const int NWORKSPACE=1000;
    const double RECOMBINATION_REDSHIFT=1e3;
    const double AGE_AT_RECOMBINATION=0.37*1e-3;/*in Gyr ( recombination = 0.37 Myr)*/
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(NWORKSPACE);
    gsl_function F;
    double dummy=0.0;
    double result=0.0,error=0.0;

    F.function = &agefunc;
    F.params = &dummy;
    gsl_integration_qags (&F, z, RECOMBINATION_REDSHIFT, 0, 1e-7,NWORKSPACE,w,&result, &error);
    result *=  9.77813/LITTLE_H;

    result += AGE_AT_RECOMBINATION;

    gsl_integration_workspace_free (w);
    return result;
}

double agefunc(double z,void *params)
{
    (void) params;
    return 1.0/(epeebles(z)*(1.0+z));
}


double get_comoving_distance(const double zlow, const double z)
{
    const int NWORKSPACE=1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(NWORKSPACE);
    gsl_function F;
    double dummy=0.0;
    double result=0.0,error=0.0;

    if(comoving_distance_func(0.0, NULL) < 0) {
        return -1.0;
    }
    
    F.function = &comoving_distance_func;
    F.params = &dummy;
    gsl_integration_qag (&F, zlow, z, 0, 1e-7, GSL_INTEG_GAUSS51, NWORKSPACE,w,&result, &error);

    gsl_integration_workspace_free (w);
    const double smallh=1.0;
    const double Dh = SPEED_OF_LIGHT*0.01/smallh ;// c/(100) -> in units of little h^-1 Mpc
    return Dh*result;
}


double comoving_distance_func(const double z, void *params)
{
    (void) params;
    return 1.0/epeebles(z);
}



double epeebles(const double z)
{
    if(cosmology_initialized != 1) {
        XRETURN(cosmology_initialized == 1, -1.0, "Cosmology needs to be initialized before calling the co-moving distance routines.\n"
                "initialize cosmology by calling the function `init_cosmology' in cosmology_params.c\n");
        return -1.0;
    }
    const double Omegak = 1.0 - OMEGA_M - OMEGA_L;
    const double ez = sqrt(OMEGA_M*(1.0+z)*(1.0+z)*(1.0+z) + Omegak *(1+z) +  OMEGA_L);
    return ez;
}
