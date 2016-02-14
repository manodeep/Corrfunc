/* File: set_cosmology.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include "set_cosmology.h"

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
    result *=  9.77813/PARAMS.COSMO->h100;

    result += AGE_AT_RECOMBINATION;

    gsl_integration_workspace_free (w);
    return result;
}

double agefunc(double z,void *params)
{
    return 1.0/(epeebles(z)*(1.0+z));
}


double get_comoving_distance(const double z)
{
    const int NWORKSPACE=1000;
    const double RECOMBINATION_REDSHIFT=1e3;
    const double AGE_AT_RECOMBINATION=0.37*1e-3;/*in Gyr ( recombination = 0.37 Myr)*/
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(NWORKSPACE);
    gsl_function F;
    double dummy=0.0;
    double result=0.0,error=0.0;

    F.function = &comoving_distance_func;
    F.params = &dummy;
    gsl_integration_qags (&F, 0.0, z, 0, 1e-7,NWORKSPACE,w,&result, &error);
    /* result *=  9.77813/PARAMS.COSMO->h100; */
    /* result += AGE_AT_RECOMBINATION; */

    gsl_integration_workspace_free (w);
    return result;
}


double comoving_distance_func(const double z, void *params)
{
    return epeebles(z);
}



double epeebles(const double z)
{
    assert(cosmology_initialized==1 && "Cosmology needs to be initialized before calling the distance routines");
    double ez = sqrt(OMEGA_M*(1.0+z)*(1.0+z)*(1.0+z) +  OMEGA_L);/*assumes flat Universe with only matter and lambda*/
    return ez;
}
