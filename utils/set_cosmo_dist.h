/* File: set_cosmo_dist.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#ifndef _SET_COSMO_DIST_H_
#define _SET_COSMO_DIST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_integration.h"

#ifdef __cplusplus
extern "C" {
#endif


#define COSMO_DIST_SIZE 10000
#define MAX_REDSHIFT_FOR_COSMO_DIST     1.0
#define SPEED_OF_LIGHT 299800.0



    int set_cosmo_dist(const double zmax,const int max_size,double *zc,double *dc,const int lasdamas_cosmology);

#ifdef __cplusplus
}
#endif


#endif
