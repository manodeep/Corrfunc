#ifndef _SET_COSMO_DIST_H_
#define _SET_COSMO_DIST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define COSMO_DIST_SIZE 10000
#define MAX_REDSHIFT_FOR_COSMO_DIST     1.0 //used to be 0.2 for Consuelo
#define SPEED_OF_LIGHT 299800.0

int set_cosmo_dist(const double zmax,const int max_size,double *zc,double *dc,const int lasdamas_cosmology);

#endif
