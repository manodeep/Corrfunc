#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

#ifdef __cplusplus
extern "C" {
#endif
	

#include "cosmology_params.h"

double get_age(const double z);
double agefunc(double z,void *params);
double get_comoving_distance(const double z);
double comoving_distance_func(const double z, void *params);
double epeebles(const double z);

	
#ifdef __cplusplus
}
#endif
		
