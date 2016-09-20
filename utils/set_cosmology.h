/* File: set_cosmology.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


#include "cosmology_params.h"

    double get_age(const double z);
    double agefunc(double z,void *params);
    double get_comoving_distance(const double zlow, const double z);
    double comoving_distance_func(const double z, void *params);
    double epeebles(const double z);


#ifdef __cplusplus
}
#endif
