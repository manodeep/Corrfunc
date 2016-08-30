/* File: cosmology_params.h */
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

    extern int cosmology_initialized;
    extern double OMEGA_M;
    extern double OMEGA_B;
    extern double OMEGA_L;
    extern double HUBBLE;
    extern double LITTLE_H;
    extern double SIGMA_8;
    extern double NS;
    extern int active_cosmology;

    int init_cosmology(const int lasdamas_cosmology)__attribute__((warn_unused_result));

#ifdef __cplusplus
}
#endif
