/* File: cellarray.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" //definition of DOUBLE depending on Makefile flag -DOUBLE_PREC

#define NLATMAX   100      /* maximum grid dimension in X-Y plane */

typedef struct{
    DOUBLE *x;
    DOUBLE *y;
    DOUBLE *z;
    int64_t nelements;//Here the xyz positions will be stored in their individual pointers. More amenable to sorting -> used by wp and xi
} cellarray;

typedef struct{
    DOUBLE *pos;
    int64_t nelements;
} cellarray_nvec;//Here the xyz positions will be stored as pos[x[NVEC],y{NVEC],z[NVEC],x[NVEC]...]. Note amenable to easy sorting -> used by xi_of_r and vpf


typedef struct cellarray_index cellarray_index;
    
/* This cellarray will avoid duplicating the particle positions */
struct cellarray_index{
    int64_t start;
    int64_t nelements;
    cellarray_index **ngb_cells;
    DOUBLE *xwrap;
    DOUBLE *ywrap;
    DOUBLE *zwrap;
    int64_t num_ngb;
};

#ifdef __cplusplus
}
#endif
