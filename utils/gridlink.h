/* File: gridlink.h */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h"
#include "cellarray.h"

cellarray * gridlink(const int64_t np,
										 const DOUBLE *x,const DOUBLE *y,const DOUBLE *z,
										 const DOUBLE xmin, const DOUBLE xmax,
										 const DOUBLE ymin, const DOUBLE ymax,
										 const DOUBLE zmin, const DOUBLE zmax,
										 const DOUBLE max_x_size,
										 const DOUBLE max_y_size,
										 const DOUBLE max_z_size,
										 const int xbin_refine_factor,
										 const int ybin_refine_factor,
										 const int zbin_refine_factor,
										 int *nlattice_x,
										 int *nlattice_y,
										 int *nlattice_z) __attribute__((warn_unused_result)); 
										 

cellarray_nvec * gridlink_nvec(const int64_t np,
															 const DOUBLE *x,const DOUBLE *y,const DOUBLE *z,
															 const DOUBLE xmin, const DOUBLE xmax,
															 const DOUBLE ymin, const DOUBLE ymax,
															 const DOUBLE zmin, const DOUBLE zmax,
															 const DOUBLE max_x_size,
															 const DOUBLE max_y_size,
															 const DOUBLE max_z_size,
															 const int xbin_refine_factor,
															 const int ybin_refine_factor,
															 const int zbin_refine_factor,
															 int *nlattice_x,
															 int *nlattice_y,
															 int *nlattice_z) __attribute__((warn_unused_result)); 

	
void get_max_min(const int64_t ND1, const DOUBLE * restrict X1, const DOUBLE * restrict Y1, const DOUBLE * restrict Z1,
				 DOUBLE *min_x, DOUBLE *min_y, DOUBLE *min_z, DOUBLE *max_x, DOUBLE *max_y, DOUBLE *max_z);

	
#ifdef __cplusplus
}
#endif
