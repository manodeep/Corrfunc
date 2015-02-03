/* File: gridlink.h */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include "cellarray.h"


void get_max_min_data(const int64_t ND1, const DOUBLE * restrict cz, 
					  DOUBLE *min_cz, DOUBLE *max_cz
#ifdef LINK_IN_DEC
					  ,const DOUBLE * restrict dec, 
					  DOUBLE *min_dec, DOUBLE *max_dec
#endif

#ifdef LINK_IN_RA
					  ,const DOUBLE * restrict ra,
					  DOUBLE *min_ra, DOUBLE *max_ra
#endif
					  );






//For DDrppi, when neither linking in dec or ra -> only linking in cz
cellarray_mocks *gridlink1D(const int64_t np,const DOUBLE czmin,const DOUBLE czmax, const DOUBLE rcell,
						   const DOUBLE *theta, const DOUBLE *phi, const DOUBLE *cz, 
						   int *ngrid,int *max_in_cell,
						   const int zbin_refine_factor) __attribute__((warn_unused_result));

//For DDrppi, when linking in cz + dec
cellarray_mocks ** gridlink2D(const int64_t np,
						const DOUBLE czmin, const DOUBLE czmax, const DOUBLE rcell,
						const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE rpmax,
						const DOUBLE *x1,const DOUBLE *y1,const DOUBLE *z1, const DOUBLE *cz,const DOUBLE *dec,
						int *ngrid_cz,
						int **ngrid_declination,
						int *max_in_cell,
						const int rbin_refine_factor,
						const int zbin_refine_factor) __attribute__((warn_unused_result));
  

//For DDrppi, when linking in cz, dec + ra
cellarray_mocks *** gridlink3D(const int64_t np,
							   const DOUBLE czmin,const DOUBLE czmax,const DOUBLE rcell,
							   const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE rpmax,
							   const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1,
							   const DOUBLE * restrict cz,
							   const DOUBLE * restrict dec,
							   int *ngrid_cz,
							   int **ngrid_declination,
							   const DOUBLE * restrict phi, 
							   const DOUBLE phi_min,const DOUBLE phi_max,
							   int ***ngrid_phi,
							   int *max_in_cell,
							   const int zbin_refine_factor,
							   const int rbin_refine_factor,
							   const int phibin_refine_factor) __attribute__((warn_unused_result));


//For DDtheta, when linking in dec
  cellarray_mocks * gridlink1D_theta(const int64_t np,
							   const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE thetamax,
							   const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1, const DOUBLE  * restrict dec,
							   int *ngrid_declination,
							   int *max_in_cell,
							   const int rbin_refine_factor) __attribute__((warn_unused_result));

  
//For DDtheta, when linking in dec + ra
cellarray_mocks **gridlink2D_theta(const int64_t np,
							 const DOUBLE dec_min, const DOUBLE dec_max,const DOUBLE thetamax,
							 const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1,
							 const DOUBLE * restrict dec,
							 int *ngrid_declination,
							 const DOUBLE * restrict phi, const DOUBLE phi_min,const DOUBLE phi_max,
							 int **ngrid_phi,
							 int *max_in_cell,
							 const int rbin_refine_factor,
							 const int phibin_refine_factor) __attribute__((warn_unused_result));
  

	
#ifdef __cplusplus
}
#endif
