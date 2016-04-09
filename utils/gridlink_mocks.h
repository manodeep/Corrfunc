/* File: gridlink_mocks.h */
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

#include "cellarray_mocks.h"


    void get_max_min_data(const int64_t ND1, const DOUBLE * restrict cz,
                          DOUBLE *min_cz, DOUBLE *max_cz
#ifdef LINK_IN_DEC
                          ,const DOUBLE * restrict dec,
                          DOUBLE *min_dec, DOUBLE *max_dec
#endif
                          );



    //For DDrppi, when neither linking in dec or ra -> only linking in cz
    cellarray_mocks *gridlink1D(const int64_t np,const DOUBLE czmin,const DOUBLE czmax, const DOUBLE rcell,
                                const DOUBLE *theta, const DOUBLE *phi, const DOUBLE *cz,
                                int *ngrid,int *max_in_cell,
                                const int zbin_refine_factor) __attribute__((warn_unused_result));

    //For DDrppi, when linking in cz + dec
    cellarray_mocks **gridlink2D(const int64_t np,
                                 const DOUBLE czmin, const DOUBLE czmax, const DOUBLE rcell,
                                 const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE rpmax,
                                 const DOUBLE *cz,const DOUBLE *dec, const DOUBLE *ra,
                                 int *ngrid_cz,
                                 int **ngrid_declination,
                                 int *max_in_cell,
                                 const int rbin_refine_factor,
                                 const int zbin_refine_factor) __attribute__((warn_unused_result));


    //For DDrppi, when linking in cz, dec + ra
    cellarray_mocks *** gridlink3D(const int64_t np,
                                   const DOUBLE czmin,const DOUBLE czmax,const DOUBLE rcell,
                                   const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE rpmax,
                                   const DOUBLE * restrict cz,
                                   const DOUBLE * restrict dec,
                                   const DOUBLE * restrict phi,
                                   int *ngrid_cz,
                                   int **ngrid_declination,
                                   const DOUBLE phi_min,const DOUBLE phi_max,
                                   int ***ngrid_phi,
                                   int *max_in_cell,
                                   const int phibin_refine_factor,
                                   const int rbin_refine_factor,
                                   const int zbin_refine_factor) __attribute__((warn_unused_result));


    //For DDtheta, when linking in dec
    cellarray * gridlink1D_theta(const int64_t np,
                                 const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE thetamax,
                                 const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1, const DOUBLE  * restrict dec,
                                 int *ngrid_declination,
                                 int *max_in_cell,
                                 const int rbin_refine_factor) __attribute__((warn_unused_result));


    //For DDtheta, when linking in dec + ra
    cellarray **gridlink2D_theta(const int64_t np,
                                 const DOUBLE dec_min, const DOUBLE dec_max,const DOUBLE thetamax,
                                 const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1,
                                 const DOUBLE * restrict dec,
                                 int *ngrid_declination,
                                 const DOUBLE * restrict phi, const DOUBLE phi_min,const DOUBLE phi_max,
                                 int **ngrid_phi,
                                 int *max_in_cell,
                                 const int rbin_refine_factor,
                                 const int phibin_refine_factor) __attribute__((warn_unused_result));



    //Following definitions are directly from the theory side gridlink.h -> included to avoid namespace collisions (cellarray struct is defined differently in theory side)

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
    
    
#ifdef __cplusplus
}
#endif
