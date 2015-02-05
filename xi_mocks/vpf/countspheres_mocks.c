#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "cellarray.h"
#include "function_precision.h"

int count_neighbors(const DOUBLE xcen,const DOUBLE ycen,const DOUBLE zcen,const DOUBLE smin,const DOUBLE inv_rcube,const DOUBLE rmax,
					const int ngrid,const cellarray *lattice, const int nthreshold, const int bin_refine_factor)
{
  int numngb=0;
  const cellarray *cellstruct;
  const DOUBLE rmax_sqr = (DOUBLE) (rmax*rmax);
  int ix = (int)(ngrid*(xcen-smin)*inv_rcube);
  int iy = (int)(ngrid*(ycen-smin)*inv_rcube);
  int iz = (int)(ngrid*(zcen-smin)*inv_rcube);
  if(ix > ngrid-1) ix--;
  if(iy > ngrid-1) iy--;
  if(iz > ngrid-1) iz--;
  
  assert(ix >= 0 && ix < ngrid && "x-position is inside limits");
  assert(iy >= 0 && iy < ngrid && "y-position is inside limits");
  assert(iz >= 0 && iz < ngrid && "z-position is inside limits");
  
  const int min_ix = ix - bin_refine_factor < 0 ?             0:ix - bin_refine_factor;
  const int max_ix = ix + bin_refine_factor > ngrid-1 ? ngrid-1:ix + bin_refine_factor;
  for(int iix=min_ix;iix<=max_ix;iix++) {
    const DOUBLE newxpos = xcen;
    const int min_iy = iy - bin_refine_factor < 0 ?             0:iy - bin_refine_factor;
    const int max_iy = iy + bin_refine_factor > ngrid-1 ? ngrid-1:iy + bin_refine_factor;
    
    for(int iiy=min_iy;iiy<=max_iy;iiy++) {
      const DOUBLE newypos = ycen;
      const int min_iz = iz - bin_refine_factor < 0 ?             0:iz - bin_refine_factor;
      const int max_iz = iz + bin_refine_factor > ngrid-1 ? ngrid-1:iz + bin_refine_factor;
      
      for(int iiz=min_iz;iiz<=max_iz;iiz++) {
		DOUBLE newzpos = zcen;
		const int64_t index=iix*ngrid*ngrid + iiy*ngrid + iiz;
		cellstruct = &(lattice[index]);
		const DOUBLE *x2 = cellstruct->x;
		const DOUBLE *y2 = cellstruct->y;
		const DOUBLE *z2 = cellstruct->z;
#if  __INTEL_COMPILER
#pragma simd reduction(+:numngb) vectorlengthfor(DOUBLE)
#endif
		for(int i=0;i<cellstruct->nelements;i++) {
		  const DOUBLE dx=x2[i]-newxpos;
		  const DOUBLE dy=y2[i]-newypos;
		  const DOUBLE dz=z2[i]-newzpos;
		  const DOUBLE r2 = dx*dx + dy*dy + dz*dz;
		  if (r2 < rmax_sqr) numngb++;
		}
		
		if(numngb > nthreshold)
		  return numngb;
      }
    }
  }
  return numngb;
}
