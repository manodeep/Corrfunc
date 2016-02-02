/* File: gridlink.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "gridlink.h"
#include "utils.h"
#include "function_precision.h"
#include "cellarray.h"

#define MEMORY_INCREASE_FAC   1.2

double get_binsize(const double xmin,const double xmax, const double rmax, const int refine_factor, const int max_ncells, int *nlattice)  __attribute__((warn_unused_result));

double get_binsize(const double xmin,const double xmax, const double rmax, const int refine_factor, const int max_ncells, int *nlattice)
{
  double xdiff = xmax-xmin;
  int nmesh=(int)(refine_factor*xdiff/rmax) ;
#ifdef PERIODIC
  if (nmesh<(2*refine_factor+1))  {
    fprintf(stderr,"linklist> ERROR:  nlattice = %d is so small that with periodic wrapping the same cells will be counted twice ....exiting\n",nmesh) ;
    exit(EXIT_FAILURE) ;
  }
#endif
  
  if (nmesh>max_ncells)  nmesh=max_ncells;
  double xbinsize = xdiff/nmesh;
  *nlattice = nmesh;
  return xbinsize;
}

void get_max_min(const int64_t ND1, const DOUBLE * restrict X1, const DOUBLE * restrict Y1, const DOUBLE * restrict Z1,
				 DOUBLE *min_x, DOUBLE *min_y, DOUBLE *min_z, DOUBLE *max_x, DOUBLE *max_y, DOUBLE *max_z)
{
  DOUBLE xmin = *min_x, ymin = *min_y, zmin=*min_z;
  DOUBLE xmax = *max_x, ymax = *max_y, zmax=*max_z;
	
  for(int64_t i=0;i<ND1;i++) {
    if(X1[i] < xmin) xmin=X1[i];
    if(Y1[i] < ymin) ymin=Y1[i];
    if(Z1[i] < zmin) zmin=Z1[i];


    if(X1[i] > xmax) xmax=X1[i];
    if(Y1[i] > ymax) ymax=Y1[i];
    if(Z1[i] > zmax) zmax=Z1[i];
  }
	*min_x=xmin;*min_y=ymin;*min_z=zmin;
	*max_x=xmax;*max_y=ymax;*max_z=zmax;
}


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
										 int *nlattice_z)
{
  cellarray *lattice=NULL;
  int ix,iy,iz;
  int nmesh_x,nmesh_y,nmesh_z;
  int64_t *nallocated=NULL;
  DOUBLE xdiff,ydiff,zdiff;
  DOUBLE cell_volume,box_volume;
  DOUBLE xbinsize,ybinsize,zbinsize;
  int64_t expected_n=0;
  int64_t totncells;

#ifndef SILENT	
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
#endif
	
  xbinsize = get_binsize(xmin,xmax,max_x_size,xbin_refine_factor, NLATMAX, &nmesh_x);
  ybinsize = get_binsize(ymin,ymax,max_y_size,ybin_refine_factor, NLATMAX, &nmesh_y);
  zbinsize = get_binsize(zmin,zmax,max_z_size,zbin_refine_factor, NLATMAX, &nmesh_z);
  
  totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;

  xdiff = xmax-xmin;
  ydiff = ymax-ymin;
  zdiff = zmax-zmin;
  
  cell_volume=xbinsize*ybinsize*zbinsize;
  box_volume=xdiff*ydiff*zdiff;
  expected_n=(int64_t)(np*cell_volume/box_volume*MEMORY_INCREASE_FAC);
  expected_n=expected_n < NVEC ? NVEC:expected_n;
	while((expected_n % NVEC) != 0)
		expected_n++;
	
#ifndef SILENT	
  fprintf(stderr,"In %s> Running with [nmesh_x, nmesh_y, nmesh_z]  = %d,%d,%d. ",__FUNCTION__,nmesh_x,nmesh_y,nmesh_z);
#endif	
  lattice    = (cellarray *) my_malloc(sizeof(cellarray), totncells);
  nallocated = (int64_t *)       my_malloc(sizeof(*nallocated)      , totncells);

  /*
	Allocate memory for each of the fields in cellarray. Since we haven't processed the data yet, 
	expected_n is a reasonable guess as to the number of points in the cell. 
   */
  for (int64_t index=0;index<totncells;index++) {
		const size_t memsize=sizeof(DOUBLE);
		lattice[index].x = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
		lattice[index].y = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
		lattice[index].z = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
		lattice[index].nelements=0;
		nallocated[index] = expected_n;
  }

  DOUBLE xinv=1.0/xbinsize;
  DOUBLE yinv=1.0/ybinsize;
  DOUBLE zinv=1.0/zbinsize;

  for (int64_t i=0;i<np;i++)  {
    ix=(int)((x[i]-xmin)*xinv) ;
    iy=(int)((y[i]-ymin)*yinv) ;
    iz=(int)((z[i]-zmin)*zinv) ;
    if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
    if (iy>nmesh_y-1)  iy--;
    if (iz>nmesh_z-1)  iz--;
	if(! ( ix >= 0 && ix < nmesh_x && iy >=0 && iy < nmesh_y && iz >= 0 && iz < nmesh_z)) {
	  fprintf(stderr,"Problem with i = %"PRId64" x = %lf y = %lf z = %lf \n",i,x[i],y[i],z[i]);
	  fprintf(stderr,"ix = %d iy = %d iz = %d\n",ix,iy,iz);
	}
	assert(x[i] >= xmin && x[i] <= xmax && "x-position is within limits");
	assert(y[i] >= ymin && y[i] <= ymax && "y-position is within limits");
	assert(z[i] >= zmin && z[i] <= zmax && "z-position is within limits");
	
	assert(ix >= 0 && ix < nmesh_x && "ix is in range");
    assert(iy >= 0 && iy < nmesh_y && "iy is in range");
    assert(iz >= 0 && iz < nmesh_z && "iz is in range");

	int64_t index = ix*nmesh_y*nmesh_z + iy*nmesh_z + iz;

    if(lattice[index].nelements == nallocated[index]) {
      expected_n = nallocated[index]*MEMORY_INCREASE_FAC;

	  //In case expected_n is 1 or MEMORY_INCREASE_FAC is 1. 
	  //This way, we only increase by a very few particles 
	  // at a time. Smaller memory footprint
      while(expected_n <= nallocated[index] || ((expected_n % NVEC) != 0))
				expected_n++;

			const size_t memsize=sizeof(DOUBLE);
			lattice[index].x = my_realloc(lattice[index].x ,memsize,expected_n,"lattice.x");
			lattice[index].y = my_realloc(lattice[index].y ,memsize,expected_n,"lattice.y");
			lattice[index].z = my_realloc(lattice[index].z ,memsize,expected_n,"lattice.z");
			nallocated[index] = expected_n;
    }
    assert(lattice[index].nelements < nallocated[index] && "Ensuring that number of particles in a cell doesn't corrupt memory");
		const int64_t ipos = lattice[index].nelements;
		lattice[index].x[ipos] = x[i];
		lattice[index].y[ipos] = y[i];
		lattice[index].z[ipos] = z[i];
		lattice[index].nelements++;
	}
  free(nallocated);

  //You can free the extra memory reserved by the mallocs by looping over totncells and doing a realloc(lattice[index].x,sizeof(DOUBLE),lattice[index].nelements,"lattice.x")
  
  *nlattice_x=nmesh_x;
  *nlattice_y=nmesh_y;
  *nlattice_z=nmesh_z;
#ifndef SILENT	
  gettimeofday(&t1,NULL);
  fprintf(stderr," Time taken = %6.2lf sec\n",ADD_DIFF_TIME(t0,t1));
#endif  
  return lattice;
}

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
															 int *nlattice_z)
{
  cellarray_nvec *lattice=NULL;
  int ix,iy,iz;
  int nmesh_x,nmesh_y,nmesh_z;
  int64_t *nallocated=NULL;
  DOUBLE xdiff,ydiff,zdiff;
  DOUBLE cell_volume,box_volume;
  DOUBLE xbinsize,ybinsize,zbinsize;
  int64_t expected_n=0;
  int64_t totncells;

#ifndef SILENT	
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
#endif
	
  xbinsize = get_binsize(xmin,xmax,max_x_size,xbin_refine_factor, NLATMAX, &nmesh_x);
  ybinsize = get_binsize(ymin,ymax,max_y_size,ybin_refine_factor, NLATMAX, &nmesh_y);
  zbinsize = get_binsize(zmin,zmax,max_z_size,zbin_refine_factor, NLATMAX, &nmesh_z);
  
  totncells = (int64_t) nmesh_x * (int64_t) nmesh_y * (int64_t) nmesh_z;

  xdiff = xmax-xmin;
  ydiff = ymax-ymin;
  zdiff = zmax-zmin;
  
  cell_volume=xbinsize*ybinsize*zbinsize;
  box_volume=xdiff*ydiff*zdiff;
  expected_n=(int64_t)(np*cell_volume/box_volume*MEMORY_INCREASE_FAC);
  expected_n=expected_n < NVEC ? NVEC:expected_n;
	while((expected_n % NVEC) != 0)
		expected_n++;
	
#ifndef SILENT	
  fprintf(stderr,"In %s> Running with [nmesh_x, nmesh_y, nmesh_z]  = %d,%d,%d. ",__FUNCTION__,nmesh_x,nmesh_y,nmesh_z);
#endif	
  lattice    = (cellarray_nvec *) my_malloc(sizeof(cellarray_nvec), totncells);
  nallocated = (int64_t *)       my_malloc(sizeof(*nallocated)      , totncells);

  /*
	Allocate memory for each of the fields in cellarray_nvec. Since we haven't processed the data yet, 
	expected_n is a reasonable guess as to the number of points in the cell. 
   */
  for (int64_t index=0;index<totncells;index++) {
		const size_t memsize=3*sizeof(DOUBLE);
		lattice[index].pos = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
		lattice[index].nelements=0;
		nallocated[index] = expected_n;
  }

  DOUBLE xinv=1.0/xbinsize;
  DOUBLE yinv=1.0/ybinsize;
  DOUBLE zinv=1.0/zbinsize;

  for (int64_t i=0;i<np;i++)  {
    ix=(int)((x[i]-xmin)*xinv) ;
    iy=(int)((y[i]-ymin)*yinv) ;
    iz=(int)((z[i]-zmin)*zinv) ;
    if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
    if (iy>nmesh_y-1)  iy--;
    if (iz>nmesh_z-1)  iz--;
	if(! ( ix >= 0 && ix < nmesh_x && iy >=0 && iy < nmesh_y && iz >= 0 && iz < nmesh_z)) {
	  fprintf(stderr,"Problem with i = %"PRId64" x = %lf y = %lf z = %lf \n",i,x[i],y[i],z[i]);
	  fprintf(stderr,"ix = %d iy = %d iz = %d\n",ix,iy,iz);
	}
	assert(x[i] >= xmin && x[i] <= xmax && "x-position is within limits");
	assert(y[i] >= ymin && y[i] <= ymax && "y-position is within limits");
	assert(z[i] >= zmin && z[i] <= zmax && "z-position is within limits");
	
	assert(ix >= 0 && ix < nmesh_x && "ix is in range");
    assert(iy >= 0 && iy < nmesh_y && "iy is in range");
    assert(iz >= 0 && iz < nmesh_z && "iz is in range");

	int64_t index = ix*nmesh_y*nmesh_z + iy*nmesh_z + iz;

    if(lattice[index].nelements == nallocated[index]) {
      expected_n = nallocated[index]*MEMORY_INCREASE_FAC;

	  //In case expected_n is 1 or MEMORY_INCREASE_FAC is 1. 
	  //This way, we only increase by a very few particles 
	  // at a time. Smaller memory footprint
      while(expected_n <= nallocated[index] || ((expected_n % NVEC) != 0))
				expected_n++;

			const size_t memsize=3*sizeof(DOUBLE);
			lattice[index].pos = my_realloc(lattice[index].pos ,memsize,expected_n,"lattice.pos");
			nallocated[index] = expected_n;
    }
    assert(lattice[index].nelements < nallocated[index] && "Ensuring that number of particles in a cell doesn't corrupt memory");
		const int num_nvec_bunch = lattice[index].nelements/NVEC;
		const size_t xoffset = num_nvec_bunch * NVEC * 3;
		const size_t yoffset = xoffset + NVEC;
		const size_t zoffset = xoffset + 2*NVEC;
		const int ipos=lattice[index].nelements % NVEC;
		DOUBLE *xpos = &(lattice[index].pos[xoffset]);
		DOUBLE *ypos = &(lattice[index].pos[yoffset]);
		DOUBLE *zpos = &(lattice[index].pos[zoffset]);
		xpos[ipos] = x[i];
		ypos[ipos] = y[i];
		zpos[ipos] = z[i];
		lattice[index].nelements++;
	}
  free(nallocated);

  //You can free the extra memory reserved by the mallocs by looping over totncells and doing a realloc(lattice[index].x,sizeof(DOUBLE),lattice[index].nelements,"lattice.x")
  
  *nlattice_x=nmesh_x;
  *nlattice_y=nmesh_y;
  *nlattice_z=nmesh_z;
#ifndef SILENT	
  gettimeofday(&t1,NULL);
  fprintf(stderr," Time taken = %6.2lf sec\n",ADD_DIFF_TIME(t0,t1));
#endif  
  return lattice;
}

