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
#include <stdbool.h>

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

void free_cellarray(cellarray *lattice, const int64_t totncells)
{
    for(int64_t i=0;i<totncells;i++) {
        free(lattice[i].x);
        free(lattice[i].y);
        free(lattice[i].z);
    }
    free(lattice);
}

void free_cellarray_nvec(cellarray_nvec *lattice, const int64_t totncells)
{
    for(int64_t i=0;i<totncells;i++) {
        free(lattice[i].pos);
    }
    free(lattice);
}

void free_cellarray_index(cellarray_index *lattice, const int64_t totncells, const int free_wraps)
{
    
    for(int64_t i=0;i<totncells;i++){
        if(lattice[i].nelements == 0) continue;

        if(free_wraps == 1) {
            free(lattice[i].xwrap);
            free(lattice[i].ywrap);
            free(lattice[i].zwrap);
        }
        free(lattice[i].ngb_cells);
    }
    free(lattice);
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
        int ix=(int)((x[i]-xmin)*xinv) ;
        int iy=(int)((y[i]-ymin)*yinv) ;
        int iz=(int)((z[i]-zmin)*zinv) ;
        if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
        if (iy>nmesh_y-1)  iy--;
        if (iz>nmesh_z-1)  iz--;
        XASSERT(x[i] >= xmin && x[i] <= xmax,
                "x[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, x[i], xmin, xmax);
        XASSERT(y[i] >= ymin && y[i] <= ymax,
                "y[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, y[i], ymin, ymax);
        XASSERT(z[i] >= zmin && z[i] <= zmax,
                "z[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, z[i], zmin, zmax);

        XASSERT(ix >= 0 && ix < nmesh_x, "ix=%d must be within [0,%d)\n", ix, nmesh_x);
        XASSERT(iy >= 0 && iy < nmesh_y, "iy=%d must be within [0,%d)\n", iy, nmesh_y);
        XASSERT(iz >= 0 && iz < nmesh_z, "iz=%d must be within [0,%d)\n", iz, nmesh_z);

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
        XASSERT(lattice[index].nelements < nallocated[index],
                ANSI_COLOR_RED"BUG: lattice[%"PRId64"].nelements = %"PRId64" must be less than allocated memory = %"PRId64 ANSI_COLOR_RESET"\n",
                index, lattice[index].nelements, nallocated[index]);

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
        int ix=(int)((x[i]-xmin)*xinv) ;
        int iy=(int)((y[i]-ymin)*yinv) ;
        int iz=(int)((z[i]-zmin)*zinv) ;
        if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
        if (iy>nmesh_y-1)  iy--;
        if (iz>nmesh_z-1)  iz--;
        if(! ( ix >= 0 && ix < nmesh_x && iy >=0 && iy < nmesh_y && iz >= 0 && iz < nmesh_z)) {
            fprintf(stderr,"Problem with i = %"PRId64" x = %lf y = %lf z = %lf \n",i,x[i],y[i],z[i]);
            fprintf(stderr,"ix = %d iy = %d iz = %d\n",ix,iy,iz);
        }
        XASSERT(x[i] >= xmin && x[i] <= xmax,
                "x[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, x[i], xmin, xmax);
        XASSERT(y[i] >= ymin && y[i] <= ymax,
                "y[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, y[i], ymin, ymax);
        XASSERT(z[i] >= zmin && z[i] <= zmax,
                "z[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, z[i], zmin, zmax);

        XASSERT(ix >= 0 && ix < nmesh_x, "ix=%d must be within [0,%d)\n", ix, nmesh_x);
        XASSERT(iy >= 0 && iy < nmesh_y, "iy=%d must be within [0,%d)\n", iy, nmesh_y);
        XASSERT(iz >= 0 && iz < nmesh_z, "iz=%d must be within [0,%d)\n", iz, nmesh_z);

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
        XASSERT(lattice[index].nelements < nallocated[index],
                ANSI_COLOR_RED"BUG: lattice[%"PRId64"].nelements = %"PRId64" must be less than allocated memory = %"PRId64 ANSI_COLOR_RESET"\n",
                index, lattice[index].nelements, nallocated[index]);

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


/* Need SGLIB to simultaneously sort the particles */
#include "sglib.h"

struct cellarray_index * gridlink_index(const int64_t np,
                                 DOUBLE *x, DOUBLE *y, DOUBLE *z,
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

    int nmesh_x,nmesh_y,nmesh_z;
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
    struct cellarray_index *lattice  = (struct cellarray_index *) my_malloc(sizeof(*lattice), totncells);
    for(int64_t i=0;i<totncells;i++) {
      lattice[i].start = LONG_MAX;//Hoping to crash the code if I incorrectly try to access this!
      lattice[i].nelements = 0;
    }

    int64_t *cell_index = my_malloc(sizeof(*cell_index), np);//each particle needs to get a cell index.

    /*
      Allocate memory for each of the fields in struct cellarray_index. Since we haven't processed the data yet,
      expected_n is a reasonable guess as to the number of points in the cell.
    */

    const DOUBLE xinv=1.0/xbinsize;
    const DOUBLE yinv=1.0/ybinsize;
    const DOUBLE zinv=1.0/zbinsize;

    for (int64_t i=0;i<np;i++)  {
        int ix=(int)((x[i]-xmin)*xinv) ;
        int iy=(int)((y[i]-ymin)*yinv) ;
        int iz=(int)((z[i]-zmin)*zinv) ;
        if (ix>nmesh_x-1)  ix--;    /* this shouldn't happen, but . . . */
        if (iy>nmesh_y-1)  iy--;
        if (iz>nmesh_z-1)  iz--;
        XASSERT(x[i] >= xmin && x[i] <= xmax,
                "x[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, x[i], xmin, xmax);
        XASSERT(y[i] >= ymin && y[i] <= ymax,
                "y[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, y[i], ymin, ymax);
        XASSERT(z[i] >= zmin && z[i] <= zmax,
                "z[%"PRId64"] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, z[i], zmin, zmax);

        XASSERT(ix >= 0 && ix < nmesh_x, "ix=%d must be within [0,%d)\n", ix, nmesh_x);
        XASSERT(iy >= 0 && iy < nmesh_y, "iy=%d must be within [0,%d)\n", iy, nmesh_y);
        XASSERT(iz >= 0 && iz < nmesh_z, "iz=%d must be within [0,%d)\n", iz, nmesh_z);

        const int64_t index = ix*nmesh_y*nmesh_z + (int64_t) iy*nmesh_z + (int64_t) iz;
        cell_index[i] = index;
    }

    /* Now sort the particles based on cell_index */
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,x,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,y,i,j);                   \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,z,i,j);\
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64_t,cell_index,i,j)}

    SGLIB_ARRAY_QUICK_SORT(int64_t, cell_index, np, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

    int64_t start_cell = cell_index[0];
    lattice[start_cell].start=0;
    lattice[start_cell].nelements=1;
    for(int64_t i=1;i<np;i++) {
        const int64_t icell = cell_index[i];
        if(icell != start_cell) {
            lattice[icell].start = i;
            lattice[icell].nelements = 1;
            start_cell = icell;
        } else {
            lattice[icell].nelements++;
        }
    }
    free(cell_index);
    
    *nlattice_x=nmesh_x;
    *nlattice_y=nmesh_y;
    *nlattice_z=nmesh_z;

#ifndef SILENT
    gettimeofday(&t1,NULL);
    fprintf(stderr," Time taken = %6.2lf sec\n",ADD_DIFF_TIME(t0,t1));
#endif
    return lattice;
}


void assign_ngb_cells(struct cellarray_index *lattice1, struct cellarray_index *lattice2, const int64_t totncells,
                      const int xbin_refine_factor, const int ybin_refine_factor, const int zbin_refine_factor,
                      const int nmesh_x, const int nmesh_y, const int nmesh_z,
                      const DOUBLE xdiff, const DOUBLE ydiff, const DOUBLE zdiff, 
                      const int autocorr, const int periodic)
{
    /* Now figure out the neighbouring cells*/
    //First create a giant list of cells opened
    bool *opened=NULL;
    if(autocorr == 1) {
        opened = my_malloc(sizeof(bool), totncells * totncells);
        for(int64_t i=0;i < totncells * totncells; i++) {
            opened[i] = false;
        }
    }

    const int64_t nx_ngb = 2*xbin_refine_factor + 1;
    const int64_t ny_ngb = 2*ybin_refine_factor + 1;
    const int64_t nz_ngb = (autocorr == 0) ? 2*zbin_refine_factor + 1: zbin_refine_factor+1;
    const int64_t max_ngb_cells = nx_ngb * ny_ngb * nz_ngb;

    for(int64_t icell=0;icell<totncells;icell++) {
        struct cellarray_index *first = &(lattice1[icell]);
        if(first->nelements == 0) continue;
        const int iz = icell % nmesh_z;
        const int ix = icell / (nmesh_y * nmesh_z );
        const int iy = (icell - iz - ix*nmesh_z*nmesh_y)/nmesh_z;
        XASSERT(icell == (ix * nmesh_y * nmesh_z + iy * nmesh_z + (int64_t) iz),
                ANSI_COLOR_RED"BUG: Index reconstruction is wrong. icell = %"PRId64" reconstructed index = %"PRId64 ANSI_COLOR_RESET"\n",
                icell, (ix * nmesh_y * nmesh_z + iy * nmesh_z + (int64_t) iz));
        

        first->num_ngb = 0;
        if(periodic == 1) {
            first->xwrap = my_malloc(sizeof(*(first->xwrap)), max_ngb_cells);
            first->ywrap = my_malloc(sizeof(*(first->ywrap)), max_ngb_cells);
            first->zwrap = my_malloc(sizeof(*(first->zwrap)), max_ngb_cells);
        } else {
            first->xwrap = NULL;
            first->ywrap = NULL;
            first->zwrap = NULL;
        }
        first->ngb_cells = my_malloc(sizeof(*(first->ngb_cells)) , max_ngb_cells);

        for(int iix=-xbin_refine_factor;iix<=xbin_refine_factor;iix++){
            const int periodic_ix = (ix + iix + nmesh_x) % nmesh_x;
            const int non_periodic_ix = ix + iix;
            const int iiix = (periodic == 1) ? periodic_ix:non_periodic_ix;
            if(iiix < 0 || iiix >= nmesh_x) continue;
            const DOUBLE off_xwrap = ((ix + iix) >= 0) && ((ix + iix) < nmesh_x) ? 0.0: ((ix+iix) < 0 ? xdiff:-xdiff);
            
            for(int iiy=-ybin_refine_factor;iiy<=ybin_refine_factor;iiy++) {
                const int periodic_iy = (iy + iiy + nmesh_y) % nmesh_y;
                const int non_periodic_iy = iy + iiy;
                const int iiiy = (periodic == 1) ? periodic_iy:non_periodic_iy;
                if(iiiy < 0 || iiiy >= nmesh_y) continue;
                const DOUBLE off_ywrap = ((iy + iiy) >= 0) && ((iy + iiy) < nmesh_y) ? 0.0: ((iy+iiy) < 0 ? ydiff:-ydiff);
                const int start_iz = (autocorr == 0) ? -zbin_refine_factor:0;
                for(int64_t iiz=start_iz;iiz<=zbin_refine_factor;iiz++){
                    const int periodic_iz = (iz + iiz + nmesh_z) % nmesh_z;
                    const int non_periodic_iz = iz + iiz;
                    const int iiiz = (periodic == 1) ? periodic_iz:non_periodic_iz;
                    if(iiiz < 0 || iiiz >= nmesh_z) continue;

                    const DOUBLE off_zwrap = ((iz + iiz) >= 0) && ((iz + iiz) < nmesh_z) ? 0.0: ((iz+iiz) < 0 ? zdiff:-zdiff);
                    const int64_t icell2 = iiiz + (int64_t) nmesh_z*iiiy + nmesh_z*nmesh_y*iiix;

                    //For cases where we are not double-counting (i.e., wp and xi), the same-cell
                    //must always be evaluated. In all other cases, (i.e., where double-counting is occurring)
                    //is used, include that in the ngb_cells! The interface is a lot cleaner in the double-counting
                    //kernels in that case!
                    if(autocorr == 1 && icell2 == icell) {
                        continue;
                    }
                    const int64_t index1 = icell2 * totncells + icell;
                    const int64_t index2 = icell * totncells + icell2;
                    
                    if(autocorr == 1 && opened[index1] == true) {
                        continue;
                    }

                    const int64_t ngb_index = first->num_ngb;
                    XASSERT(ngb_index < max_ngb_cells,"ngb index = %"PRId64" should be less than max_ngb = %"PRId64"\n", ngb_index, max_ngb_cells);
                    first->ngb_cells[ngb_index] = &(lattice2[icell2]);

                    //Note the xwrap/ywraps do not have memory allocated for them in the
                    //non-periodic case. 
                    if(periodic == 1) {
                        first->xwrap[ngb_index] = off_xwrap;
                        first->ywrap[ngb_index] = off_ywrap;
                        first->zwrap[ngb_index] = off_zwrap;
                    } 
                    first->num_ngb++;

                    if(autocorr == 1) {
                        opened[index1] = true;
                        opened[index2] = true;
                    }
                }
            }
        }

    }
    if(autocorr == 1) {
        free(opened);
    }
    
}    

