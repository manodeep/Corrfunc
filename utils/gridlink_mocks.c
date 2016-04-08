/* File: gridlink_mocks.c */
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
#include "function_precision.h"
#include "cellarray_mocks.h"
#include "utils.h"
#include "gridlink_mocks.h"

#define MEMORY_INCREASE_FAC   1.1


void get_max_min_data(const int64_t ND1, const DOUBLE * restrict cz,
                      DOUBLE *min_cz, DOUBLE *max_cz
#ifdef LINK_IN_DEC
                      ,const DOUBLE * restrict dec,
                      DOUBLE *min_dec, DOUBLE *max_dec
#endif

                      )
{
    DOUBLE czmin = *min_cz, czmax=*max_cz;
#ifdef LINK_IN_DEC
    DOUBLE dec_min = *min_dec, dec_max = *max_dec;
#endif

    for(int64_t i=0;i<ND1;i++) {
        if(cz[i] < czmin) czmin=cz[i];
        if(cz[i] > czmax) czmax=cz[i];

#ifdef LINK_IN_DEC
        if(dec[i] < dec_min) dec_min=dec[i];
        if(dec[i] > dec_max) dec_max=dec[i];
#endif

    }
    *min_cz=czmin;
    *max_cz=czmax;

#ifdef LINK_IN_DEC
    *min_dec = dec_min;
    *max_dec = dec_max;
#endif

}



cellarray_mocks *gridlink1D(const int64_t np,const DOUBLE czmin,const DOUBLE czmax, const DOUBLE pimax,
                            const DOUBLE *dec, const DOUBLE *ra, const DOUBLE *cz,
                            int *ngrid,int *max_in_cell,
                            const int zbin_refine_factor)
{
    int nmesh,max_n;
    const DOUBLE czdiff = (czmax-czmin);
    XASSERT(czdiff > 0.0, "There needs to be some width in cz for the data. Current width = %"DOUBLE_FORMAT"\n", czdiff);
    const DOUBLE inv_czdiff=1.0/czdiff;

    /* Instead of directly assigning to int nmesh via
       truncation, I have chosen this two-step process to
       avoid int overflow (and negative nmesh) for very small
       values of pimax (or some combination of the three
       variables in the next line).
    */
    const DOUBLE this_nmesh = zbin_refine_factor*czdiff/pimax ;
    nmesh = this_nmesh > NLATMAX ? NLATMAX:(int) this_nmesh;
    *ngrid=nmesh ;
    const int64_t totncells = nmesh;


#ifndef SILENT
    struct timeval t0,t1;
    gettimeofday(&t0,NULL);
#endif


    //Because we are only binning in 1-D, I have left exoected_n as a 64 bit integer.
    //However, both nallocated and elements in the lattice structure are ints.
    int64_t expected_n=(int64_t)((np/(double) nmesh)  *MEMORY_INCREASE_FAC);
    expected_n=expected_n < NVEC ? NVEC:expected_n;
    while(expected_n % NVEC != 0) {
        expected_n++;
    }

    cellarray_mocks *lattice = my_malloc(sizeof(cellarray_mocks), totncells);

    /*
      Allocate memory for each of the fields in cellarray. Since we haven't processed the data yet,
      expected_n is a reasonable guess as to the number of points in the cell.
    */
    for (int64_t index=0;index<totncells;index++) {
        const size_t memsize=3*sizeof(DOUBLE);
        lattice[index].pos = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
        lattice[index].nelements=0;
        lattice[index].nallocated = expected_n;
    }

    max_n=0;
    /*---Loop-over-particles-and-build-grid-arrays----*/
    for(int64_t i=0;i<np;i++) {
        int iz = (int)(nmesh*(cz[i]-czmin)*inv_czdiff) ;
        if (iz >= nmesh) iz--;
        XASSERT(iz >= 0 && iz < nmesh, "iz (cz bin index) = %d must be within [0,%d)\n", iz, nmesh);

        const int64_t index = iz;
        if(lattice[index].nelements == lattice[index].nallocated) {
            expected_n = lattice[index].nallocated*MEMORY_INCREASE_FAC;

            //In case expected_n is 1 or MEMORY_INCREASE_FAC is 1.
            //This way, we only increase by a very few particles
            // at a time. Smaller memory footprint
            while(expected_n <= lattice[index].nelements || expected_n % NVEC != 0) {
                expected_n++;
            }

            const size_t memsize=3*sizeof(DOUBLE);
            lattice[index].pos  = my_realloc(lattice[index].pos ,memsize,expected_n,"lattice.pos");
            lattice[index].nallocated = expected_n;
        }
        XASSERT(lattice[index].nelements < lattice[index].nallocated,
                ANSI_COLOR_RED"BUG: lattice[%"PRId64"].nelements = %d must be less than allocated memory = %d" ANSI_COLOR_RESET"\n",
                index, lattice[index].nelements, lattice[index].nallocated);

        /*
          The particles are stored like this:
          x[NVEC], y{NVEC], z[NVEC], x[NVEC],y[NVEC].....

          So first thing is to find how many xyz chunks have been written (num_nvec_bunch)
          For each chunk of these xyz quantities, we have an offset of 3*NVEC elements for x.
          y is further offset by another NVEC elements (the set of x's in the current chunk)
          And so on.

          After that, we need to figure out the position within the NVEC block for each of
          xyz -- given by ipos. Now we have all the pieces we need to insert the new element.

        */

        const int num_nvec_bunch = lattice[index].nelements/NVEC;
        const int xoffset  = num_nvec_bunch * NVEC * 3;
        const int yoffset  = xoffset + NVEC;
        const int zoffset  = xoffset + 2*NVEC;
        /*      const int czoffset = xoffset + 3*NVEC;           */

        const int ipos=lattice[index].nelements % NVEC;

        DOUBLE *xpos  = &(lattice[index].pos[xoffset]);
        DOUBLE *ypos  = &(lattice[index].pos[yoffset]);
        DOUBLE *zpos  = &(lattice[index].pos[zoffset]);
        /*      DOUBLE *czpos = &(lattice[index].pos[czoffset]); */

        xpos[ipos]  = cz[i]*COSD(dec[i])*COSD(ra[i]);
        ypos[ipos]  = cz[i]*COSD(dec[i])*SIND(ra[i]);
        zpos[ipos]  = cz[i]*SIND(dec[i]);
        /*      czpos[ipos] = cz[i]; */

        lattice[index].nelements++;
        if(lattice[index].nelements > max_n) {
            max_n = lattice[index].nelements;
        }
    }
    *max_in_cell = max_n;
#ifndef SILENT
    gettimeofday(&t1,NULL);
    fprintf(stderr,"%s> Allocating %0.2g (MB) memory for the lattice, expected_n = %"PRId64" nmesh_cz = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,(3*4)*expected_n*nmesh/(1024.*1024.),expected_n,
            nmesh,np,ADD_DIFF_TIME(t0,t1));
#endif
    return lattice;
}

#ifdef LINK_IN_DEC

cellarray_mocks **gridlink2D(const int64_t np,
                             const DOUBLE czmin, const DOUBLE czmax, const DOUBLE pimax,
                             const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE rpmax,
                             const DOUBLE *cz,const DOUBLE *dec, const DOUBLE *ra,
                             int *ngrid_cz,
                             int **ngrid_declination,
                             int *max_in_cell,
                             const int rbin_refine_factor,
                             const int zbin_refine_factor)

{
    int nmesh_cz;
    int max_nmesh_dec;
    int expected_n,max_n;
    size_t totnbytes=0;
    int *ngrid_dec = NULL;

    const DOUBLE czdiff = czmax-czmin;
    const DOUBLE dec_diff = dec_max-dec_min;
    const DOUBLE inv_dec_diff = 1.0/dec_diff;

    cellarray_mocks **lattice=NULL;
    int assigned_n=0;

#ifndef SILENT
    struct timeval t0,t1;
    gettimeofday(&t0,NULL);
#endif
    XASSERT(czdiff > 0.0, "There needs to be some width in cz for the data. Current width = %"DOUBLE_FORMAT"\n", czdiff);
    XASSERT(pimax > 0.0, "Minimum los separation = %"DOUBLE_FORMAT" must be non-zero\n", pimax);
    XASSERT(dec_diff > 0.0, "All of the points can not be at the same declination. Declination difference = %"DOUBLE_FORMAT" must be non-zero\n", dec_diff);
    XASSERT(MEMORY_INCREASE_FAC >= 1.0, "Memory increase factor = %lf must be >=1 \n",MEMORY_INCREASE_FAC);

    //Written this way to work around INT overflows
    const DOUBLE this_nmesh = zbin_refine_factor*czdiff/pimax ;
    nmesh_cz = this_nmesh > NLATMAX ? NLATMAX:(int) this_nmesh;

    *ngrid_cz=nmesh_cz ;
    const DOUBLE cz_binsize = czdiff/nmesh_cz;
    const DOUBLE inv_cz_binsize = 1.0/cz_binsize;

    *ngrid_declination = my_malloc(sizeof(*ngrid_dec),nmesh_cz);
    ngrid_dec = *ngrid_declination;

    /* Find the max. number of declination cells that can be */
    max_nmesh_dec=0;
    for(int i=0;i<nmesh_cz;i++) {
        const int min_iz = (i - zbin_refine_factor) < 0  ? 0:i-zbin_refine_factor;
        const DOUBLE dmin = czmin + 0.5*(min_iz+i)*cz_binsize;//Really I am taking the average of the left edges for the cz bins corresponding to i (= czmin + i*cz_binsize) and min_iz (= czmin + min_iz*cz_binsize).
        const DOUBLE dec_cell = ASIN(rpmax/(2*dmin))*2.0*INV_PI_OVER_180;
        XASSERT(dec_cell > 0.0, "Declination binsize=%"DOUBLE_FORMAT" in cz bin = %d must be positive\n", dec_cell, i);
        const DOUBLE this_nmesh_dec = dec_diff*rbin_refine_factor/dec_cell;
        const int nmesh_dec = this_nmesh_dec > NLATMAX ? NLATMAX:(int) this_nmesh_dec;
        if(nmesh_dec > max_nmesh_dec) max_nmesh_dec = nmesh_dec;
        ngrid_dec[i]=nmesh_dec ;
    }

    //We need to create a matrix with nmesh_cz and max_nmesh_dec row/columns
    //The crude estimate of the average number of points per cell
    expected_n=(int)( (np/(DOUBLE) (nmesh_cz*max_nmesh_dec)) *MEMORY_INCREASE_FAC);
    expected_n = expected_n < NVEC ? NVEC:expected_n;
    while(expected_n % NVEC != 0) {
        expected_n++;
    }

    /*---Allocate-and-initialize-grid-arrays----------*/
    lattice = (cellarray_mocks **) matrix_malloc(sizeof(cellarray_mocks),nmesh_cz,max_nmesh_dec); //This allocates extra and is wasteful
    totnbytes += nmesh_cz*max_nmesh_dec*sizeof(cellarray_mocks);
    for(int i=0;i<nmesh_cz;i++) {
        const int nmesh_dec = ngrid_dec[i];
        for(int j=0;j<nmesh_dec;j++) {
            const size_t memsize = 3*sizeof(DOUBLE);
            lattice[i][j].pos    = my_malloc(memsize,expected_n);
            lattice[i][j].nelements=0;
            lattice[i][j].nallocated=expected_n;
            totnbytes += memsize*expected_n;
        }
    }


    max_n = 0;
    /*---Loop-over-particles-and-build-grid-arrays----*/
    for(int i=0;i<np;i++) {
        int iz = (int)((cz[i]-czmin)*inv_cz_binsize) ;
        if (iz >= nmesh_cz) iz--;
        XASSERT(iz >= 0 && iz < nmesh_cz, "iz (cz bin index) = %d must be within [0,%d)\n", iz, nmesh_cz);
        XASSERT(dec[i] >= dec_min && dec[i] <= dec_max,
                "dec[%d] = %"DOUBLE_FORMAT" must be within [%"DOUBLE_FORMAT",%"DOUBLE_FORMAT"]\n",
                i, dec[i], dec_min, dec_max);
        int idec = (int)(ngrid_dec[iz]*(dec[i]-dec_min)*inv_dec_diff);
        if(idec >= ngrid_dec[iz]) idec--;
        XASSERT(idec >= 0 && idec < ngrid_dec[iz],
                "idec (dec bin index) = %d must be within [0,%d) for iz (zbin index) = %d\n",
                idec, ngrid_dec[iz], iz);

        if(lattice[iz][idec].nelements == lattice[iz][idec].nallocated) {
            expected_n = lattice[iz][idec].nallocated*MEMORY_INCREASE_FAC;
            while(expected_n <= lattice[iz][idec].nelements || expected_n % NVEC != 0) {
                expected_n++;
            }
            const size_t memsize = 3*sizeof(DOUBLE);
            lattice[iz][idec].pos  = my_realloc(lattice[iz][idec].pos ,memsize,expected_n,"lattice.pos");
            lattice[iz][idec].nallocated = expected_n;
        }
        XASSERT(lattice[iz][idec].nelements < lattice[iz][idec].nallocated,
                ANSI_COLOR_RED"BUG: lattice[%d][%d].nelements = %d must be less than allocated memory = %d" ANSI_COLOR_RESET"\n",
                iz, idec, lattice[iz][idec].nelements, lattice[iz][idec].nallocated);

        const int num_nvec_bunch = lattice[iz][idec].nelements/NVEC;
        const size_t xoffset  = num_nvec_bunch * NVEC * 3;
        const size_t yoffset  = xoffset + NVEC;
        const size_t zoffset  = xoffset + 2*NVEC;
        /*  const size_t czoffset = xoffset + 3*NVEC;            */

        const int ipos=lattice[iz][idec].nelements % NVEC;

        DOUBLE *xpos  = &(lattice[iz][idec].pos[xoffset]);
        DOUBLE *ypos  = &(lattice[iz][idec].pos[yoffset]);
        DOUBLE *zpos  = &(lattice[iz][idec].pos[zoffset]);
        /*  DOUBLE *czpos = &(lattice[iz][idec].pos[czoffset]); */

        xpos[ipos]  = cz[i]*COSD(dec[i])*COSD(ra[i]);
        ypos[ipos]  = cz[i]*COSD(dec[i])*SIND(ra[i]);
        zpos[ipos]  = cz[i]*SIND(dec[i]);
        /*  czpos[ipos] = cz[i]; */
        lattice[iz][idec].nelements++;
        if(lattice[iz][idec].nelements > max_n) {
            max_n = lattice[iz][idec].nelements;
        }
        assigned_n++;
    }
    *max_in_cell = max_n;
#ifndef SILENT
    gettimeofday(&t1,NULL);
    fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d nmesh_cz = %d max_nmesh_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,nmesh_cz,max_nmesh_dec,np,
            ADD_DIFF_TIME(t0,t1));
#endif

    return lattice;
}





cellarray * gridlink1D_theta(const int64_t np,
                             const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE thetamax,
                             const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1, const DOUBLE  * restrict dec,
                             int *ngrid_declination,
                             int *max_in_cell,
                             const int rbin_refine_factor)

{
    int expected_n,max_n;
    size_t totnbytes=0;

    const DOUBLE dec_diff = dec_max-dec_min;
    const DOUBLE inv_dec_diff = 1.0/dec_diff;

    cellarray *lattice=NULL;
    int assigned_n=0;

#ifndef SILENT
    struct timeval t0,t1;
    gettimeofday(&t0,NULL);
#endif

    XASSERT(thetamax > 0.0, "Minimum angular separation = %"DOUBLE_FORMAT" must be positive\n", thetamax);
    XASSERT(dec_diff > 0.0, "All of the points can not be at the same declination. Declination difference = %"DOUBLE_FORMAT" must be non-zero\n", dec_diff);
    XASSERT(MEMORY_INCREASE_FAC >= 1.0, "Memory increase factor = %lf must be >=1 \n",MEMORY_INCREASE_FAC);
    
    /* Find the max. number of declination cells that can be */
    const DOUBLE this_ngrid_dec = dec_diff*rbin_refine_factor/thetamax;
    const int ngrid_dec = this_ngrid_dec > NLATMAX ? NLATMAX:(int) this_ngrid_dec;

    *ngrid_declination=ngrid_dec;

    expected_n=(int)( (np/(DOUBLE) (ngrid_dec)) *MEMORY_INCREASE_FAC);
    expected_n = expected_n < NVEC ? NVEC:expected_n;
    while(expected_n % NVEC != 0) {
        expected_n++;
    }
    totnbytes += ngrid_dec*sizeof(cellarray);

    /*---Allocate-and-initialize-grid-arrays----------*/
    lattice = (cellarray *) my_malloc(sizeof(cellarray),ngrid_dec); //This allocates extra and is wasteful
    for(int j=0;j<ngrid_dec;j++) {
        const size_t memsize = 3*sizeof(DOUBLE);
        lattice[j].pos     = my_malloc(memsize,expected_n);
        lattice[j].nelements=0;
        lattice[j].nallocated=expected_n;
        totnbytes += memsize*expected_n;
    }


    max_n = 0;
    /*---Loop-over-particles-and-build-grid-arrays----*/
    for(int i=0;i<np;i++) {
        int idec = (int)(ngrid_dec*(dec[i]-dec_min)*inv_dec_diff);
        if(idec >=ngrid_dec) idec--;
        XASSERT(idec >= 0 && idec < ngrid_dec,
                "idec (dec bin index) = %d must be within [0, %d)", idec, ngrid_dec);
        if(lattice[idec].nelements == lattice[idec].nallocated) {
            expected_n = lattice[idec].nallocated*MEMORY_INCREASE_FAC;
            while(expected_n <= lattice[idec].nelements || expected_n % NVEC != 0){
                expected_n++;
            }

            const size_t memsize = 3*sizeof(DOUBLE);
            lattice[idec].pos  = my_realloc(lattice[idec].pos, memsize,expected_n,"lattice.pos");
            lattice[idec].nallocated = expected_n;
        }
        XASSERT(lattice[idec].nelements < lattice[idec].nallocated,
                ANSI_COLOR_RED"BUG: lattice[%d].nelements = %d must be less than allocated memory = %d" ANSI_COLOR_RESET"\n",
                idec, lattice[idec].nelements, lattice[idec].nallocated);

        const int num_nvec_bunch = lattice[idec].nelements/NVEC;
        const size_t xoffset = num_nvec_bunch * NVEC * 3;
        const size_t yoffset = xoffset + NVEC;
        const size_t zoffset = xoffset + 2*NVEC;
        const int ipos=lattice[idec].nelements % NVEC;
        DOUBLE *xpos = &(lattice[idec].pos[xoffset]);
        DOUBLE *ypos = &(lattice[idec].pos[yoffset]);
        DOUBLE *zpos = &(lattice[idec].pos[zoffset]);

        xpos[ipos]  = x1[i];
        ypos[ipos]  = y1[i];
        zpos[ipos]  = z1[i];

        lattice[idec].nelements++;
        if(lattice[idec].nelements > max_n)
            max_n = lattice[idec].nelements;
        assigned_n++;
    }
    *max_in_cell = max_n;
#ifndef SILENT
    gettimeofday(&t1,NULL);
    fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d ngrid_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,ngrid_dec,np,
            ADD_DIFF_TIME(t0,t1));
#endif
    /* fprintf(stderr,"np = %d assigned_n = %d\n",np,assigned_n); */
    return lattice;
}



#ifdef LINK_IN_RA

cellarray_mocks *** gridlink3D(const int64_t np,
                               const DOUBLE czmin,const DOUBLE czmax,const DOUBLE pimax,
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
                               const int zbin_refine_factor)


{
    int nmesh_cz;
    const DOUBLE czdiff = (czmax-czmin);
    int expected_n,max_n;
    size_t totnbytes=0;
    int *ngrid_dec = NULL;

    const DOUBLE dec_diff = dec_max-dec_min;
    DOUBLE inv_dec_diff = 1.0/dec_diff;
    int max_nmesh_dec;
    const DOUBLE phi_diff = phi_max - phi_min;
    DOUBLE inv_phi_diff = 1.0/phi_diff;
    DOUBLE phi_cell=0.0;
    DOUBLE *dec_binsizes=NULL;
    int **ngrid_ra=NULL;

    cellarray_mocks ***lattice=NULL;
#ifndef SILENT
    struct timeval t0,t1;
    gettimeofday(&t0,NULL);
#endif

    XASSERT(czdiff > 0.0, "There needs to be some width in cz for the data. Current width = %"DOUBLE_FORMAT"\n", czdiff);
    XASSERT(pimax > 0.0, "Minimum los separation = %"DOUBLE_FORMAT" must be non-zero\n", pimax);
    XASSERT(dec_diff > 0.0, "All of the points can not be at the same declination. Declination difference = %"DOUBLE_FORMAT" must be non-zero\n", dec_diff);
    XASSERT(phi_diff > 0.0, "All of the points can not be at the same RA. RA difference = %"DOUBLE_FORMAT" must be non-zero\n", phi_diff);
    XASSERT(MEMORY_INCREASE_FAC >= 1.0, "Memory increase factor = %lf must be >=1 \n",MEMORY_INCREASE_FAC);


    const DOUBLE this_nmesh = zbin_refine_factor*czdiff/pimax ;
    nmesh_cz = this_nmesh > NLATMAX ? NLATMAX:(int) this_nmesh;
    *ngrid_cz=nmesh_cz ;

    const DOUBLE cz_binsize = czdiff/nmesh_cz;
    const DOUBLE inv_cz_binsize = 1.0/cz_binsize;

    *ngrid_declination = my_malloc(sizeof(*ngrid_dec),nmesh_cz);
    ngrid_dec = *ngrid_declination;

    dec_binsizes=my_malloc(sizeof(*dec_binsizes),nmesh_cz);

    //First find the max. number of dec bins (max_nmesh_dec)
    max_nmesh_dec = 0;
    for(int iz=0;iz<nmesh_cz;iz++) {
        const int min_iz = iz-zbin_refine_factor < 0 ? 0:iz-zbin_refine_factor;
        const DOUBLE dmin = czmin + 0.5*(min_iz+iz)*cz_binsize;
        const DOUBLE dec_cell =  ASIN(rpmax/(2*dmin))*2.0*INV_PI_OVER_180;// \sigma = 2*arcsin(C/2) -> 2*arcsin( (rpmax/d2min) /2)
        XASSERT(dec_cell > 0.0, "Declination binsize=%"DOUBLE_FORMAT" in cz bin = %d must be positive\n", dec_cell, iz);
        dec_binsizes[iz] = dec_cell;
        const DOUBLE this_nmesh_dec = dec_diff*rbin_refine_factor/dec_cell;
        const int nmesh_dec = this_nmesh_dec > NLATMAX ? NLATMAX:(int) this_nmesh_dec;
        if(nmesh_dec > max_nmesh_dec) max_nmesh_dec = nmesh_dec;
        ngrid_dec[iz]=nmesh_dec ;
    }

    /*---Allocate-and-initialize-grid-arrays----------*/
    *ngrid_phi = (int **) matrix_malloc(sizeof(int), nmesh_cz, max_nmesh_dec);
    ngrid_ra = *ngrid_phi;

    //Now find the maximum number of ra cells
    int max_nmesh_phi=0;
    for(int iz=0;iz<nmesh_cz;iz++) {
        const int nmesh_dec = ngrid_dec[iz];
        DOUBLE dec_binsize=dec_diff/nmesh_dec;
        int min_iz = iz-zbin_refine_factor < 0 ? 0:iz-zbin_refine_factor;
        const DOUBLE dec_cell = dec_binsizes[min_iz];
        const DOUBLE costhetamax=COSD(dec_cell);
        for(int idec=0;idec<nmesh_dec;idec++) {
            int max_idec;
            DOUBLE this_min_dec;
            DOUBLE this_dec = dec_min + idec*dec_binsize;
            if(this_dec > 0) {
                max_idec = idec + rbin_refine_factor >= nmesh_dec ? nmesh_dec-1:idec+rbin_refine_factor;
                this_min_dec = dec_min + (max_idec+1)*dec_binsize;//upper limit for that dec-bin
            } else {
                max_idec = idec - rbin_refine_factor < 0 ? 0:idec-rbin_refine_factor;
                this_min_dec = dec_min + max_idec*dec_binsize;//lower limit for that dec-bin
            }

            /* fprintf(stderr,"min_iz=%d,idec=%d max_idec = %d nmesh_dec = %d ngrid_dec[min_iz] = %d costhetamax = %lf cos(dec_binsize) = %lf \n" */
            /*        ,min_iz,idec,max_idec,nmesh_dec,ngrid_dec[min_iz],costhetamax,COSD(dec_binsize)); */

            phi_cell = 120.0;
            if( (90.0 - FABS(this_min_dec) ) > 1.0) { //make sure min_dec is not close to the pole (within 1 degree)
                DOUBLE sin_min_dec = SIND(this_min_dec),cos_min_dec=COSD(this_min_dec);
                phi_cell = ACOS((costhetamax - sin_min_dec*sin_min_dec)/(cos_min_dec*cos_min_dec))*INV_PI_OVER_180;
                /* phi_cell *= rbin_refine_factor;//My logic does not work - but multiplying with rbin_refine_factor sorts out the problem */
                /* phi_cell *= 1.2;//Still needs a fudge-factor */
                if(!(phi_cell > 0.0)) {
                    /* DOUBLE tmp3 = (costhetamax - sin_min_dec*sin_min_dec)/(cos_min_dec*cos_min_dec); */
                    /* fprintf(stderr,"ERROR: idec = %d max_idec = %d nmesh_dec = %d this_min_dec = %lf dec_cell = %lf phi_cell = %lf is negative. thetamax = %lf sin_min_dec = %lf cos_min_dec = %lf tmp3 = %lf \n", */
                    /*    idec,max_idec,nmesh_dec,this_min_dec,dec_cell,phi_cell,dec_cell,sin_min_dec,cos_min_dec,tmp3); */
                    phi_cell = 120.0;
                }
            }
            XASSERT(phi_cell > 0.0, "RA binsize=%"DOUBLE_FORMAT" in [cz,dec] bin = (%d,%d) must be positive\n", phi_cell, iz, idec);

            phi_cell = phi_cell > 120.0 ? 120.0:phi_cell;
            /* fprintf(stderr,"iz = %4d idec = %4d dec_cell = %6.3lf dec_binsize=%6.2lf this_dec = %6.2lf phi_cell = %7.2lf \n",iz,idec,dec_cell,dec_binsize,dec_min + idec*dec_binsize,phi_cell); */
            const DOUBLE this_nmesh_ra = phi_diff*phibin_refine_factor/phi_cell;
            int nmesh_ra = this_nmesh_ra > NLATMAX ? this_nmesh_ra:(int) this_nmesh_ra;
            if(nmesh_ra < (2*phibin_refine_factor + 1)) {
                nmesh_ra = 2*phibin_refine_factor + 1;
                fprintf(stderr,"%s> Using sub-optimal RA binning to ensure correct functioning of the code\n",__FUNCTION__);
            }
            if(nmesh_ra > max_nmesh_phi) max_nmesh_phi = nmesh_ra;
            ngrid_ra[iz][idec] = nmesh_ra;
        }
    }

    //Allocate the lattice structure
    expected_n=(int)( (np/(DOUBLE) (nmesh_cz*max_nmesh_dec*max_nmesh_phi)) *MEMORY_INCREASE_FAC);
    expected_n = expected_n < NVEC ? NVEC:expected_n;

    //But we are going to store NVEC's each.
    while((expected_n % NVEC) != 0) {
        expected_n++;
    }
    lattice = (cellarray_mocks ***) volume_malloc(sizeof(cellarray_mocks),nmesh_cz,max_nmesh_dec,max_nmesh_phi); //This allocates extra and is wasteful
    totnbytes += nmesh_cz*max_nmesh_dec*max_nmesh_phi*sizeof(cellarray_mocks);


    //Now allocate the memory for positions inside each cell.
    for(int iz=0;iz<nmesh_cz;iz++) {
        const int nmesh_dec = ngrid_dec[iz];
        for(int idec=0;idec<nmesh_dec;idec++) {
            const int nmesh_ra = ngrid_ra[iz][idec];
            for(int ira=0;ira<nmesh_ra;ira++) {
                const size_t memsize=3*sizeof(DOUBLE);//4 pointers for x/y/z
                lattice[iz][idec][ira].pos = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
                lattice[iz][idec][ira].nelements=0;
                lattice[iz][idec][ira].nallocated=expected_n;
                totnbytes += memsize*expected_n;
            }
        }
    }

    max_n = 0;
    /*---Loop-over-particles-and-build-grid-arrays----*/
    for(int i=0;i<np;i++) {
        if(cz[i] >=czmin && cz[i] <= czmax) {
            int iz = (int)((cz[i]-czmin)*inv_cz_binsize) ;
            if (iz >= nmesh_cz) iz--;
            XASSERT(iz >= 0 && iz < nmesh_cz, "iz (cz bin index) = %d must be within [0,%d)\n", iz, nmesh_cz);
            int idec = (int)(ngrid_dec[iz]*(dec[i]-dec_min)*inv_dec_diff);
            if(idec >= ngrid_dec[iz]) idec--;
            XASSERT(idec >=0 && idec < ngrid_dec[iz],
                    "Declination index for particle position = %d must be within [0, %d) for cz-bin = %d\n",
                    idec, ngrid_dec[iz], iz);
            
            
            int ira = (int) (ngrid_ra[iz][idec]*(phi[i]-phi_min)*inv_phi_diff);
            if(ira >= ngrid_ra[iz][idec]) ira--;
            XASSERT(ira >=0 && ira < ngrid_ra[iz][idec],
                    "RA index for particle position = %d must be within [0, %d) for (cz, dec) bin indices = (%d, %d)\n",
                    ira, ngrid_ra[iz][idec], iz, idec);
            
            if(lattice[iz][idec][ira].nelements == lattice[iz][idec][ira].nallocated) {
                expected_n = lattice[iz][idec][ira].nallocated*MEMORY_INCREASE_FAC;
                while(expected_n <= lattice[iz][idec][ira].nelements || expected_n % NVEC != 0) {
                    expected_n++;
                }

                const size_t memsize=3*sizeof(DOUBLE);
                lattice[iz][idec][ira].pos  = my_realloc(lattice[iz][idec][ira].pos ,memsize,expected_n,"lattice.pos");
                lattice[iz][idec][ira].nallocated = expected_n;
            }
            XASSERT(lattice[iz][idec][ira].nelements < lattice[iz][idec][ira].nallocated,
                    ANSI_COLOR_RED"BUG: lattice[%d][%d][%d].nelements = %d must be less than allocated memory = %d" ANSI_COLOR_RESET"\n",
                    iz, idec, ira, lattice[iz][idec][ira].nelements, lattice[iz][idec][ira].nallocated);

            
            const int num_nvec_bunch = lattice[iz][idec][ira].nelements/NVEC;
            const size_t xoffset  = num_nvec_bunch * NVEC * 3;
            const size_t yoffset  = xoffset + NVEC;
            const size_t zoffset  = xoffset + 2*NVEC;
            /*          const size_t czoffset = xoffset + 3*NVEC;            */

            const int ipos=lattice[iz][idec][ira].nelements % NVEC;

            DOUBLE *xpos  = &(lattice[iz][idec][ira].pos[xoffset]);
            DOUBLE *ypos  = &(lattice[iz][idec][ira].pos[yoffset]);
            DOUBLE *zpos  = &(lattice[iz][idec][ira].pos[zoffset]);
            /*          DOUBLE *czpos = &(lattice[iz][idec][ira].pos[czoffset]); */

            xpos[ipos]  = cz[i]*COSD(dec[i])*COSD(phi[i]);
            ypos[ipos]  = cz[i]*COSD(dec[i])*SIND(phi[i]);
            zpos[ipos]  = cz[i]*SIND(dec[i]);
            /*          czpos[ipos] = cz[i]; */

            lattice[iz][idec][ira].nelements++;
            if(lattice[iz][idec][ira].nelements > max_n) {
                max_n = lattice[iz][idec][ira].nelements;
            }
        }
    }
    free(dec_binsizes);
    *max_in_cell = max_n;
#ifndef SILENT
    gettimeofday(&t1,NULL);
    fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d (max_n = %d) nmesh_cz = %d max_nmesh_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,max_n,nmesh_cz,max_nmesh_dec,np, ADD_DIFF_TIME(t0,t1));
#endif
    return lattice;
}




cellarray ** gridlink2D_theta(const int64_t np,
                              const DOUBLE dec_min, const DOUBLE dec_max,const DOUBLE thetamax,
                              const DOUBLE * restrict x1,const DOUBLE * restrict y1,const DOUBLE * restrict z1,
                              const DOUBLE * restrict dec,
                              int *ngrid_declination,
                              const DOUBLE * restrict phi, const DOUBLE phi_min,const DOUBLE phi_max,
                              int **ngrid_phi,
                              int *max_in_cell,
                              const int rbin_refine_factor,
                              const int phibin_refine_factor)
{
    int expected_n,max_n;
    size_t totnbytes=0;
    const DOUBLE dec_diff = dec_max-dec_min;
    const DOUBLE phi_diff = phi_max - phi_min;
    int *ngrid_ra=NULL;

    XASSERT(thetamax > 0.0, "Minimum angular separation = %"DOUBLE_FORMAT" must be positive\n", thetamax);
    XASSERT(dec_diff > 0.0, "All of the points can not be at the same declination. Declination difference = %"DOUBLE_FORMAT" must be non-zero\n", dec_diff);
    XASSERT(phi_diff > 0.0, "All of the points can not be at the same RA. RA difference = %"DOUBLE_FORMAT" must be non-zero\n", phi_diff);
    XASSERT(MEMORY_INCREASE_FAC >= 1.0, "Memory increase factor = %lf must be >=1 \n",MEMORY_INCREASE_FAC);

    const DOUBLE inv_dec_diff = 1.0/dec_diff;
    const DOUBLE inv_phi_diff = 1.0/phi_diff;

    cellarray **lattice=NULL;
#ifndef SILENT
    struct timeval t0,t1;
    gettimeofday(&t0,NULL);
#endif

    const DOUBLE this_ngrid_dec = dec_diff*rbin_refine_factor/thetamax;
    const int ngrid_dec = this_ngrid_dec > NLATMAX ? NLATMAX:(int) this_ngrid_dec;
    *ngrid_declination=ngrid_dec;

    DOUBLE dec_binsize=dec_diff/ngrid_dec;
    XASSERT(NLATMAX >= (2*phibin_refine_factor + 1),
            "NLATMAX = %d needs to be larger than the minimum required number of ra cells = %d\n",
            NLATMAX, 2*phibin_refine_factor + 1);

    *ngrid_phi = my_malloc(sizeof(*ngrid_ra),ngrid_dec);
    ngrid_ra = *ngrid_phi;

    const DOUBLE costhetamax=COSD(thetamax);
    const DOUBLE max_phi_cell = 120.0;
    int max_nmesh_phi=0;
    for(int idec=0;idec<ngrid_dec;idec++) {
        DOUBLE this_min_dec;
        DOUBLE this_dec = dec_min + idec*dec_binsize;
        if(this_dec > 0) {
            int max_idec = idec + rbin_refine_factor >= ngrid_dec ? ngrid_dec-1:idec+rbin_refine_factor;
            this_min_dec = dec_min + (max_idec+1)*dec_binsize;//upper limit for that dec-bin
        } else {
            int max_idec = idec - rbin_refine_factor < 0 ? 0:idec-rbin_refine_factor;
            this_min_dec = dec_min + max_idec*dec_binsize;//lower limit for that dec-bin
        }

        DOUBLE phi_cell = max_phi_cell;
        /* if(!(i==0 || i == 1 || i == ngrid_dec-2 || i == ngrid_dec-1)) { */
        if( (90.0 - ABS(this_min_dec) ) > 1.0) { //make sure min_dec is not close to the pole (within 1 degree)-> divide by zero happens the cosine term
            DOUBLE sin_min_dec = SIND(this_min_dec),cos_min_dec=COSD(this_min_dec);
            phi_cell = ACOS((costhetamax - sin_min_dec*sin_min_dec)/(cos_min_dec*cos_min_dec))*INV_PI_OVER_180;
            /* phi_cell *= rbin_refine_factor;//My logic does not work - but multiplying with rbin_refine_factor sorts out the problem */
            if(!(phi_cell > 0.0)) {
                /* DOUBLE tmp3 = (costhetamax - sin_min_dec*sin_min_dec)/(cos_min_dec*cos_min_dec); */
                /* fprintf(stderr,"ERROR: this_min_dec = %20.16lf phi_cell = %lf is negative. thetamax = %lf sin_min_dec = %lf cos_min_dec = %lf tmp3 = %lf \n",this_min_dec,phi_cell,thetamax,sin_min_dec,cos_min_dec,tmp3); */
                phi_cell = max_phi_cell;
            }
        }
        XASSERT(phi_cell > 0.0, "RA binsize=%"DOUBLE_FORMAT" in declination bin = %d must be positive\n", phi_cell, idec);
        phi_cell = phi_cell > max_phi_cell ? max_phi_cell:phi_cell;
        const DOUBLE this_nmesh_ra = phi_diff*phibin_refine_factor/phi_cell;
        int nmesh_ra = this_nmesh_ra > NLATMAX ? this_nmesh_ra:(int) this_nmesh_ra;
        if(nmesh_ra < (2*phibin_refine_factor + 1)) {
            nmesh_ra = 2*phibin_refine_factor + 1;
            fprintf(stderr,"%s> Using sub-optimal RA binning to ensure correct functioning of the code\n",__FUNCTION__);
        }
        /* fprintf(stderr,"idec = %d nmesh_ra = %d max_nmesh_phi = %d thetamax = %lf phi_diff = %lf phi_cell = %lf phi_cell/thetamax=%lf\n",idec,nmesh_ra,max_nmesh_phi,thetamax,phi_diff,phi_cell,phi_cell/thetamax); */
        if(nmesh_ra > max_nmesh_phi) max_nmesh_phi = nmesh_ra;
        ngrid_ra[idec] = nmesh_ra;
    }

    expected_n=(int)( (np/(DOUBLE) (ngrid_dec*max_nmesh_phi)) *MEMORY_INCREASE_FAC);
    expected_n = expected_n < NVEC ? NVEC:expected_n;
    while(expected_n % NVEC != 0) {
        expected_n++;
    }
    totnbytes += ngrid_dec*max_nmesh_phi*sizeof(cellarray);


    /*---Allocate-and-initialize-grid-arrays----------*/
    lattice = (cellarray **) matrix_malloc(sizeof(cellarray),ngrid_dec,max_nmesh_phi);
    for(int idec=0;idec<ngrid_dec;idec++) {
        const int nmesh_ra = ngrid_ra[idec];
        for(int ira=0;ira<nmesh_ra;ira++) {
            const size_t memsize=3*sizeof(DOUBLE);
            lattice[idec][ira].pos     = my_malloc(memsize,expected_n);
            lattice[idec][ira].nelements=0;
            lattice[idec][ira].nallocated=expected_n;
            totnbytes += memsize*expected_n;
        }
    }

    max_n = 0;
    /*---Loop-over-particles-and-build-grid-arrays----*/
    for(int i=0;i<np;i++) {
        int idec = (int)(ngrid_dec*(dec[i]-dec_min)*inv_dec_diff);
        if(idec >= ngrid_dec) idec--;
        
        XASSERT(idec >=0 && idec < ngrid_dec,
                "Declination index for particle position = %d must be within [0, %d)\n",
                idec, ngrid_dec);

        int ira  = (int)(ngrid_ra[idec]*(phi[i]-phi_min)*inv_phi_diff);
        if(ira >=ngrid_ra[idec]) ira--;
        XASSERT(ira >=0 && ira < ngrid_ra[idec],
                "RA index for particle position = %d must be within [0, %d) for declination bin = %d\n",
                idec, ngrid_ra[idec], idec);

        if(lattice[idec][ira].nelements == lattice[idec][ira].nallocated) {
            expected_n = lattice[idec][ira].nallocated*MEMORY_INCREASE_FAC;
            while(expected_n <= lattice[idec][ira].nelements || expected_n % NVEC != 0) {
                expected_n++;
            }

            const size_t memsize=3*sizeof(DOUBLE);
            lattice[idec][ira].pos = my_realloc(lattice[idec][ira].pos ,memsize,expected_n,"lattice.pos");
            lattice[idec][ira].nallocated = expected_n;
        }
        XASSERT(lattice[idec][ira].nelements < lattice[idec][ira].nallocated,
                ANSI_COLOR_RED"BUG: lattice[%d][%d].nelements = %d must be less than allocated memory = %d" ANSI_COLOR_RESET"\n",
                idec, ira, lattice[idec][ira].nelements, lattice[idec][ira].nallocated);
        
        const int num_nvec_bunch = lattice[idec][ira].nelements/NVEC;
        const size_t xoffset = num_nvec_bunch * NVEC * 3;
        const size_t yoffset = xoffset + NVEC;
        const size_t zoffset = xoffset + 2*NVEC;
        const int ipos=lattice[idec][ira].nelements % NVEC;
        DOUBLE *xpos = &(lattice[idec][ira].pos[xoffset]);
        DOUBLE *ypos = &(lattice[idec][ira].pos[yoffset]);
        DOUBLE *zpos = &(lattice[idec][ira].pos[zoffset]);

        xpos[ipos]  = x1[i];
        ypos[ipos]  = y1[i];
        zpos[ipos]  = z1[i];
        lattice[idec][ira].nelements++;
        if(lattice[idec][ira].nelements > max_n)
            max_n = lattice[idec][ira].nelements;
    }
    *max_in_cell = max_n;
#ifndef SILENT
    gettimeofday(&t1,NULL);
    fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d ngrid_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__, totnbytes/(1024*1024.),expected_n,ngrid_dec,np,
            ADD_DIFF_TIME(t0,t1));
#endif
    return lattice;
}



#endif//LINK_IN_RA
#endif//LINK_IN_DEC

//The following is copy-pasted from the theory-side gridlink.c (used by ../xi_mocks/vpf/countspheres_mocks.c for computing the VPF on mocks)
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
    DOUBLE xdiff,ydiff,zdiff;
    DOUBLE cell_volume,box_volume;
    DOUBLE xbinsize,ybinsize,zbinsize;
    int expected_n=0;
    int64_t totncells;
    size_t totnbytes=0;

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
    expected_n=(int)(np*cell_volume/box_volume*MEMORY_INCREASE_FAC);
    expected_n=expected_n < NVEC ? NVEC:expected_n;
    while((expected_n % NVEC) != 0)
        expected_n++;

    lattice    = (cellarray *) my_malloc(sizeof(cellarray), totncells);
    int *nallocated = (int *) my_malloc(sizeof(*nallocated), totncells);

    totnbytes += sizeof(cellarray)*totncells;
    totnbytes += 3*sizeof(DOUBLE)*totncells;
    /*
      Allocate memory for each of the fields in cellarray. Since we haven't processed the data yet,
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
                ANSI_COLOR_RED"BUG: lattice[%"PRId64"].nelements = %d must be less than allocated memory = %d" ANSI_COLOR_RESET"\n",
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
    fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice. [nmesh_x, nmesh_y, nmesh_z]  = %d,%d,%d. Time taken = %6.2lf sec\n",__FUNCTION__,totnbytes/(1024*1024.0),nmesh_x,nmesh_y,nmesh_z,ADD_DIFF_TIME(t0,t1));
#endif
    return lattice;
}
