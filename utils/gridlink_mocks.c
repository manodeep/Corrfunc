/* PROGRAM gridlink1D

   --- gridlink np rmin rmax rcell z &ngrid &gridinit &gridlist
   --- Creates a 1D grid and places particles into it via a linked
   --- list.  Similar to gridlink.c, but in 1D.

   ---inputs---
      * np = number of particles
      * rmin,rmax = particles are located in a box running from 
                    (rmin,rmin,rmin) to (rmax,rmax,rmax).
      * rcell = size of a single grid cell 
      * z = array of particle coordinate that determines grid
   ---outputs---
      * ngrid = dimension of grid - computed from rmin,rmax,rcell
      * grid = 1D grid array where each cell contains the index of the 
                  first particle in that cell.
      * gridlist = array of length np containing linked list.
   -------------------------------------------------------------
      If cell (iz) contains N particles with indices j1,j2,j3,...,jN, 
      then: j1 = grid[iz], j2 = gridlist[j1], j3 = gridlist[j2],...,
      jN = gridlist[j<N-1>], and gridlist[jN] = -1.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "gridlink_mocks.h"
#include "utils.h"
#include "function_precision.h"
#include "cellarray_mocks.h"


#define MEMORY_INCREASE_FAC   1.1


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
					  )
{
  DOUBLE czmin = *min_cz, czmax=*max_cz;
#ifdef LINK_IN_DEC
  DOUBLE dec_min = *min_dec, dec_max = *max_dec;
#endif

#ifdef LINK_IN_RA
  DOUBLE ra_min = *min_ra, ra_max = *max_ra;
#endif
        
  for(int64_t i=0;i<ND1;i++) {
    if(cz[i] < czmin) czmin=cz[i];
    if(cz[i] > czmax) czmax=cz[i];

#ifdef LINK_IN_DEC
    if(dec[i] < dec_min) dec_min=dec[i];
    if(dec[i] > dec_max) dec_max=dec[i];
#endif

#ifdef LINK_IN_RA
    if(ra[i] < ra_min) ra_min=ra[i];
    if(ra[i] > ra_max) ra_max=ra[i];
#endif

  }
  *min_cz=czmin;
  *max_cz=czmax;

#ifdef LINK_IN_DEC
  *min_dec = dec_min;
  *max_dec = dec_max;
#endif

#ifdef LINK_IN_RA
  *min_ra = ra_min;
  *max_ra = ra_max;
#endif  

}



cellarray_mocks *gridlink1D(const int64_t np,const DOUBLE czmin,const DOUBLE czmax, const DOUBLE rcell,
						   const DOUBLE *theta, const DOUBLE *phi, const DOUBLE *cz, 
						   int *ngrid,int *max_in_cell,
						   const int zbin_refine_factor)
{
  int nmesh,iz,max_n;
  double sdiff = czmax-czmin;
  assert(sdiff > 0.0);
  double inv_sdiff=1.0/sdiff;
  nmesh = (int)(zbin_refine_factor*sdiff/rcell) ;
  if(nmesh>NLATMAX) nmesh=NLATMAX ;
  *ngrid=nmesh ;
  const int64_t totncells = nmesh;

  int64_t expected_n=(int64_t)((np/(double) nmesh)  *MEMORY_INCREASE_FAC);
  expected_n=expected_n <= 1 ? 2:expected_n;
#ifndef SILENT
  fprintf(stderr,"%s> Allocating %0.2g (MB) memory for the lattice, expected_n = %"PRId64" nmesh = %d np=%"PRId64" \n",__FUNCTION__,(3*4)*expected_n*nmesh/(1024.*1024.),expected_n,nmesh,np);
#endif
  cellarray_mocks *lattice = my_malloc(sizeof(cellarray_mocks), totncells);

  /*
	Allocate memory for each of the fields in cellarray. Since we haven't processed the data yet, 
	expected_n is a reasonable guess as to the number of points in the cell. 
  */
  for (int64_t index=0;index<totncells;index++) {
	lattice[index].x  = my_malloc(sizeof(DOUBLE),expected_n);//This allocates extra and is wasteful
	lattice[index].y  = my_malloc(sizeof(DOUBLE),expected_n);
	lattice[index].z  = my_malloc(sizeof(DOUBLE),expected_n);
	lattice[index].cz = my_malloc(sizeof(DOUBLE),expected_n);
	//allocate new fields in cellarray here (if you are creating a custom correlation function)
	lattice[index].nelements=0;
	lattice[index].nallocated = expected_n;
  }

  max_n=0;
  /*---Loop-over-particles-and-build-grid-arrays----*/
  for(int64_t i=0;i<np;i++) {
    iz = (int)(nmesh*(cz[i]-czmin)*inv_sdiff) ;
    if (iz >= nmesh) iz--;
    assert(iz >= 0 && iz < nmesh && "cz is inside bounds");
	
	const int64_t index = iz;
    if(lattice[index].nelements == lattice[index].nallocated) {
      expected_n = lattice[index].nallocated*MEMORY_INCREASE_FAC;

	  //In case expected_n is 1 or MEMORY_INCREASE_FAC is 1. 
	  //This way, we only increase by a very few particles 
	  // at a time. Smaller memory footprint
      if(expected_n == lattice[index].nallocated)
		expected_n += 3;

      lattice[index].x  = my_realloc(lattice[index].x ,sizeof(DOUBLE),expected_n,"lattice.x");
      lattice[index].y  = my_realloc(lattice[index].y ,sizeof(DOUBLE),expected_n,"lattice.y");
      lattice[index].z  = my_realloc(lattice[index].z ,sizeof(DOUBLE),expected_n,"lattice.z");
      lattice[index].cz = my_realloc(lattice[index].cz ,sizeof(DOUBLE),expected_n,"lattice.cz");      
      lattice[index].nallocated = expected_n;
    }
    assert(lattice[index].nelements < lattice[index].nallocated && "Ensuring that number of particles in a cell doesn't corrupt memory");
	//Index is the 1-D index for the 3-D cell. 
	//ipos is the ipos'th particle in that 3-D cell.
    int64_t ipos=lattice[index].nelements;
    lattice[index].x[ipos]  = cz[i]*COS(theta[i]*PI_OVER_180)*COS(phi[i]*PI_OVER_180);
    lattice[index].y[ipos]  = cz[i]*COS(theta[i]*PI_OVER_180)*SIN(phi[i]*PI_OVER_180) ;
    lattice[index].z[ipos]  = cz[i]*SIN(theta[i]*PI_OVER_180);
    lattice[index].cz[ipos] = cz[i];
    lattice[index].nelements++;
    if(lattice[index].nelements > max_n) {
	  max_n = lattice[index].nelements;
	}
  }
  *max_in_cell = max_n;
  return lattice;
}

#ifdef LINK_IN_DEC

cellarray_mocks **gridlink2D(const int64_t np,
						   const DOUBLE czmin, const DOUBLE czmax, const DOUBLE rcell,
						   const DOUBLE dec_min,const DOUBLE dec_max,const DOUBLE rpmax,
						   const DOUBLE *x1,const DOUBLE *y1,const DOUBLE *z1, const DOUBLE *cz,const DOUBLE *dec,
						   int *ngrid_cz,
						   int **ngrid_declination,
						   int *max_in_cell,
						   const int rbin_refine_factor,
						   const int zbin_refine_factor)

{
  int nmesh_cz,iz ;
  const DOUBLE dcz = czmax-czmin;
  /* DOUBLE inv_dcz = 1.0/dcz; */
  cellarray_mocks *tmp;

  int expected_n,index,max_n;
  size_t totnbytes=0;
  int *ngrid_dec = NULL;
  
  const DOUBLE dec_diff = dec_max-dec_min;
  DOUBLE inv_dec_diff = 1.0/dec_diff;
  DOUBLE cz_binsize,inv_cz_binsize;
  DOUBLE dec_cell=0.0,d2min=0.0;
  int nmesh_dec,idec,max_nmesh_dec;
  cellarray_mocks **lattice=NULL;
  int assigned_n=0;
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
  
  assert(dcz > 0.0);
  assert(rcell > 0.0);
  assert(dec_diff > 0.0);

  assert(MEMORY_INCREASE_FAC >= 1.0);

  nmesh_cz = (int)(dcz*zbin_refine_factor/rcell) ;
  if(nmesh_cz>NLATMAX) nmesh_cz=NLATMAX ;
  *ngrid_cz=nmesh_cz ;
  cz_binsize = dcz/nmesh_cz;
  inv_cz_binsize = 1.0/cz_binsize;
  fprintf(stderr,"nmesh_cz = %d\n",nmesh_cz);


  *ngrid_declination = my_malloc(sizeof(*ngrid_dec),nmesh_cz);
  ngrid_dec = *ngrid_declination;

  /* Find the max. number of declination cells that can be */
  DOUBLE min_dec_cell  = asin(rpmax/(2*czmax))*2.0*INV_PI_OVER_180;
  max_nmesh_dec = (int)(dec_diff*rbin_refine_factor/min_dec_cell) ;
  if(max_nmesh_dec > NLATMAX) max_nmesh_dec = NLATMAX;

  expected_n=(int)( (np/(DOUBLE) (nmesh_cz*max_nmesh_dec)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n <=10 ? 10:expected_n;
  totnbytes += nmesh_cz*max_nmesh_dec*sizeof(cellarray_mocks);
  
  /*---Allocate-and-initialize-grid-arrays----------*/
  lattice = (cellarray_mocks **) matrix_malloc(sizeof(cellarray_mocks),nmesh_cz,max_nmesh_dec); //This allocates extra and is wasteful
  for(int i=0;i<nmesh_cz;i++) {
    {
      int min_iz = (i - zbin_refine_factor) < 0  ? 0:i-zbin_refine_factor;
      d2min = czmin + 0.5*(min_iz+i)*cz_binsize;
    }
    dec_cell = asin(rpmax/(2*d2min))*2.0*INV_PI_OVER_180;
    assert(dec_cell > 0.0);
    nmesh_dec = (int)(dec_diff*rbin_refine_factor/dec_cell) ;
    if(nmesh_dec>NLATMAX)
      nmesh_dec=NLATMAX ;
    if( !(nmesh_dec > 0 && nmesh_dec <= max_nmesh_dec)) {
      fprintf(stderr,"ERROR: dec_cell = %lf czmax=%lf d2min = %lf nmesh_dec = %d max_nmesh_dec = %d\n",dec_cell,czmax,d2min,nmesh_dec,max_nmesh_dec);
    }
    assert(nmesh_dec > 0 && nmesh_dec <= max_nmesh_dec && "Number of declination cells within bounds");
    ngrid_dec[i]=nmesh_dec ;
    for(int j=0;j<nmesh_dec;j++) {
      tmp = &(lattice[i][j]);
      tmp->x     = my_malloc(sizeof(*(tmp->x)),expected_n);
      tmp->y     = my_malloc(sizeof(*(tmp->y)),expected_n);
      tmp->z     = my_malloc(sizeof(*(tmp->z)),expected_n);
      tmp->cz    = my_malloc(sizeof(*(tmp->cz)),expected_n);
      tmp->nelements=0;
      tmp->nallocated=expected_n;
      totnbytes += (sizeof(*(tmp->x)) + sizeof(*(tmp->y)) + sizeof(*(tmp->z)) )*expected_n;
    }
    /* fprintf(stderr,"ngrid_dec[%d] = %d nmesh_dec = %d\n",i,ngrid_dec[i],nmesh_dec); */
  }
  
  max_n = 0;
  /*---Loop-over-particles-and-build-grid-arrays----*/
  for(int i=0;i<np;i++) {
    iz = (int)((cz[i]-czmin)*inv_cz_binsize) ;
    if (iz >= nmesh_cz) iz--;
    assert(iz >=0 && iz < nmesh_cz && "cz position is within bounds");
    assert(dec[i] >= dec_min && dec[i] <= dec_max && "Declination within bounds");
    idec = (int)(ngrid_dec[iz]*(dec[i]-dec_min)*inv_dec_diff);
    if(idec >= ngrid_dec[iz]) idec--;
    assert(idec >=0 && idec < ngrid_dec[iz] && "Declination index within range");
    tmp = &(lattice[iz][idec]);
    if(tmp->nelements == tmp->nallocated) {
      expected_n = tmp->nallocated*MEMORY_INCREASE_FAC;
      tmp->x  = my_realloc(tmp->x ,sizeof(*(tmp->x)),expected_n,"lattice.x");
      tmp->y  = my_realloc(tmp->y ,sizeof(*(tmp->y)),expected_n,"lattice.y");
      tmp->z  = my_realloc(tmp->z ,sizeof(*(tmp->z)),expected_n,"lattice.z");
      tmp->cz = my_realloc(tmp->cz ,sizeof(*(tmp->cz)),expected_n,"lattice.cz");
      tmp->nallocated = expected_n;
    }
    index=tmp->nelements;
    tmp->x[index]  = x1[i];
    tmp->y[index]  = y1[i];
    tmp->z[index]  = z1[i];
    tmp->cz[index] = cz[i];
    tmp->nelements++;
    if(tmp->nelements > max_n)
      max_n = tmp->nelements;
    assigned_n++;
  }
  *max_in_cell = max_n;
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d nmesh_cz = %d max_nmesh_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,nmesh_cz,max_nmesh_dec,np,
	  ADD_DIFF_TIME(t0,t1));
  /* fprintf(stderr,"np = %d assigned_n = %d\n",np,assigned_n); */
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
  int ngrid_dec = 0;
  
  const DOUBLE dec_diff = dec_max-dec_min;
  DOUBLE inv_dec_diff = 1.0/dec_diff;
  DOUBLE dec_cell;
  /* int idec; */
  cellarray *lattice=NULL;
  int assigned_n=0;
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
  
  assert(thetamax > 0.0);
  assert(dec_diff > 0.0);
  assert(MEMORY_INCREASE_FAC >= 1.0);


  /* Find the max. number of declination cells that can be */
  dec_cell  = thetamax; 
  ngrid_dec = (int)(dec_diff*rbin_refine_factor/dec_cell) ;
  if(ngrid_dec >= NLATMAX)
    ngrid_dec = NLATMAX;
  
  *ngrid_declination=ngrid_dec;

  expected_n=(int)( (np/(DOUBLE) (ngrid_dec)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n <=1 ? 2:expected_n;
  totnbytes += ngrid_dec*sizeof(cellarray);
  
  /*---Allocate-and-initialize-grid-arrays----------*/
  lattice = (cellarray *) my_malloc(sizeof(cellarray),ngrid_dec); //This allocates extra and is wasteful
  for(int j=0;j<ngrid_dec;j++) {
    lattice[j].x     = my_malloc(sizeof(DOUBLE),expected_n);
    lattice[j].y     = my_malloc(sizeof(DOUBLE),expected_n);
    lattice[j].z     = my_malloc(sizeof(DOUBLE),expected_n);
    /* lattice[j].cz    = my_malloc(sizeof(*(tmp->cz)),expected_n); */
    lattice[j].nelements=0;
    lattice[j].nallocated=expected_n;
    totnbytes += 3*sizeof(DOUBLE)*expected_n;
  }

  
  max_n = 0;
  /*---Loop-over-particles-and-build-grid-arrays----*/
  for(int i=0;i<np;i++) {
    int idec = (int)(ngrid_dec*(dec[i]-dec_min)*inv_dec_diff);
    if(idec >=ngrid_dec) idec--;
    assert(idec >=0 && idec < ngrid_dec && "Declination is within bounds");
    if(lattice[idec].nelements == lattice[idec].nallocated) {
      expected_n = lattice[idec].nallocated*MEMORY_INCREASE_FAC;
			while(expected_n <= lattice[idec].nelements)
				expected_n +=5;
			
      lattice[idec].x  = my_realloc(lattice[idec].x ,sizeof(DOUBLE),expected_n,"lattice.x");
      lattice[idec].y  = my_realloc(lattice[idec].y ,sizeof(DOUBLE),expected_n,"lattice.y");
      lattice[idec].z  = my_realloc(lattice[idec].z ,sizeof(DOUBLE),expected_n,"lattice.z");
      lattice[idec].nallocated = expected_n;
    }
    const int ipos=lattice[idec].nelements;
		if( ! (lattice[idec].nallocated > ipos ) ) {
			fprintf(stderr,"ERROR: About to crash. i = %d idec = %d nelements = %d  nallocated = %d expected_n = %d \n",i,idec,lattice[idec].nelements,lattice[idec].nallocated,expected_n);
		}
		assert(lattice[idec].nallocated > ipos && "Enough memory has been allocated to assign particles");
    lattice[idec].x[ipos]  = x1[i];
    lattice[idec].y[ipos]  = y1[i];
    lattice[idec].z[ipos]  = z1[i];
    lattice[idec].nelements++;
    if(lattice[idec].nelements > max_n)
      max_n = lattice[idec].nelements;
    assigned_n++;
  }
  *max_in_cell = max_n;
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d ngrid_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,ngrid_dec,np,
	  ADD_DIFF_TIME(t0,t1));
  /* fprintf(stderr,"np = %d assigned_n = %d\n",np,assigned_n); */
  return lattice;
}



#ifdef LINK_IN_RA

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
							   const int phibin_refine_factor)

{
  int nmesh_cz;
  const DOUBLE dcz = czmax-czmin;
  cellarray_mocks *tmp;
  int expected_n,index,max_n;
  size_t totnbytes=0;
  int *ngrid_dec = NULL;
  
  const DOUBLE dec_diff = dec_max-dec_min;
  DOUBLE inv_dec_diff = 1.0/dec_diff;
  DOUBLE cz_binsize,inv_cz_binsize;
  DOUBLE dec_cell=0.0,d2min=0.0;
  int nmesh_dec,max_nmesh_dec;
  const DOUBLE phi_diff = phi_max - phi_min;
  DOUBLE inv_phi_diff = 1.0/phi_diff;
  DOUBLE phi_cell=0.0;
  DOUBLE *dec_binsizes=NULL;
  int **ngrid_ra=NULL;
  
  cellarray_mocks ***lattice=NULL;
  int assigned_n=0;
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
  
  assert(dcz > 0.0);
  assert(rcell > 0.0);
  assert(dec_diff > 0.0);

  assert(MEMORY_INCREASE_FAC >= 1.0);

  nmesh_cz = (int)(dcz*zbin_refine_factor/rcell) ;
  if(nmesh_cz>NLATMAX) nmesh_cz=NLATMAX ;
  *ngrid_cz=nmesh_cz ;
  cz_binsize = dcz/nmesh_cz;
  inv_cz_binsize = 1.0/cz_binsize;

  *ngrid_declination = my_malloc(sizeof(*ngrid_dec),nmesh_cz);
  ngrid_dec = *ngrid_declination;

  /* Find the max. number of declination cells that can be */
  DOUBLE min_dec_cell  = asin(rpmax/(2*czmax))*2.0*INV_PI_OVER_180;
  max_nmesh_dec = (int)(dec_diff*rbin_refine_factor/min_dec_cell) ;
  if(max_nmesh_dec > NLATMAX) max_nmesh_dec = NLATMAX;
  DOUBLE thetamax=dec_diff/max_nmesh_dec;

  dec_binsizes=my_malloc(sizeof(*dec_binsizes),nmesh_cz);
  DOUBLE min_phi_cell = thetamax;
  int max_nmesh_phi = (int) (phi_diff*phibin_refine_factor/min_phi_cell) ;
  if(max_nmesh_phi > NLATMAX) max_nmesh_phi = NLATMAX;
  
  expected_n=(int)( (np/(DOUBLE) (nmesh_cz*max_nmesh_dec*max_nmesh_phi)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n <=10 ? 10:expected_n;
  totnbytes += nmesh_cz*max_nmesh_dec*max_nmesh_phi*sizeof(cellarray_mocks);

  /*---Allocate-and-initialize-grid-arrays----------*/
  *ngrid_phi = (int **) matrix_malloc(sizeof(int), nmesh_cz, max_nmesh_dec);
  ngrid_ra = *ngrid_phi;
  lattice = (cellarray_mocks ***) volume_malloc(sizeof(cellarray_mocks),nmesh_cz,max_nmesh_dec,max_nmesh_phi); //This allocates extra and is wasteful
  DOUBLE costhetamax=COSD(thetamax);
  for(int iz=0;iz<nmesh_cz;iz++) {
    {
      int min_iz = iz-zbin_refine_factor < 0 ? 0:iz-zbin_refine_factor;
      d2min = czmin + 0.5*(min_iz+iz)*cz_binsize;
    }
    /* d2min = czmin + iz*cz_binsize; */
    dec_cell = asin(rpmax/(2*d2min))*2.0*INV_PI_OVER_180;// \sigma = 2*arcsin(C/2) -> 2*arcsin( (rpmax/d2min) /2)
    assert(dec_cell > 0.0);
    dec_binsizes[iz] = dec_cell;
    nmesh_dec = (int)(dec_diff*rbin_refine_factor/dec_cell) ;
    if(nmesh_dec > NLATMAX) nmesh_dec = NLATMAX;
    assert(nmesh_dec <= max_nmesh_dec && "Number of Declination cells is within bounds");
    ngrid_dec[iz]=nmesh_dec ;
  }
  for(int iz=0;iz<nmesh_cz;iz++) {
    nmesh_dec = ngrid_dec[iz];
    DOUBLE dec_binsize=dec_diff/nmesh_dec;
    int min_iz = iz-zbin_refine_factor < 0 ? 0:iz-zbin_refine_factor;
    /* DOUBLE max_dec_binsize = dec_diff/ngrid_dec[min_iz]; */
    /* dec_cell = rbin_refine_factor*dec_diff/nmesh_dec; */
    /* fprintf(stderr,"ngrid_dec[%03d] = %03d dec_cell = %lf \n",iz,nmesh_dec,dec_cell); */
    /* costhetamax=COSD(max_dec_binsize); */
    dec_cell = dec_binsizes[min_iz];
    costhetamax=COSD(dec_cell);
    for(int idec=0;idec<nmesh_dec;idec++) {
      int nmesh_ra,max_idec;
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
      /* 	      ,min_iz,idec,max_idec,nmesh_dec,ngrid_dec[min_iz],costhetamax,COSD(dec_binsize)); */

      phi_cell = 120.0;
      /* if(!(max_idec==0 || max_idec==1 || max_idec == (nmesh_dec-2) ||max_idec == (nmesh_dec-1))) { */
      /* if(!(max_idec==0 || max_idec == (nmesh_dec-1))) { */
      if( (90.0 - fabs(this_min_dec) ) > 1.0) { //make sure min_dec is not close to the pole (within 1 degree)
	DOUBLE tmp1 = SIND(this_min_dec),tmp2=COSD(this_min_dec);
	phi_cell = acos((costhetamax - tmp1*tmp1)/(tmp2*tmp2))*INV_PI_OVER_180;
	/* phi_cell *= rbin_refine_factor;//My logic does not work - but multiplying with rbin_refine_factor sorts out the problem */
	/* phi_cell *= 1.2;//Still needs a fudge-factor */
	if(!(phi_cell > 0.0)) {
	  /* DOUBLE tmp3 = (costhetamax - tmp1*tmp1)/(tmp2*tmp2); */
	  /* fprintf(stderr,"ERROR: idec = %d max_idec = %d nmesh_dec = %d this_min_dec = %lf dec_cell = %lf phi_cell = %lf is negative. thetamax = %lf tmp1 = %lf tmp2 = %lf tmp3 = %lf \n", */
	  /* 	  idec,max_idec,nmesh_dec,this_min_dec,dec_cell,phi_cell,dec_cell,tmp1,tmp2,tmp3); */
	  phi_cell = 120.0;
	}
      }
      assert(phi_cell > 0.0 && "RA bin-width is positive");
      phi_cell = phi_cell > 120.0 ? 120.0:phi_cell;
      /* fprintf(stderr,"iz = %4d idec = %4d dec_cell = %6.3lf dec_binsize=%6.2lf this_dec = %6.2lf phi_cell = %7.2lf \n",iz,idec,dec_cell,dec_binsize,dec_min + idec*dec_binsize,phi_cell); */
      nmesh_ra = (int) (phi_diff*phibin_refine_factor/phi_cell);
      if(nmesh_ra > NLATMAX)
	nmesh_ra = NLATMAX;
      assert(nmesh_ra <= max_nmesh_phi && "Number of RA cells in within bounds");
      ngrid_ra[iz][idec] = nmesh_ra;
      for(int ira=0;ira<nmesh_ra;ira++) {
	tmp = &(lattice[iz][idec][ira]);
	tmp->x     = my_malloc(sizeof(DOUBLE),expected_n);
	tmp->y     = my_malloc(sizeof(DOUBLE),expected_n);
	tmp->z     = my_malloc(sizeof(DOUBLE),expected_n);
	tmp->cz    = my_malloc(sizeof(DOUBLE),expected_n);
	tmp->nelements=0;
	tmp->nallocated=expected_n;
	totnbytes += 4*sizeof(DOUBLE)*expected_n;
      }
    }
    /* fprintf(stderr,"ngrid_dec[%d] = %d nmesh_dec = %d\n",i,ngrid_dec[i],nmesh_dec); */
  }
  
  max_n = 0;
  /*---Loop-over-particles-and-build-grid-arrays----*/
  for(int i=0;i<np;i++) {
    if(cz[i] >=czmin && cz[i] <= czmax) {
      int iz = (int)((cz[i]-czmin)*inv_cz_binsize) ;
      if (iz >= nmesh_cz) iz--;
      assert(iz >=0 && iz < nmesh_cz && "cz (particle) position is within bounds");
      /* assert(dec[i] >= dec_min && dec[i] <= dec_max); */
      int idec = (int)(ngrid_dec[iz]*(dec[i]-dec_min)*inv_dec_diff);
      if(idec >= ngrid_dec[iz]) idec--;
      assert(idec >=0 && idec < ngrid_dec[iz] && "Dec (particle) position within bounds");
      int ira = (int) (ngrid_ra[iz][idec]*(phi[i]-phi_min)*inv_phi_diff);
      if(ira >= ngrid_ra[iz][idec]) ira--;
      assert(ira >=0 && ira < ngrid_ra[iz][idec] && "RA (particle) position within bounds");
      tmp = &(lattice[iz][idec][ira]);
      if(tmp->nelements == tmp->nallocated) {
	expected_n = tmp->nallocated*MEMORY_INCREASE_FAC;
	tmp->x  = my_realloc(tmp->x ,sizeof(*(tmp->x)),expected_n,"lattice.x");
	tmp->y  = my_realloc(tmp->y ,sizeof(*(tmp->y)),expected_n,"lattice.y");
	tmp->z  = my_realloc(tmp->z ,sizeof(*(tmp->z)),expected_n,"lattice.z");
	tmp->cz = my_realloc(tmp->cz ,sizeof(*(tmp->cz)),expected_n,"lattice.cz");
	tmp->nallocated = expected_n;
      }
      index=tmp->nelements;
      tmp->x[index]  = x1[i];
      tmp->y[index]  = y1[i];
      tmp->z[index]  = z1[i];
      tmp->cz[index] = cz[i];
      tmp->nelements++;
      if(tmp->nelements > max_n) {
	max_n = tmp->nelements;
      }
      assigned_n++;
    }
  }
  free(dec_binsizes);
  *max_in_cell = max_n;
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d (max_n = %d) nmesh_cz = %d max_nmesh_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,max_n,nmesh_cz,max_nmesh_dec,np,
	  ADD_DIFF_TIME(t0,t1));
  /* fprintf(stderr,"np = %d assigned_n = %d\n",np,assigned_n); */
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
  int expected_n,index,max_n;
  size_t totnbytes=0;
  int ngrid_dec = 0;
  
  const DOUBLE dec_diff = dec_max-dec_min;
  DOUBLE inv_dec_diff = 1.0/dec_diff;
  DOUBLE dec_cell=0.0;
  /* int idec; */
  int *ngrid_ra=NULL;

  const DOUBLE phi_diff = phi_max - phi_min;
  DOUBLE inv_phi_diff = 1.0/phi_diff;
  DOUBLE phi_cell=0.0;
  /* int ira; */
  
  cellarray **lattice=NULL;
  int assigned_n=0;
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
  
  assert(thetamax > 0.0);
  assert(dec_diff > 0.0);
  assert(phi_diff > 0.0);
  assert(MEMORY_INCREASE_FAC >= 1.0);
  

  /* Find the max. number of declination cells that can be */
  dec_cell  = thetamax; 
  ngrid_dec = (int)(dec_diff*rbin_refine_factor/dec_cell) ;
  if(ngrid_dec > NLATMAX)
    ngrid_dec = NLATMAX;
  
  /* assert(ngrid_dec <= NLATMAX); */
  *ngrid_declination=ngrid_dec;
  DOUBLE dec_binsize=dec_diff/ngrid_dec;

	assert(NLATMAX >= (2*phibin_refine_factor + 1) && "NLATMAX needs to be larger than the minimum required number of ra cells for correct code function");
  DOUBLE min_phi_cell = thetamax;
  DOUBLE max_phi_cell = phi_diff/(2*phibin_refine_factor + 1);
  int max_nmesh_phi = (int) (phi_diff*phibin_refine_factor/min_phi_cell) ;
  if(max_nmesh_phi > NLATMAX) max_nmesh_phi = NLATMAX;
  max_nmesh_phi = max_nmesh_phi < (2*phibin_refine_factor + 1) ? (2*phibin_refine_factor+1):max_nmesh_phi;
  /* fprintf(stderr,"phi_diff = %lf thetamax = %lf min_phi_cell = %lf max_nmesh_phi = %d \n",phi_diff,thetamax,min_phi_cell,max_nmesh_phi); */
  
  *ngrid_phi = my_malloc(sizeof(*ngrid_ra),ngrid_dec);
  ngrid_ra = *ngrid_phi;

  expected_n=(int)( (np/(DOUBLE) (ngrid_dec*max_nmesh_phi)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n <=10 ? 10:expected_n;
  totnbytes += ngrid_dec*max_nmesh_phi*sizeof(cellarray);
  
  
  /*---Allocate-and-initialize-grid-arrays----------*/
  lattice = (cellarray **) matrix_malloc(sizeof(cellarray),ngrid_dec,max_nmesh_phi); 
  DOUBLE costhetamax=COSD(thetamax);
  for(int idec=0;idec<ngrid_dec;idec++) {
    /* DOUBLE this_min_dec = dec_min + i*dec_binsize; */
    DOUBLE this_min_dec;
    DOUBLE this_dec = dec_min + idec*dec_binsize;
    if(this_dec > 0) {
      int max_idec = idec + rbin_refine_factor >= ngrid_dec ? ngrid_dec-1:idec+rbin_refine_factor;
      this_min_dec = dec_min + (max_idec+1)*dec_binsize;//upper limit for that dec-bin
    } else {
      int max_idec = idec - rbin_refine_factor < 0 ? 0:idec-rbin_refine_factor;
      this_min_dec = dec_min + max_idec*dec_binsize;//lower limit for that dec-bin
    }

    phi_cell = max_phi_cell;
    /* if(!(i==0 || i == 1 || i == ngrid_dec-2 || i == ngrid_dec-1)) { */
    if( (90.0 - fabs(this_min_dec) ) > 1.0) { //make sure min_dec is not close to the pole (within 1 degree)      
      DOUBLE tmp1 = SIND(this_min_dec),tmp2=COSD(this_min_dec);
      phi_cell = acos((costhetamax - tmp1*tmp1)/(tmp2*tmp2))*INV_PI_OVER_180;
      /* phi_cell *= rbin_refine_factor;//My logic does not work - but multiplying with rbin_refine_factor sorts out the problem */
      if(!(phi_cell > 0.0)) {
		/* DOUBLE tmp3 = (costhetamax - tmp1*tmp1)/(tmp2*tmp2); */
		/* fprintf(stderr,"ERROR: this_min_dec = %20.16lf phi_cell = %lf is negative. thetamax = %lf tmp1 = %lf tmp2 = %lf tmp3 = %lf \n",this_min_dec,phi_cell,thetamax,tmp1,tmp2,tmp3); */
				phi_cell = max_phi_cell;
      }
    }
    assert(phi_cell > 0.0);
    phi_cell = phi_cell > max_phi_cell ? max_phi_cell:phi_cell;
		
    int nmesh_ra = (int) (phi_diff*phibin_refine_factor/phi_cell);
    if(nmesh_ra > NLATMAX)
      nmesh_ra = NLATMAX;

    nmesh_ra = nmesh_ra < (2*phibin_refine_factor + 1) ? (2*phibin_refine_factor+1):nmesh_ra;
    /* fprintf(stderr,"idec = %d nmesh_ra = %d max_nmesh_phi = %d thetamax = %lf phi_diff = %lf phi_cell = %lf phi_cell/thetamax=%lf\n",idec,nmesh_ra,max_nmesh_phi,thetamax,phi_diff,phi_cell,phi_cell/thetamax); */
    assert(nmesh_ra <= max_nmesh_phi);

    ngrid_ra[idec] = nmesh_ra;
    
    for(int ira=0;ira<nmesh_ra;ira++) {
      lattice[idec][ira].x     = my_malloc(sizeof(DOUBLE),expected_n);
      lattice[idec][ira].y     = my_malloc(sizeof(DOUBLE),expected_n);
      lattice[idec][ira].z     = my_malloc(sizeof(DOUBLE),expected_n);
      lattice[idec][ira].nelements=0;
      lattice[idec][ira].nallocated=expected_n;
      totnbytes += 3*sizeof(DOUBLE)*expected_n;
    }
  }

  
  max_n = 0;
  /*---Loop-over-particles-and-build-grid-arrays----*/
  for(int i=0;i<np;i++) {
    int idec = (int)(ngrid_dec*(dec[i]-dec_min)*inv_dec_diff);
    if(idec >= ngrid_dec) idec--;
    assert(idec >=0 && idec < ngrid_dec);
    int ira  = (int)(ngrid_ra[idec]*(phi[i]-phi_min)*inv_phi_diff);
    if(ira >=ngrid_ra[idec]) ira--;
    assert(ira >=0 && ira < ngrid_ra[idec]);
    if(lattice[idec][ira].nelements == lattice[idec][ira].nallocated) {
      expected_n = lattice[idec][ira].nallocated*MEMORY_INCREASE_FAC;
      lattice[idec][ira].x  = my_realloc(lattice[idec][ira].x ,sizeof(DOUBLE),expected_n,"lattice.x");
      lattice[idec][ira].y  = my_realloc(lattice[idec][ira].y ,sizeof(DOUBLE),expected_n,"lattice.y");
      lattice[idec][ira].z  = my_realloc(lattice[idec][ira].z ,sizeof(DOUBLE),expected_n,"lattice.z");
      /* lattice[idec][ira].cz = my_realloc(lattice[idec][ira].cz ,sizeof(*(lattice[idec][ira].cz)),expected_n,"lattice.cz"); */
      lattice[idec][ira].nallocated = expected_n;
    }
    index=lattice[idec][ira].nelements;
    lattice[idec][ira].x[index]  = x1[i];
    lattice[idec][ira].y[index]  = y1[i];
    lattice[idec][ira].z[index]  = z1[i];
    lattice[idec][ira].nelements++;
    if(lattice[idec][ira].nelements > max_n)
      max_n = lattice[idec][ira].nelements;
    assigned_n++;
  }
  *max_in_cell = max_n;
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d ngrid_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__, totnbytes/(1024*1024.),expected_n,ngrid_dec,np,
	  ADD_DIFF_TIME(t0,t1));
  /* fprintf(stderr,"np = %d assigned_n = %d\n",np,assigned_n); */
  return lattice;
}



#endif


#endif

