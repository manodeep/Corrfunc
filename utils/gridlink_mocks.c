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



cellarray_mocks *gridlink1D(const int64_t np,const DOUBLE czmin,const DOUBLE czmax, const DOUBLE rcell,
														const DOUBLE *dec, const DOUBLE *ra, const DOUBLE *cz,
														int *ngrid,int *max_in_cell,
														const int zbin_refine_factor)
{
	int nmesh,iz,max_n;
	const DOUBLE sdiff = (czmax-czmin);
  assert(sdiff > 0.0 && "There needs to be some depth to the data");
  const DOUBLE inv_sdiff=1.0/sdiff;
  nmesh = (int)(zbin_refine_factor*sdiff/rcell) ;
  if(nmesh>NLATMAX) nmesh=NLATMAX ;
  *ngrid=nmesh ;
  const int64_t totncells = nmesh;

  int64_t expected_n=(int64_t)((np/(double) nmesh)  *MEMORY_INCREASE_FAC);
  expected_n=expected_n < NVEC ? NVEC:expected_n;
	while(expected_n % NVEC != 0) {
		expected_n++;
	}
	
#ifndef SILENT
  fprintf(stderr,"%s> Allocating %0.2g (MB) memory for the lattice, expected_n = %"PRId64" nmesh = %d np=%"PRId64" \n",__FUNCTION__,(3*4)*expected_n*nmesh/(1024.*1024.),expected_n,nmesh,np);
#endif
  cellarray_mocks *lattice = my_malloc(sizeof(cellarray_mocks), totncells);

  /*
	Allocate memory for each of the fields in cellarray. Since we haven't processed the data yet,
	expected_n is a reasonable guess as to the number of points in the cell.
  */
  for (int64_t index=0;index<totncells;index++) {
		const size_t memsize=4*sizeof(DOUBLE);
		lattice[index].pos = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
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
			while(expected_n <= lattice[index].nelements || expected_n % NVEC != 0) {
				expected_n++;
			}
			


			const size_t memsize=4*sizeof(DOUBLE);
			lattice[index].pos  = my_realloc(lattice[index].pos ,memsize,expected_n,"lattice.pos");
			lattice[index].nallocated = expected_n;
		}
		assert(lattice[index].nallocated > lattice[index].nelements && "Making sure memory access if fine");
		const int num_nvec_bunch = lattice[index].nelements/NVEC;
		const size_t xoffset  = num_nvec_bunch * NVEC * 4;
		const size_t yoffset  = xoffset + NVEC;
		const size_t zoffset  = xoffset + 2*NVEC;
		const size_t czoffset = xoffset + 3*NVEC;			
		
		const int ipos=lattice[index].nelements % NVEC;
		
		DOUBLE *xpos  = &(lattice[index].pos[xoffset]);
		DOUBLE *ypos  = &(lattice[index].pos[yoffset]);
		DOUBLE *zpos  = &(lattice[index].pos[zoffset]);
		DOUBLE *czpos = &(lattice[index].pos[czoffset]);
		
		xpos[ipos]  = cz[i]*COSD(dec[i])*COSD(ra[i]);
		ypos[ipos]  = cz[i]*COSD(dec[i])*SIND(ra[i]);
		zpos[ipos]  = cz[i]*SIND(dec[i]);
		czpos[ipos] = cz[i];

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
														 const DOUBLE *cz,const DOUBLE *dec, const DOUBLE *ra,
														 int *ngrid_cz,
														 int **ngrid_declination,
														 int *max_in_cell,
														 const int rbin_refine_factor,
														 const int zbin_refine_factor)
	
{
  int nmesh_cz,iz ;
  const DOUBLE dcz = czmax-czmin;
  /* DOUBLE inv_dcz = 1.0/dcz; */

  int expected_n,max_n;
  size_t totnbytes=0;
  int *ngrid_dec = NULL;
  
  const DOUBLE dec_diff = dec_max-dec_min;
  DOUBLE inv_dec_diff = 1.0/dec_diff;
  DOUBLE cz_binsize,inv_cz_binsize;
  DOUBLE dec_cell=0.0,d2min=0.0;
  int nmesh_dec,idec,max_nmesh_dec;
  cellarray_mocks **lattice=NULL;
  int assigned_n=0;

#ifndef SILENT
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
#endif
	
  assert(dcz > 0.0 && "There has to be some depth to the survey");
  assert(rcell > 0.0 && "Minimum separation has to be non-zero");
  assert(dec_diff > 0.0 && "All of the points can not be at the same declination");

  assert(MEMORY_INCREASE_FAC >= 1.0 && "Memory increase factor must be >=1 ");

  nmesh_cz = (int)(dcz*zbin_refine_factor/rcell) ;
  if(nmesh_cz>NLATMAX) nmesh_cz=NLATMAX ;
  *ngrid_cz=nmesh_cz ;
  cz_binsize = dcz/nmesh_cz;
  inv_cz_binsize = 1.0/cz_binsize;
  /* fprintf(stderr,"nmesh_cz = %d\n",nmesh_cz); */


  *ngrid_declination = my_malloc(sizeof(*ngrid_dec),nmesh_cz);
  ngrid_dec = *ngrid_declination;

  /* Find the max. number of declination cells that can be */
  DOUBLE min_dec_cell  = ASIN(rpmax/(2*czmax))*2.0*INV_PI_OVER_180;
  max_nmesh_dec = (int)(dec_diff*rbin_refine_factor/min_dec_cell) ;
  if(max_nmesh_dec > NLATMAX) max_nmesh_dec = NLATMAX;

  expected_n=(int)( (np/(DOUBLE) (nmesh_cz*max_nmesh_dec)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n < NVEC ? NVEC:expected_n;
	while(expected_n % NVEC != 0) {
		expected_n++;
	}
  totnbytes += nmesh_cz*max_nmesh_dec*sizeof(cellarray_mocks);
  
  /*---Allocate-and-initialize-grid-arrays----------*/
  lattice = (cellarray_mocks **) matrix_malloc(sizeof(cellarray_mocks),nmesh_cz,max_nmesh_dec); //This allocates extra and is wasteful
  for(int i=0;i<nmesh_cz;i++) {
    {
      int min_iz = (i - zbin_refine_factor) < 0  ? 0:i-zbin_refine_factor;
      d2min = czmin + 0.5*(min_iz+i)*cz_binsize;
    }
    dec_cell = ASIN(rpmax/(2*d2min))*2.0*INV_PI_OVER_180;
    assert(dec_cell > 0.0 && "Declination binsize is non-zero");
    nmesh_dec = (int)(dec_diff*rbin_refine_factor/dec_cell) ;
    if(nmesh_dec>NLATMAX)nmesh_dec=NLATMAX ;
		
    if( !(nmesh_dec > 0 && nmesh_dec <= max_nmesh_dec)) {
      fprintf(stderr,"ERROR: dec_cell = %lf czmax=%lf d2min = %lf nmesh_dec = %d max_nmesh_dec = %d\n",dec_cell,czmax,d2min,nmesh_dec,max_nmesh_dec);
    }
    assert(nmesh_dec > 0 && nmesh_dec <= max_nmesh_dec && "Number of declination cells within bounds");
    ngrid_dec[i]=nmesh_dec ;
    for(int j=0;j<nmesh_dec;j++) {
			const size_t memsize = 4*sizeof(DOUBLE);
      lattice[i][j].pos     = my_malloc(memsize,expected_n);
      lattice[i][j].nelements=0;
      lattice[i][j].nallocated=expected_n;
      totnbytes += memsize*expected_n;
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
    if(lattice[iz][idec].nelements == lattice[iz][idec].nallocated) {
      expected_n = lattice[iz][idec].nallocated*MEMORY_INCREASE_FAC;
			while(expected_n <= lattice[iz][idec].nelements || expected_n % NVEC != 0) {
				expected_n++;
			}
			const size_t memsize = 4*sizeof(DOUBLE);
      lattice[iz][idec].pos  = my_realloc(lattice[iz][idec].pos ,memsize,expected_n,"lattice.pos");
      lattice[iz][idec].nallocated = expected_n;
    }
		assert(lattice[iz][idec].nallocated > lattice[iz][idec].nelements && "Making sure memory access is fine");
		const int num_nvec_bunch = lattice[iz][idec].nelements/NVEC;
		const size_t xoffset  = num_nvec_bunch * NVEC * 4;
		const size_t yoffset  = xoffset + NVEC;
		const size_t zoffset  = xoffset + 2*NVEC;
		const size_t czoffset = xoffset + 3*NVEC;			
		
		const int ipos=lattice[iz][idec].nelements % NVEC;
		
		DOUBLE *xpos  = &(lattice[iz][idec].pos[xoffset]);
		DOUBLE *ypos  = &(lattice[iz][idec].pos[yoffset]);
		DOUBLE *zpos  = &(lattice[iz][idec].pos[zoffset]);
		DOUBLE *czpos = &(lattice[iz][idec].pos[czoffset]);
		
		xpos[ipos]  = cz[i]*COSD(dec[i])*COSD(ra[i]);
		ypos[ipos]  = cz[i]*COSD(dec[i])*SIND(ra[i]);
		zpos[ipos]  = cz[i]*SIND(dec[i]);
		czpos[ipos] = cz[i];

		lattice[iz][idec].nelements++;
    if(lattice[iz][idec].nelements > max_n)
      max_n = lattice[iz][idec].nelements;
    assigned_n++;
  }
  *max_in_cell = max_n;
#ifndef SILENT
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d nmesh_cz = %d max_nmesh_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,nmesh_cz,max_nmesh_dec,np,
	  ADD_DIFF_TIME(t0,t1));
#endif
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
  const DOUBLE inv_dec_diff = 1.0/dec_diff;
  DOUBLE dec_cell;
  /* int idec; */
  cellarray *lattice=NULL;
  int assigned_n=0;

#ifndef SILENT
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
#endif
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
    assert(idec >=0 && idec < ngrid_dec && "Declination is within bounds");
    if(lattice[idec].nelements == lattice[idec].nallocated) {
      expected_n = lattice[idec].nallocated*MEMORY_INCREASE_FAC;
			while(expected_n <= lattice[idec].nelements || expected_n % NVEC != 0){
				expected_n++;
			}

			const size_t memsize = 3*sizeof(DOUBLE);
      lattice[idec].pos  = my_realloc(lattice[idec].pos, memsize,expected_n,"lattice.pos");
      lattice[idec].nallocated = expected_n;
    }
		assert(lattice[idec].nallocated > lattice[idec].nelements && "Enough memory has been allocated to assign particles");
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
															 const int zbin_refine_factor)
															 

{
  int nmesh_cz;
  const DOUBLE dcz = (czmax-czmin);
  int expected_n,max_n;
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
#ifndef SILENT	
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
#endif  
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
  DOUBLE min_dec_cell  = ASIN(rpmax/(2*czmax))*2.0*INV_PI_OVER_180;
  max_nmesh_dec = (int)(dec_diff*rbin_refine_factor/min_dec_cell) ;
  if(max_nmesh_dec > NLATMAX) max_nmesh_dec = NLATMAX;
  DOUBLE thetamax=dec_diff/max_nmesh_dec;

  dec_binsizes=my_malloc(sizeof(*dec_binsizes),nmesh_cz);
  DOUBLE min_phi_cell = thetamax;
  int max_nmesh_phi = (int) (phi_diff*phibin_refine_factor/min_phi_cell) ;
  if(max_nmesh_phi > NLATMAX) max_nmesh_phi = NLATMAX;
  
  expected_n=(int)( (np/(DOUBLE) (nmesh_cz*max_nmesh_dec*max_nmesh_phi)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n < NVEC ? NVEC:expected_n;

	//But we are going to store NVEC's each. 
	while((expected_n % NVEC) != 0) {
		expected_n++;
	}
	
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
    dec_cell =  ASIN(rpmax/(2*d2min))*2.0*INV_PI_OVER_180;// \sigma = 2*arcsin(C/2) -> 2*arcsin( (rpmax/d2min) /2)
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
      if( (90.0 - FABS(this_min_dec) ) > 1.0) { //make sure min_dec is not close to the pole (within 1 degree)
				DOUBLE tmp1 = SIND(this_min_dec),tmp2=COSD(this_min_dec);
				phi_cell = ACOS((costhetamax - tmp1*tmp1)/(tmp2*tmp2))*INV_PI_OVER_180;
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
				const size_t memsize=4*sizeof(DOUBLE);//4 pointers for x/y/z/cz
				lattice[iz][idec][ira].pos = my_malloc(memsize,expected_n);//This allocates extra and is wasteful
				lattice[iz][idec][ira].nelements=0;
				lattice[iz][idec][ira].nallocated=expected_n;
				totnbytes += memsize*expected_n;
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
      if(lattice[iz][idec][ira].nelements == lattice[iz][idec][ira].nallocated) {
				expected_n = lattice[iz][idec][ira].nallocated*MEMORY_INCREASE_FAC;
				while(expected_n <= lattice[iz][idec][ira].nelements || expected_n % NVEC != 0) {
					expected_n++;
				}
				
				const size_t memsize=4*sizeof(DOUBLE);
				lattice[iz][idec][ira].pos  = my_realloc(lattice[iz][idec][ira].pos ,memsize,expected_n,"lattice.pos");
				lattice[iz][idec][ira].nallocated = expected_n;
      }
			assert(lattice[iz][idec][ira].nallocated > lattice[iz][idec][ira].nelements && "Making sure memory access if fine");
			const int num_nvec_bunch = lattice[iz][idec][ira].nelements/NVEC;
			const size_t xoffset  = num_nvec_bunch * NVEC * 4;
			const size_t yoffset  = xoffset + NVEC;
			const size_t zoffset  = xoffset + 2*NVEC;
			const size_t czoffset = xoffset + 3*NVEC;			

			const int ipos=lattice[iz][idec][ira].nelements % NVEC;

			DOUBLE *xpos  = &(lattice[iz][idec][ira].pos[xoffset]);
			DOUBLE *ypos  = &(lattice[iz][idec][ira].pos[yoffset]);
			DOUBLE *zpos  = &(lattice[iz][idec][ira].pos[zoffset]);
			DOUBLE *czpos = &(lattice[iz][idec][ira].pos[czoffset]);

			xpos[ipos]  = cz[i]*COSD(dec[i])*COSD(phi[i]);
			ypos[ipos]  = cz[i]*COSD(dec[i])*SIND(phi[i]);
			zpos[ipos]  = cz[i]*SIND(dec[i]);
			czpos[ipos] = cz[i];
			
      lattice[iz][idec][ira].nelements++;
      if(lattice[iz][idec][ira].nelements > max_n) {
				max_n = lattice[iz][idec][ira].nelements;
      }
      assigned_n++;
    }
  }
  free(dec_binsizes);
  *max_in_cell = max_n;
#ifndef SILENT	
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d (max_n = %d) nmesh_cz = %d max_nmesh_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__,totnbytes/(1024*1024.),expected_n,max_n,nmesh_cz,max_nmesh_dec,np, ADD_DIFF_TIME(t0,t1));
#endif	
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
  int expected_n,max_n;
  size_t totnbytes=0;
  int ngrid_dec = 0;

	assert(thetamax > 0 && "Angular separation must be non-zero");
	
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
#ifndef SILENT
  struct timeval t0,t1;
  gettimeofday(&t0,NULL);
#endif
	
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
  /* fprintf(stderr,"phi_diff = %lf thetamax = %lf min_phi_cell = %lf max_nmesh_phi = %d 2*phi_bin_refine_factor +1 = %d\n",phi_diff,thetamax,min_phi_cell,max_nmesh_phi,(2*phibin_refine_factor+1)); */
  
  *ngrid_phi = my_malloc(sizeof(*ngrid_ra),ngrid_dec);
  ngrid_ra = *ngrid_phi;

  expected_n=(int)( (np/(DOUBLE) (ngrid_dec*max_nmesh_phi)) *MEMORY_INCREASE_FAC);
  expected_n = expected_n < NVEC ? NVEC:expected_n;
	while(expected_n % NVEC != 0) {
		expected_n++;
	}
	
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
    if( (90.0 - ABS(this_min_dec) ) > 1.0) { //make sure min_dec is not close to the pole (within 1 degree)-> divide by zero happens the cosine term
      DOUBLE tmp1 = SIND(this_min_dec),tmp2=COSD(this_min_dec);
      phi_cell = ACOS((costhetamax - tmp1*tmp1)/(tmp2*tmp2))*INV_PI_OVER_180;
      /* phi_cell *= rbin_refine_factor;//My logic does not work - but multiplying with rbin_refine_factor sorts out the problem */
      if(!(phi_cell > 0.0)) {
				/* DOUBLE tmp3 = (costhetamax - tmp1*tmp1)/(tmp2*tmp2); */
				/* fprintf(stderr,"ERROR: this_min_dec = %20.16lf phi_cell = %lf is negative. thetamax = %lf tmp1 = %lf tmp2 = %lf tmp3 = %lf \n",this_min_dec,phi_cell,thetamax,tmp1,tmp2,tmp3); */
				phi_cell = max_phi_cell;
      }
    }
    assert(phi_cell > 0.0 && "Making sure that RA bin-width is non-zero");
    phi_cell = phi_cell > max_phi_cell ? max_phi_cell:phi_cell;
		
    int nmesh_ra = (int) (phi_diff*phibin_refine_factor/phi_cell);
    if(nmesh_ra > NLATMAX)
      nmesh_ra = NLATMAX;

		if(nmesh_ra < (2*phibin_refine_factor + 1)) {
			nmesh_ra = 2*phibin_refine_factor + 1;
			fprintf(stderr,"%s> Using sub-optimal RA binning to ensure correct functioning of the code\n",__FUNCTION__);
		}
    /* fprintf(stderr,"idec = %d nmesh_ra = %d max_nmesh_phi = %d thetamax = %lf phi_diff = %lf phi_cell = %lf phi_cell/thetamax=%lf\n",idec,nmesh_ra,max_nmesh_phi,thetamax,phi_diff,phi_cell,phi_cell/thetamax); */
    assert(nmesh_ra <= max_nmesh_phi);

    ngrid_ra[idec] = nmesh_ra;
    
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
    assert(idec >=0 && idec < ngrid_dec);
    int ira  = (int)(ngrid_ra[idec]*(phi[i]-phi_min)*inv_phi_diff);
    if(ira >=ngrid_ra[idec]) ira--;
    assert(ira >=0 && ira < ngrid_ra[idec]);
    if(lattice[idec][ira].nelements == lattice[idec][ira].nallocated) {
      expected_n = lattice[idec][ira].nallocated*MEMORY_INCREASE_FAC;
			while(expected_n <= lattice[idec][ira].nelements || expected_n % NVEC != 0) {
				expected_n++;
			}

      const size_t memsize=3*sizeof(DOUBLE);
			lattice[idec][ira].pos = my_realloc(lattice[idec][ira].pos ,memsize,expected_n,"lattice.pos");
      /* lattice[idec][ira].cz = my_realloc(lattice[idec][ira].cz ,sizeof(*(lattice[idec][ira].cz)),expected_n,"lattice.cz"); */
      lattice[idec][ira].nallocated = expected_n;
    }
		assert(lattice[idec][ira].nallocated > lattice[idec][ira].nelements && "Making sure memory access if fine");
		/* index=lattice[idec][ira].nelements; */

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
    assigned_n++;
  }
  *max_in_cell = max_n;
#ifndef SILENT
  gettimeofday(&t1,NULL);
  fprintf(stderr,"%s> Allocated %0.2g (MB) memory for the lattice, expected_n = %d ngrid_dec = %d np=%"PRId64". Time taken = %6.2lf sec \n",__FUNCTION__, totnbytes/(1024*1024.),expected_n,ngrid_dec,np,
	  ADD_DIFF_TIME(t0,t1));
#endif
  /* fprintf(stderr,"np = %d assigned_n = %d\n",np,assigned_n); */
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


