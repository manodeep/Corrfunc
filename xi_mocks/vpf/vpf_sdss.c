/* PROGRAM vpf_sdss

   --- vpf_sdss rmax nbin nc volume SDSSfile RANDfile > output
   --- compute the void probability function for SDSS galaxies

      * rmax = maximum radius (in h^-1 Mpc)
      * nbin = number of radii (evenly spaced in r)
      * nc = number of centers to place (does not count rejected centers)
      * volume = volume of sample (in Mpc^3/h^3)
      * galaxy file, (ascii format - contains: ra,dec,cz)
      * random file, (ascii format - contains: ra,dec,cz)
      > output: <R P0 P1 P2 Navg Nvar>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_interp.h>

#include "proto.h"
#include "utils.h"
#include "cellarray.h"
#include "set_cosmo_dist.h"
#include "progressbar.h"
#include "ftread.h"


#ifdef USE_AVX
#include "avx_calls.h"
#endif

void Printhelp(void) ;

int main(int argc, char *argv[])
{
  int i;
  /*---Arguments-------------------------*/
  int nbin,nc ;
  DOUBLE volume,rmax ;

  /*---Particle-distribution-variables---*/
  int64_t Ngal,Nran;
  DOUBLE theta,phi,cz,dc;
  DOUBLE *xgal,*ygal,*zgal,*xran,*yran,*zran;
  /*---Gridlink-variables----------------*/
  int ngrid;

  //*---Measurement-----------------------*/
  int itry,isucceed;
  /* int maxbin,p,Nnbrs_gal,Nnbrs_ran,*nbrs_gal,*nbrs_ran; */
  int Nnbrs_ran,*counts ;
  double xcen,ycen,zcen,rcube,Nexpected ;
  double *R,*P0,*P1,*P2,*Navg,*Nvar ;
  /* double cz2comovingdist(double) ; */
  struct timeval t0,t1,tstart;
  double loop_time=0.0,geometry_time=0.0,conversion_time=0.0;
	int lasdamas_cosmology=1;
  const char argnames[][100]={"rmax","nbin","ncenters","volume","galaxy file","galaxy file-format","randoms file","randoms file-format","centers file","cosmology flag ->1,2 (1->LD)"};
  int nargs=sizeof(argnames)/(sizeof(char)*100);

  gettimeofday(&tstart,NULL);
  /*---Read-arguments-----------------------------------*/
  if(argc < (nargs+1)) {
    Printhelp() ;
    fprintf(stderr,"\nFound: %d parameters\n ",argc-1);
    for(i=1;i<argc;i++) {
      if(i <= nargs)
				fprintf(stderr,"\t\t %25s = `%s' \n",argnames[i-1],argv[i]);
      else
				fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
    }
    if(i <= nargs) {
      fprintf(stderr,"\nMissing required parameters \n");
      for(i=argc;i<=nargs;i++)
				fprintf(stderr,"\t\t %s = `?'\n",argnames[i-1]);
    }
    return EXIT_FAILURE;
  }
  sscanf(argv[1],"%lf",&rmax);
  sscanf(argv[2],"%d",&nbin);
  sscanf(argv[3],"%d",&nc);
  sscanf(argv[4],"%"DOUBLE_FORMAT,&volume);

  fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
  fprintf(stderr,"\n\t\t -------------------------------------\n");
  for(i=1;i<argc;i++) {
    if(i <= nargs) {
      fprintf(stderr,"\t\t %-25s = %s \n",argnames[i-1],argv[i]);
    }  else {
      fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
    }
  }
  fprintf(stderr,"\t\t -------------------------------------\n");
  assert(rmax > 0 && "rmax must be > 0");
  assert(nbin > 0 && "Number of bins must be > 0");
  assert(nc > 0 && "Number of spheres must be > 0");
  assert(volume > 0 && "Mock volume must be > 0");
  lasdamas_cosmology = atoi(argv[10]);
  assert(lasdamas_cosmology >=1 && lasdamas_cosmology <= 2 && "Cosmology flag has to be 1 (LasDamas) or 2 (Planck)");
  
  const char *galaxy_file=argv[5];
  const char *galaxy_file_format=argv[6];
  const char *random_file=argv[7];
  const char *random_file_format=argv[8];
  FILE *fpcen = fopen(argv[9],"r");
  int generate_randoms = 1;
  int64_t num_centers_in_file=0;
  if(fpcen != NULL) {
    double rr = 0.0;
    fscanf(fpcen,"%*lf %*lf %*lf %lf",&rr);
    num_centers_in_file = getnumlines(argv[9],'#');
    if( rr >= rmax && num_centers_in_file >= nc) {
      generate_randoms = 0;
      rewind(fpcen);
    } else {
      fclose(fpcen);
      num_centers_in_file = 0;
      fpcen = my_fopen(argv[9],"w");
      generate_randoms = 1;
    }
  } else {
    num_centers_in_file = 0;
    fpcen = my_fopen(argv[9],"w");
    generate_randoms = 1;
  }
  fprintf(stderr,"vpf_sdss> found %"PRId64" centers (need %d centers) - generating randoms = %d\n",num_centers_in_file,nc,generate_randoms);
  
  //set up the interpolation for comoving distances
  double *zz,*ddc;
  zz=my_calloc(sizeof(*zz),COSMO_DIST_SIZE);
  ddc=my_calloc(sizeof(*ddc),COSMO_DIST_SIZE);
  int Nzdc = set_cosmo_dist(MAX_REDSHIFT_FOR_COSMO_DIST, COSMO_DIST_SIZE, zz, ddc, lasdamas_cosmology);

  gsl_interp *interpolation;
  gsl_interp_accel *accelerator;
  accelerator =  gsl_interp_accel_alloc();
  interpolation = gsl_interp_alloc (gsl_interp_linear,Nzdc);
  gsl_interp_init(interpolation, zz, ddc, Nzdc);
  
  rcube=0 ;
  gettimeofday(&t0,NULL);
  /*---Read-galaxy-data1-file----------------------------------*/
  Ngal=read_positions(galaxy_file,galaxy_file_format, sizeof(DOUBLE), 3, &xgal, &ygal, &zgal);
  for(i=0;i<Ngal;i++) {
	DOUBLE new_phi   = xgal[i];
	DOUBLE new_theta = ygal[i];
	DOUBLE new_cz    = zgal[i];
	dc = gsl_interp_eval(interpolation, zz, ddc, new_cz/SPEED_OF_LIGHT, accelerator);
	if(dc>rcube) rcube = dc;
	xgal[i] = dc*COSD(new_theta)*COSD(new_phi) ;
	ygal[i] = dc*COSD(new_theta)*SIND(new_phi) ;
	zgal[i] = dc*SIND(new_theta) ;
  }

  gettimeofday(&t1,NULL);
  fprintf(stderr,"vpf_sdss> Ngal = %d. Time to read-in galaxies and do the projections =%6.2lf sec\n",Ngal,ADD_DIFF_TIME(t1,t0)) ;

  /*---Read-random-file---------------------------------*/
  if(generate_randoms == 1) {
    gettimeofday(&t0,NULL);
	Nran = read_positions(random_file,random_file_fileformat, sizeof(DOUBLE), 3, &xran, &yran, &zran);
    gettimeofday(&t1,NULL);
    fprintf(stderr,"vpf_sdss> Nran = %d. Time to read-in randoms = %6.2lf sec\n",Nran,ADD_DIFF_TIME(t1,t0)) ;
  }    

  if (generate_randoms == 1) {
    gettimeofday(&t0,NULL);
    for(i=0;i<Nran;i++) {
      DOUBLE new_phi   = xran[i];
      DOUBLE new_theta = yran[i];
      DOUBLE new_cz    = zran[i];
      dc = gsl_interp_eval(interpolation, zz, ddc, new_cz/SPEED_OF_LIGHT, accelerator);
      
      xran[i] = dc*COSD(new_theta)*COSD(new_phi) ;
      yran[i] = dc*COSD(new_theta)*SIND(new_phi) ;
      zran[i] = dc*SIND(new_theta);
    }
    gettimeofday(&t1,NULL);
    conversion_time += ADD_DIFF_TIME(t1,t0);
  }
  free(zz);free(ddc);
  gsl_interp_free(interpolation);
  gsl_interp_accel_free(accelerator);
  
  
  /*---Expected-number-of-randoms-in-sphere-------------*/
  if(generate_randoms == 1) {
    Nexpected = (double)Nran*(4.0*M_PI*(rmax*rmax*rmax)/3.)/volume ;
    fprintf(stderr,"vpf_sdss> Expected number of randoms in sphere = %f\n",Nexpected) ;
  }

  /*---Shift-coordinates--------------------------------*/
  fprintf(stderr,"vpf_sdss> maximum distance = %f\n",rcube) ;
  rcube = rcube + 1. ; //add buffer

  gettimeofday(&t0,NULL);
  for(i=0;i<Ngal;i++) {
    xgal[i] += rcube ;
    ygal[i] += rcube ;
    zgal[i] += rcube ;
  }

  if(generate_randoms == 1) {
    for(i=0;i<Nran;i++) {
      xran[i] += rcube ;
      yran[i] += rcube ;
      zran[i] += rcube ;
    }
  }
  gettimeofday(&t1,NULL);
  conversion_time += ADD_DIFF_TIME(t1,t0);

  rcube = 2*rcube ;
  DOUBLE inv_rstep = ((DOUBLE) nbin)/rmax;
  DOUBLE inv_rcube = 1.0/rcube;
  DOUBLE rmax_sqr = rmax*rmax;
  fprintf(stderr,"vpf_sdss> cube size = %f\n",rcube) ;

  /*---Construct-grid-to-speed-up-neighbor-searching----*/
  //First create the 3-d linklist
  cellarray *lattice=NULL;//pointer to the full 3-d volume for galaxies
  cellarray *randoms_lattice=NULL;//pointer to the full 3-d volume for randoms
  cellarray *cellstruct=NULL;//to be used as a pointer to one cell
  const DOUBLE smin=0.0;
  DOUBLE smax=rcube;
  gettimeofday(&t0,NULL);
  ngrid=0 ;
  lattice = linklist_4d(Ngal,xgal,ygal,zgal,smin,smax,rmax,&ngrid);
  if(generate_randoms == 1) {
    randoms_lattice = linklist_4d(Nran,xran,yran,zran,smin,smax,rmax,&ngrid);
  }
  gettimeofday(&t1,NULL);
  fprintf(stderr,"vpf_sdss> gridlink done. ngrid = %d. Time taken to create both galaxy and random lattice = %6.2lf sec\n",ngrid,ADD_DIFF_TIME(t1,t0)) ;
  

  /*---Prepare-radial-arrays----------------------------*/
  counts = my_calloc(sizeof(*counts),nbin);
  R      = my_calloc(sizeof(*R),nbin);
  Navg   = my_calloc(sizeof(*Navg),nbin);
  Nvar   = my_calloc(sizeof(*Nvar),nbin);
  P0     = my_calloc(sizeof(*P0),nbin);
  P1     = my_calloc(sizeof(*P1),nbin);
  P2     = my_calloc(sizeof(*P1),nbin);

#ifdef USE_AVX
  AVX_FLOATS m_rupp_sqr[nbin];
  AVX_FLOATS m_rmax_sqr = AVX_SET_FLOAT(rmax_sqr);
#endif  
  
  for(int k=0;k<nbin;k++) {
    R[k] = (double)(k+1)*rmax/(double)nbin ;
#ifdef USE_AVX
		m_rupp_sqr[k] = AVX_SET_FLOAT(R[k]*R[k]);
#endif    
  }

  /*---Loop-over-random-centers-------------------------*/
  itry=0 ;
  isucceed=0 ;
  int interrupted;
  int threshold_neighbors;
  if(generate_randoms == 1) {
		threshold_neighbors = Nexpected - 3*sqrt(Nexpected);
  } else {
    threshold_neighbors = 1;//dummy value -> just to ensure that the check does not compare with uninitialized values
  }
  int ncenters_written=0;
  init_my_progressbar(nc, &interrupted);
  while(isucceed < nc && itry < Nran) {
    my_progressbar(isucceed,&interrupted);

    if((generate_randoms == 1 && isucceed > num_centers_in_file) || num_centers_in_file == 0) {
      const DOUBLE xcen = xran[itry] ;
      const DOUBLE ycen = yran[itry] ;
      const DOUBLE zcen = zran[itry] ;
      gettimeofday(&t0,NULL);
      Nnbrs_ran = count_neighbors(xcen,ycen,zcen,smin,inv_rcube,rmax,ngrid,randoms_lattice, threshold_neighbors, bin_refine_factor);
      gettimeofday(&t1,NULL);
      geometry_time += ADD_DIFF_TIME(t1,t0);
    } else {
      double rr=0.0;
      const int MAXBUFSIZE=10000;
      char buffer[MAXBUFSIZE];
      assert( fgets(buffer,MAXBUFSIZE,fpcen) != NULL); 
      int nitems = sscanf(buffer,"%lf %lf %lf %lf",&xcen,&ycen,&zcen,&rr);
      if(nitems != 4) {
				fprintf(stderr,"ERROR: nitems = %d xcen = %lf ycen = %lf zcen %lf rr = %lf\n",
								nitems,xcen,ycen,zcen,rr);
				fprintf(stderr,"buffer = `%s' \n",buffer);
      }
      assert(nitems == 4 && "Read the centers from the centers file");
      assert(rr >= rmax && "Rmax from the center file is >= rmax");
      Nnbrs_ran = threshold_neighbors + 1;
    }
    
    if(Nnbrs_ran > threshold_neighbors) {  //ignore if sphere overlaps edge
      for(int k=0;k<nbin;k++) {  //initialize counts
				counts[k] = 0 ;
      }
      
      int ix = (int)(ngrid*(xcen-smin)*inv_rcube);
      int iy = (int)(ngrid*(ycen-smin)*inv_rcube);
      int iz = (int)(ngrid*(zcen-smin)*inv_rcube);
      if(ix > ngrid-1) ix--;
      if(iy > ngrid-1) iy--;
      if(iz > ngrid-1) iz--;

      assert(ix >= 0 && ix < ngrid && "x-position is inside limits");
      assert(iy >= 0 && iy < ngrid && "y-position is inside limits");
      assert(iz >= 0 && iz < ngrid && "z-position is inside limits");
      
      gettimeofday(&t0,NULL);
      int min_ix = ix - BIN_REFINE_FACTOR < 0 ?             0:ix - BIN_REFINE_FACTOR;
      int max_ix = ix + BIN_REFINE_FACTOR > ngrid-1 ? ngrid-1:ix + BIN_REFINE_FACTOR;
      for(int iix=min_ix;iix<=max_ix;iix++) {
				DOUBLE newxpos = xcen;
#ifdef USE_AVX
				AVX_FLOATS m_newxpos = AVX_SET_FLOAT(newxpos);
#endif	
	
				int min_iy = iy - BIN_REFINE_FACTOR < 0 ?             0:iy - BIN_REFINE_FACTOR;
				int max_iy = iy + BIN_REFINE_FACTOR > ngrid-1 ? ngrid-1:iy + BIN_REFINE_FACTOR;

				for(int iiy=min_iy;iiy<=max_iy;iiy++) {
					DOUBLE newypos = ycen;
#ifdef USE_AVX
					AVX_FLOATS m_newypos = AVX_SET_FLOAT(newypos);
#endif	

					int min_iz = iz - BIN_REFINE_FACTOR < 0 ?             0:iz - BIN_REFINE_FACTOR;
					int max_iz = iz + BIN_REFINE_FACTOR > ngrid-1 ? ngrid-1:iz + BIN_REFINE_FACTOR;

					for(int iiz=min_iz;iiz<=max_iz;iiz++) {
						DOUBLE newzpos = zcen;
#ifdef USE_AVX
						AVX_FLOATS m_newzpos = AVX_SET_FLOAT(newzpos);
#endif	
						int index=iix*ngrid*ngrid + iiy*ngrid + iiz;
						cellstruct = &(lattice[index]);
						DOUBLE *x2 = cellstruct->x;
						DOUBLE *y2 = cellstruct->y;
						DOUBLE *z2 = cellstruct->z;
						int ipart;
						for(ipart=0;ipart<=(cellstruct->nelements-NVEC);ipart+=NVEC) {
#ifndef USE_AVX
							int ibin[NVEC];
#pragma simd vectorlengthfor(DOUBLE)
							for(int k=0;k<NVEC;k++) {
								DOUBLE dx,dy,dz,r;
								dx=x2[ipart+k]-newxpos;
								dy=y2[ipart+k]-newypos;
								dz=z2[ipart+k]-newzpos;
								r = SQRT(dx*dx + dy*dy + dz*dz);
								ibin[k] = (int) (r*inv_rstep);
							}
#pragma unroll(NVEC)
							for(int k=0;k<NVEC;k++) {
								if(ibin[k] < nbin) counts[ibin[k]]++;
							}

							//Here is the AVX part
#else
							AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(&x2[ipart]);
							AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(&y2[ipart]);
							AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(&z2[ipart]);
							AVX_FLOATS m_xdiff = AVX_SUBTRACT_FLOATS(m_x2,m_newxpos);
							AVX_FLOATS m_ydiff = AVX_SUBTRACT_FLOATS(m_y2,m_newypos);
							AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_z2,m_newzpos);
							AVX_FLOATS m_dist  = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_xdiff),AVX_SQUARE_FLOAT(m_ydiff));
							m_dist = AVX_ADD_FLOATS(m_dist,AVX_SQUARE_FLOAT(m_zdiff));
							AVX_FLOATS m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_rmax_sqr,_CMP_LT_OS);
							int test = AVX_TEST_COMPARISON(m_mask_left);
							if(test == 0)
								continue;
				
							for(int kbin=nbin-1;kbin>=1;kbin--) {
								AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_GE_OS);
								AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
								int test2  = AVX_TEST_COMPARISON(m_bin_mask);
								counts[kbin] += AVX_BIT_COUNT_INT(test2);
								m_mask_left = AVX_COMPARE_FLOATS(m_dist,m_rupp_sqr[kbin-1],_CMP_LT_OS);
								int test3 = AVX_TEST_COMPARISON(m_mask_left);
								if(test3 == 0) {
									break;
								} else if(kbin==1){
									counts[0] += AVX_BIT_COUNT_INT(test3);
								}
							}

#endif //endof AVX section	      
						}

						//Take care of the rest
						for(;ipart < cellstruct->nelements;ipart++) {
							DOUBLE dx,dy,dz,r2;
							int ibin;
							dx=x2[ipart]-newxpos;
							dy=y2[ipart]-newypos;
							dz=z2[ipart]-newzpos;
							r2 = (dx*dx + dy*dy + dz*dz);
							if(r2 >= rmax_sqr) continue;
							ibin = (int) (SQRT(r2)*inv_rstep);
							counts[ibin]++;
						}
					}
				}
      }
      gettimeofday(&t1,NULL);
      loop_time += ADD_DIFF_TIME(t1,t0);

      //Output the center into the file -> either 
      if((generate_randoms == 1 && isucceed > num_centers_in_file) || num_centers_in_file == 0) {
				fprintf(fpcen,"%lf \t %lf \t %lf \t %lf\n",xcen,ycen,zcen,rmax);
				ncenters_written++;
      }

      
      gettimeofday(&t0,NULL);
      /* compute cumulative counts, i.e. n1 changes from the number of galaxies
				 in shell ibin to  the number of galaxies in shell ibin or any smaller shell */
      for(int ibin=1;ibin<nbin;ibin++){
				counts[ibin]+=counts[ibin-1];
      }
      
      for(int k=0;k<nbin;k++) { //compute statistics
      	Navg[k] += (double)counts[k] ;
      	Nvar[k] += (double)(counts[k]*counts[k]) ;
      	if(counts[k]==0) P0[k]++ ;
      	if(counts[k]==1) P1[k]++ ;
      	if(counts[k]==2) P2[k]++ ;
      }
      isucceed++ ;
    }
    itry++ ;
  }
  fclose(fpcen);
  finish_myprogressbar(&interrupted);
  fprintf(stderr,"vpf_sdss> Placed %d centers out of %d trials. loop_time = %6.2lf sec geometry_time = %6.2lf sec conversion_time = %6.2lf sec\n",isucceed,itry,loop_time,geometry_time,conversion_time);
  fprintf(stderr,"vpf_sdss> num_centers_in_file = %"PRId64" ncenters_written = %d\n",num_centers_in_file,ncenters_written);
  assert(isucceed > 0 && "Placed > 0 spheres in the volume");
  if(isucceed < nc) {
    fprintf(stderr,"WARNING: Could only place `%d' out of requested `%d' spheres. Increase the random-sample size might improve the situation\n",isucceed,nc);
  }
  
  double inv_nspheres = 1.0/isucceed;
  for(int k=0;k<nbin;k++) { //compute statistics
    Navg[k] *= inv_nspheres;
    Nvar[k] = Nvar[k]*inv_nspheres - (Navg[k]*Navg[k]) ;
    P0[k] *= inv_nspheres;
    P1[k] *= inv_nspheres;
    P2[k] *= inv_nspheres;
    
    fprintf(stdout,"%4.1f  %8.3e %8.3e %8.3e  %8.2f %8.2f\n",R[k],P0[k],P1[k],P2[k],Navg[k],Nvar[k]) ;
  }	  

  free(xgal);free(ygal);free(zgal);
  if(generate_randoms == 1) {
    free(xran);free(yran);free(zran);
  }
  free(Navg);free(Nvar);free(P0);free(P1);free(P2);
  free(counts);free(R);
  int64_t totncells = ngrid*ngrid*ngrid;
  for(int64_t icell=0;icell < totncells;icell++) {
    free(lattice[icell].x);
    if(generate_randoms == 1) {
      free(randoms_lattice[icell].x);
    }
  }

  free(lattice);
  if(generate_randoms == 1) {
    free(randoms_lattice);
  }
  gettimeofday(&t1,NULL);
  fprintf(stderr,"vpf_sdss> Done. Ngal = %d. Time taken = %6.2lf sec\n",Ngal,ADD_DIFF_TIME(t1,t0));
  
  return EXIT_SUCCESS ;
}

/*---Print-help-information---------------------------*/

void Printhelp(void)
{
  fprintf(stderr,"=========================================================================\n") ;
  fprintf(stderr,"   --- vpf_sdss rmax nbin nc volume SDSSfile RANDfile > output\n") ;
  fprintf(stderr,"   --- compute the void probability function for SDSS galaxies\n") ;
  fprintf(stderr,"      * rmax = maximum radius (in h^-1 Mpc)\n") ;
  fprintf(stderr,"      * nbin = number of radii (evenly spaced in r)\n") ;
  fprintf(stderr,"      * nc = number of centers to place (does not count rejected centers)\n") ;
  fprintf(stderr,"      * volume = volume of sample (in Mpc^3/h^3)\n") ;
  fprintf(stderr,"      * galaxy file, (contains: ra,dec,cz)\n") ;
  fprintf(stderr,"      * galaxy file format (a -> ascii, f->fast-food)\n") ;
  fprintf(stderr,"      * random file, (contains: ra,dec,cz)\n") ;
  fprintf(stderr,"      * random file format (a -> ascii, f-> fast-food)\n");
  fprintf(stderr,"      * file with sphere centers (centers will be read-in if enough centers exist, otherwise centers will be output into this file)\n");
  fprintf(stderr,"      > output: <R P0 P1 P2 Navg Nvar>\n") ;
  fprintf(stderr,"=========================================================================\n") ;
}



