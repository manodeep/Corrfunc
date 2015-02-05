/* File: countpairs_data.c */
/*
                This file is a part of the Corrfunc package
                Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
                License: MIT LICENSE. See LICENSE file under the top-level
                directory at https://bitbucket.org/manodeep/corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "cellarray.h"
#include "gridlink_data.h"
#include "progressbar.h"
#include "countpairs_data.h"

#ifdef USE_AVX
#include "avx_calls.h"
#endif

#define SQR(X)         ((X) * (X))


void free_results_data(results_countpairs_data **results)
{
  if(results==NULL)
	return;

  if(*results==NULL)
	return;

  results_countpairs_data *tmp = *results;

  free(tmp->npairs);
  free(tmp->rupp);
  free(tmp->rpavg);
  free(tmp);
  tmp = NULL;
}



results_countpairs_data * countpairs_data(const int64_t N1, const DOUBLE *theta1, const DOUBLE *phi1, const DOUBLE *d1,
										  const int64_t N2, const DOUBLE *theta2, const DOUBLE *phi2, const DOUBLE *d2,
#ifdef USE_OMP  
										  const int numthreads,
#endif
										  const int autocorr,
										  const char *binfile,
										  const DOUBLE pimax)
{
  DOUBLE logrpmax,logrpmin,dlogrp,inv_dlogrp;
  DOUBLE dpi,inv_dpi;

  int zbin_refine_factor=2;
  
#ifdef LINK_IN_DEC
  int rbin_refine_factor=2;
#ifdef LINK_IN_RA
  int phibin_refine_factor=2;
#endif
#endif

#if defined(LINK_IN_RA) && !defined(LINK_IN_DEC)
  #error LINK_IN_DEC Makefile option must be enabled before LINK_IN_RA is selected
#endif
  const int npibin = (int) pimax;
  /***********************
   *initializing the  bins
   ************************/
  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins is valid");


  /*---Gridlink-variables----------------*/
  int ngrid;
  cellarray_data * restrict cellstruct;

  uint64_t npairs[(nrpbin+1)*(npibin+1)];
  DOUBLE rpavg[(nrpbin+1)*(npibin+1)];
  int index=0;
  for(int i=0; i <= nrpbin;i++) {
    for(int j=0;j <= npibin;j++) {
      /* fprintf(stderr,"rpbin = %d pibin = %d index = %d\n",i,j,index); */
      npairs[index] = 0;
      rpavg[index] = 0.0;
      index++;
    }
  }
  
  DOUBLE sqr_rpmin = rpmin*rpmin;
  DOUBLE sqr_rpmax = rpmax*rpmax;
  DOUBLE sqr_pimax = pimax*pimax;
  
  logrpmin = LOG2(rpmin) ;
  logrpmax = LOG2(rpmax) ;
  dlogrp = (logrpmax-logrpmin)/(DOUBLE)nrpbin ;
  inv_dlogrp = 1.0/dlogrp;

  dpi = pimax/(DOUBLE)npibin ;
  inv_dpi = 1.0/dpi;
  
#ifdef USE_AVX
  DOUBLE rupp_sqr[nrpbin+1];
  AVX_FLOATS m_rupp_sqr[nrpbin+1];
  rupp_sqr[0] = sqr_rpmin;
  for(int i=0;i<=nrpbin;i++) {
    rupp_sqr[i] = rupp[i]*rupp[i];
    m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
  }
  /* AVX_FLOATS m_piupp_sqr[npibin+1]; */
  /* m_piupp_sqr[0] = AVX_SET_FLOAT((DOUBLE) 0.0); */
  /* for(int i=1;i<=npibin;i++) { */
  /*   DOUBLE this_pi = i*dpi; */
  /*   m_piupp_sqr[i] = AVX_SET_FLOAT(this_pi*this_pi); */
  /* } */
#endif
    
  
  
  /*---Prepare-Data2--------------------------------*/

#ifdef LINK_IN_DEC
  DOUBLE dec_min=90.0,dec_max=-90.0;
#endif

#ifdef LINK_IN_RA
  DOUBLE ra_min=360.0,ra_max=0.0;
#endif
  

  DOUBLE d2min=1000. ;
  DOUBLE d2max=0. ;

  get_max_min_data(N1, d1, &d2min, &d2max,
#ifdef LINK_IN_DEC
				   theta1,
				   &dec_min,&dec_max,
#endif

#ifdef LINK_IN_RA
				   phi1,
				   &ra_min,&ra_max
#endif
				   );
  
  if(autocorr==0) {
	get_max_min_data(N2, d2, &d2min, &d2max,
#ifdef LINK_IN_DEC
					 theta2,
					 &dec_min,&dec_max,
#endif
					 
#ifdef LINK_IN_RA
					 phi2,
					 &ra_min,&ra_max
#endif
					 );
  }


  ngrid=0;
  int max_n;
  cellarray_data *lattice1=NULL,*lattice2=NULL;
  int64_t totncells;

#ifdef LINK_IN_DEC
  int *ngrid_dec;
  const DOUBLE dec_diff = dec_max - dec_min;
  const DOUBLE inv_dec_diff=1.0/dec_diff;
  /* fprintf(stderr,"dec_min = %lf dec_max = %lf\n",dec_min,dec_max); */
#ifndef LINK_IN_RA  
  lattice2 = gridlink2D(N2,d2min,d2max,pimax,
						dec_min,dec_max,rpmax,
						theta2,phi2,
						d2, theta2,
						&ngrid, &ngrid_dec, &max_n);
#else
  //Linking in cz, Dec, RA
  DOUBLE inv_ra_diff=1.0/(ra_max-ra_min);
  int **ngrid_ra=NULL;
  /* fprintf(stderr,"ra_min = %lf ra_max = %lf\n",ra_min,ra_max); */
  lattice2 = gridlink3D(N2,d2min,d2max,pimax,
					   dec_min,dec_max,rpmax,
					   theta2,phi2,
					   d2,
					   theta2,
					   &ngrid,
					   &ngrid_dec,
					   phi2, ra_min,ra_max,
					   &ngrid_ra,
					   &max_n);

#endif
  //Need cz_binsize for LINK_IN_DEC option
  const DOUBLE cz_binsize=(d2max-d2min)/ngrid;
#else
  //Only linking in cz
  lattice2 = gridlink1D(N2, d2min, d2max, pimax, theta2, phi2, d2, &ngrid, &max_n);
  fprintf(stderr,"countpairs> gridlink1D done. ngrid= %d max_n = %d\n",ngrid,max_n);
#endif


  int interrupted=0;
  init_my_progressbar(N1,&interrupted);
  
  /*---Loop-over-Data1-particles--------------------*/
  for(int i=0;i<N1;i++) {
    my_progressbar(i,&interrupted);
	
    const DOUBLE x1 = d1[i]*COS(theta1[i]*PI_OVER_180)*COS(phi1[i]*PI_OVER_180) ;
	const DOUBLE y1 = d1[i]*COS(theta1[i]*PI_OVER_180)*SIN(phi1[i]*PI_OVER_180) ;
    const DOUBLE z1 = d1[i]*SIN(theta1[i]*PI_OVER_180) ;

    /*---Deterpmine-central-grid-cell-of-search--------*/
    int icen = (int)(ngrid*(d1[i]-d2min)/(d2max-d2min)) ;
    if(icen<0) icen = 0 ;
    if(icen>=ngrid) icen = ngrid-1 ;
    
    int min_iz,max_iz;
    min_iz = (icen - zbin_refine_factor) <= 0 ? 0:icen - zbin_refine_factor;
    max_iz = (icen + zbin_refine_factor) >= (ngrid-1) ? (ngrid-1):icen + zbin_refine_factor;

    /*---Loop-over-surrounding-cells----------------*/
    for(int icell=min_iz;icell<=max_iz;icell++) {
#ifdef LINK_IN_DEC

      DOUBLE decpos = theta1[i];
      /* double dmin_iz = (icell < icen) ? icell:icen; */
      /* dmin_iz *= (d2max-d2min)/ngrid; */
      /* dmin_iz += d2min; */
      DOUBLE dmin_iz = d2min + ((icen + icell)*cz_binsize)*0.5;
      /* double theta = rpmax/(2.0*dmin_iz); */
      DOUBLE max_dec_sep=asin(rpmax/(2*dmin_iz))*2.0*INV_PI_OVER_180;
      int dec_limits = (int) (ceil(max_dec_sep*inv_dec_diff*ngrid_dec[icell]));
      /* fprintf(stderr,"icen = %d ngrid_dec[%d] = %d rbin_refine_factor=%d dec_limits=%d\n", */
      /* 	      icen,icell,ngrid_dec[icell],RBIN_REFINE_FACTOR,dec_limits); */
      /* int dec_limits = RBIN_REFINE_FACTOR; */
      int dec_iz = (int)(ngrid_dec[icell]*(decpos-dec_min)*inv_dec_diff);
      if(dec_iz>=ngrid_dec[icell]) dec_iz-- ;
      if(!( dec_iz >=0 && dec_iz < ngrid_dec[icell])) {
				fprintf(stderr,"icell = %d ngrid_dec[icell] = %d dec_iz = %d decpos = %lf\n",icell,ngrid_dec[icell],dec_iz,decpos);
      }
      assert(dec_iz >=0 && dec_iz < ngrid_dec[icell] && "Declination inside bounds");
      const int min_dec = (dec_iz - dec_limits) <= 0 ? 0:dec_iz - dec_limits;
      const int max_dec = (dec_iz + dec_limits) >= (ngrid_dec[icell]-1) ? (ngrid_dec[icell]-1):dec_iz + dec_limits;

	  
      for(int idec=min_dec;idec<=max_dec;idec++) {
#ifdef LINK_IN_RA	
		DOUBLE rapos = phi1[i];
		int ra_iz = (int)(ngrid_ra[icell][idec]*(rapos-ra_min)*inv_ra_diff);
		if (ra_iz >= ngrid_ra[icell][idec]) ra_iz--;
		assert(ra_iz >= 0  && ra_iz < ngrid_ra[icell][idec] && "RA position within bounds");
		int ra_limits = PHI_BIN_REFINE_FACTOR;
		for(int ira_step=-ra_limits;ira_step<=ra_limits;ira_step++) {
		  int ira = (ra_iz + ira_step + ngrid_ra[icell][idec]) % ngrid_ra[icell][idec];
		  //Linked in CZ, DEC and RA
		  cellstruct = &(lattice[icell][idec][ira]);
#else
		  //Linked in CZ + DEC
		  cellstruct = &(lattice[icell][idec]);
#endif
		  
#else
		  //LINKED only in CZ
		  cellstruct = &(lattice[icell]);
#endif      
		  const DOUBLE *x2  = cellstruct->x;
		  const DOUBLE *y2  = cellstruct->y;
		  const DOUBLE *z2  = cellstruct->z;
		  const DOUBLE *cz2 = cellstruct->cz;
		  /* gettimeofday(&t0,NULL); */
		  
		  /* DOUBLE Dpar,Dperp; */
		  /* DOUBLE parx,pary,parz,perpx,perpy,perpz,tmp; */
		  /* DOUBLE HALF=0.5; */
		  DOUBLE TWO=2.0;
		  /* int rpbin,pibin; */
		  DOUBLE sqr_d1 = d1[i]*d1[i];
		  /* int last_rpbin = (nrpbin-1)*(npibin+1); */
		  
		  
#ifndef USE_AVX
		  int64_t j;
		  for(j=0;j<=(cellstruct->nelements-NVEC);j+=NVEC) {
			/* 					  int vec_rpbins[NVEC]; */
			int vec_pibins[NVEC];
			DOUBLE Dperp[NVEC];
#pragma simd vectorlengthfor(DOUBLE)
			for(int jj=0;jj<NVEC;jj++) {
			  DOUBLE sqr_cz = cz2[j+jj]*cz2[j+jj];
			  DOUBLE tmp = (sqr_d1 - sqr_cz);
			  DOUBLE xy_costheta = x1*x2[j+jj] + y1*y2[j+jj] + z1*z2[j+jj];
			  DOUBLE tmp1 = (sqr_d1 + sqr_cz + TWO*xy_costheta);
			  DOUBLE Dpar = SQR(tmp)/tmp1;
			  vec_pibins[jj] = (Dpar >= sqr_pimax) ? npibin:(int) (SQRT(Dpar)*inv_dpi);
			  DOUBLE tmp2 = sqr_d1 + sqr_cz -TWO*xy_costheta - Dpar;
			  Dperp[jj]  = (Dpar >= sqr_pimax || tmp2 >= sqr_rpmax || tmp2 < sqr_rpmin) ? 0.0: tmp2;
			  /* 							vec_rpbins[jj] = (Dperp[jj] >= sqr_rpmax || Dperp[jj] < sqr_rpmin) ? nrpbin:(int)((LOG2(Dperp[jj])*0.5-logrpmin)*inv_dlogrp); */
			}
			
#pragma unroll(NVEC)
			for(int jj=0;jj<NVEC;jj++) {
			  
			  npairs[vec_rpbins[jj]*(npibin+1) + vec_pibins[jj]]++;
			  rpavg[vec_rpbins[jj]*(npibin+1) + vec_pibins[jj]]+= SQRT(Dperp[jj]);
			}
		  }
#else //Use AVX intrinsics
		  AVX_FLOATS m_xpos    = AVX_BROADCAST_FLOAT(&x1);
		  AVX_FLOATS m_ypos    = AVX_BROADCAST_FLOAT(&y1);
		  AVX_FLOATS m_zpos    = AVX_BROADCAST_FLOAT(&z1);
		  AVX_FLOATS m_sqr_d1  = AVX_BROADCAST_FLOAT(&sqr_d1);
		  union int8 {
			AVX_INTS m_ibin;
			int ibin[NVEC];
		  };
		  union int8 union_rpbin;
		  union int8 union_pibin;
		  
		  union float8{
			AVX_FLOATS m_Dperp;
			DOUBLE Dperp[NVEC];
		  };
		  union float8 union_mDperp;
					
		  /* AVX_FLOATS m_half = AVX_SET_FLOAT(HALF); */
		  /* AVX_FLOATS m_quarter = AVX_SET_FLOAT((DOUBLE) 0.25); */
		  AVX_FLOATS m_sqr_pimax  = AVX_SET_FLOAT(sqr_pimax);
		  AVX_FLOATS m_sqr_rpmax  = AVX_SET_FLOAT(sqr_rpmax);
		  AVX_FLOATS m_sqr_rpmin  = AVX_SET_FLOAT(sqr_rpmin);
		  /* AVX_FLOATS m_logrpmin   = AVX_SET_FLOAT(logrpmin); */
		  /* AVX_FLOATS m_inv_dlogrp = AVX_SET_FLOAT(inv_dlogrp); */
		  AVX_FLOATS m_npibin     = AVX_SET_FLOAT((DOUBLE) npibin);
		  /* AVX_FLOATS m_nrpbin     = AVX_SET_FLOAT((DOUBLE) nrpbin); */
		  AVX_FLOATS m_zero       = AVX_SET_FLOAT((DOUBLE) 0.0);
		  
		  for(j=0;j<=(cellstruct->nelements-NVEC);j+=NVEC){
			AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(&x2[j]);
			AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(&y2[j]);
			AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(&z2[j]);
			AVX_FLOATS m_cz2 = AVX_LOAD_FLOATS_UNALIGNED(&cz2[j]);
			AVX_FLOATS m_sqr_cz2 = AVX_SQUARE_FLOAT(m_cz2); 
			AVX_FLOATS m_sum_of_norms = AVX_ADD_FLOATS(m_sqr_d1,m_sqr_cz2);
			AVX_FLOATS m_inv_dpi    = AVX_SET_FLOAT(inv_dpi);	    
			
			AVX_FLOATS m_twice_xy_costheta;
			{
			  AVX_FLOATS m_tmp1 = AVX_MULTIPLY_FLOATS(m_xpos,m_x2);
			  AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_ypos,m_y2);
			  AVX_FLOATS m_tmp3 = AVX_ADD_FLOATS(m_tmp1,m_tmp2);
			  AVX_FLOATS m_tmp4 = AVX_MULTIPLY_FLOATS(m_zpos,m_z2);
			  m_twice_xy_costheta = AVX_ADD_FLOATS(m_tmp3,m_tmp4);
			  m_twice_xy_costheta = AVX_ADD_FLOATS(m_twice_xy_costheta,m_twice_xy_costheta);
			}
			
			AVX_FLOATS m_Dpar;
			{
			  /* AVX_FLOATS m_tmp  = AVX_MULTIPLY_FLOATS(m_half,AVX_SUBTRACT_FLOATS(m_sqr_d1,m_sqr_cz2)); */
			  AVX_FLOATS m_tmp  = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_sqr_d1,m_sqr_cz2));
			  AVX_FLOATS m_tmp1 = AVX_ADD_FLOATS(m_sum_of_norms,m_twice_xy_costheta);
			  /* AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_quarter,m_tmp1); */
			  m_Dpar = AVX_DIVIDE_FLOATS(m_tmp,m_tmp1);
			}
			
			/* AVX_FLOATS m_Dpar  = AVX_MULTIPLY_FLOATS(m_tmp,AVX_RECIPROCAL_FLOATS(m_tmp1)); */
			AVX_FLOATS m_Dperp;
			{
			  AVX_FLOATS m_tmp1 = AVX_SUBTRACT_FLOATS(m_sum_of_norms,m_twice_xy_costheta);
			  m_Dperp = AVX_SUBTRACT_FLOATS(m_tmp1,m_Dpar);
			}
			
			AVX_FLOATS m_mask;
			{
			  {
				AVX_FLOATS m_tmp1 = AVX_COMPARE_FLOATS(m_Dpar,m_sqr_pimax,_CMP_LT_OS);
				AVX_FLOATS m_tmp2 = AVX_COMPARE_FLOATS(m_Dperp,m_sqr_rpmax,_CMP_LT_OS);
				AVX_FLOATS m_tmp3 = AVX_COMPARE_FLOATS(m_Dperp,m_sqr_rpmin,_CMP_GE_OS);
				AVX_FLOATS m_tmp4 = AVX_BITWISE_AND(m_tmp1,m_tmp2);
				m_mask = AVX_BITWISE_AND(m_tmp3,m_tmp4);
				int test = AVX_TEST_COMPARISON(m_mask);
				if(test==0)
				  continue;
			  }
			  m_Dperp = AVX_BLEND_FLOATS_WITH_MASK(m_zero,m_Dperp,m_mask);
			  m_Dpar  = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_pimax,m_Dpar,m_mask);
			  union_mDperp.m_Dperp = AVX_BLEND_FLOATS_WITH_MASK(m_zero,AVX_SQRT_FLOAT(m_Dperp),m_mask);
			  
			  {
				AVX_FLOATS m_mask_left = AVX_COMPARE_FLOATS(m_Dperp,m_sqr_rpmax,_CMP_LT_OS);
				AVX_FLOATS m_rpbin = AVX_SET_FLOAT((DOUBLE) nrpbin);
				for(int kbin=nrpbin-1;kbin>=0;kbin--) {
				  AVX_FLOATS m_mask_low = AVX_COMPARE_FLOATS(m_Dperp,m_rupp_sqr[kbin],_CMP_GE_OS);
				  AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m_mask_low,m_mask_left);
				  AVX_FLOATS m_bin = AVX_SET_FLOAT((DOUBLE) kbin);
				  m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_bin, m_bin_mask);
				  m_mask_left = AVX_COMPARE_FLOATS(m_Dperp, m_rupp_sqr[kbin],_CMP_LT_OS);
				  int test = AVX_TEST_COMPARISON(m_mask_left);
				  if(test==0)
					break;
				}
				union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
			  }
			  
			  {
				AVX_FLOATS m_tmp1 = AVX_SQRT_FLOAT(m_Dpar);
				AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_tmp1,m_inv_dpi);
				AVX_FLOATS m_pibin = AVX_BLEND_FLOATS_WITH_MASK(m_npibin, m_tmp2, m_mask);
				union_pibin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_pibin);
			  }
			}
			
#pragma unroll(NVEC)
			for(int jj=0;jj<NVEC;jj++) {
			  npairs[union_rpbin.ibin[jj]*(npibin+1) + union_pibin.ibin[jj]]++;
			  rpavg [union_rpbin.ibin[jj]*(npibin+1) + union_pibin.ibin[jj]] += union_mDperp.Dperp[jj];
			  /* fprintf(stderr,"i=%d j=%d union_rpbin.ibin[jj] = %d union_pibin.ibin[jj] = %d\n",i,j,union_rpbin.ibin[jj],union_pibin.ibin[jj]); */
			}
		  }
#endif	//END of the AVX/NO-AVX section
		  
		  //Take care of the rest
		  for(;j<cellstruct->nelements;j++) {
			int rpbin,pibin;
			DOUBLE sqr_cz = cz2[j]*cz2[j];
			DOUBLE tmp = (sqr_d1 - sqr_cz);
			DOUBLE xy_costheta = x1*x2[j] + y1*y2[j] + z1*z2[j];
			DOUBLE tmp1 = (sqr_d1 + sqr_cz + TWO*xy_costheta);
			DOUBLE Dpar = SQR(tmp)/tmp1;
			pibin  = (Dpar >= sqr_pimax) ? npibin:(int) (SQRT(Dpar)*inv_dpi);
			DOUBLE tmp2  = sqr_d1 + sqr_cz -TWO*xy_costheta - Dpar;
			DOUBLE Dperp = (Dpar >= sqr_pimax || tmp2 >= sqr_rpmax || tmp2 < sqr_rpmin) ? 0.0:tmp2;
			rpbin  = (Dperp == 0.0) ? nrpbin:(int)((LOG2(Dperp)*0.5-logrpmin)*inv_dlogrp);
			npairs[rpbin*(npibin+1) + pibin]++;
			rpavg [rpbin*(npibin+1) + pibin]+=SQRT(Dperp);
		  }
		}
#ifdef LINK_IN_DEC
#ifdef LINK_IN_RA
      }
#endif 	
    }
#endif      
  }
  finish_myprogressbar(&interrupted);
  /* fprintf(stderr,"simd_time = %6.2lf serial_time = %6.2lf sec\n",simd_time,serial_time); */
  index=0;
  for(int i=0;i<=nrpbin;i++) {
    for(int j=0;j<=npibin;j++) {
      if(npairs[index] > 0) {
		rpavg[index] /= (DOUBLE) npairs[index] ;
      }
      index++;
    }
  }
	
  for(int64_t i=0;i < totncells;i++) {
    free(lattice1[i].x);
    free(lattice1[i].y);
    free(lattice1[i].z);
    free(lattice1[i].cz);
	if(autocorr==0) {
	  free(lattice2[i].x);
	  free(lattice2[i].y);
	  free(lattice2[i].z);
	  free(lattice2[i].cz);
	}
  }
  free(lattice1);
  if(autocorr==0) {
	free(lattice2);
  }


  //rp's are all in log2 -> convert to log10
/*   const double inv_log10=1.0/log2(10); */
/*   for(int i=0;i<nrpbin;i++) { */
/*     DOUBLE logrp = logrpmin + (DOUBLE)(i+1)*dlogrp; */
/*     for(int j=0;j<npibin;j++) { */
/*       index = i*(npibin+1) + j; */
/*       fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",npairs[index],rpavg[index],logrp*inv_log10,(j+1)*dpi); */
/*     } */
/*   } */


  //Pack in the results
  results_countpairs_data *results = my_malloc(sizeof(*results), 1);
  results->nbin   = nrpbin;
  results->npibin = npibin;
  results->pimax  = pimax;
  results->npairs = my_malloc(sizeof(uint64_t), totnbins);
  results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
  results->rpavg  = my_malloc(sizeof(DOUBLE)  , totnbins);

  for(int i=0;i<nrpbin;i++) {
	results->rupp[i] = rupp[i];
	for(int j=0;j<npibin;j++) {
	  index = i*(npibin+1) + j;
	  results->npairs[index] = npairs[index];
#ifdef OUTPUT_RPAVG
	  results->rpavg[index] = rpavg[index];
#else
	  results->rpavg[index] = 0.0;
#endif
	}
  }

  free(rupp);
}



