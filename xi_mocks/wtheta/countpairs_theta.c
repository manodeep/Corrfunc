/* PROGRAM countpairs_theta

   --- countpairs_theta logthetamin logthetamax nbin N1 theta1 phi1 N2 theta2 phi2 Pairs

   --- Counts pairs of galaxies and bins them into an array of angular separation.
   ---inputs---
      * thetamin,thetamax,nbin = binning for Pairs array.
      * N1,theta1,phi1 = coords and dimension of first dataset 
      * N2,theta2,phi2 = coords and dimension of second dataset 
      * Pairs = array containing pairs (Pairs[theta])
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gridlink_mocks.h"//function proto-type for gridlink
#include "countpairs_theta.h" //function proto-type
#include "cellarray.h" //definition of struct cellarray_mocks
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar


#ifdef USE_AVX
#include "avx_calls.h"
#endif


void free_results_countpairs_theta(results_countpairs_theta **results)
{
  if(results == NULL)
	return;
  if(*results == NULL)
	return;

  results_countpairs_theta *tmp = *results;

  free(tmp->theta_upp);
  free(tmp->npairs);
  free(tmp->theta_avg);
  free(tmp);
  tmp = NULL;
}


void check_ra_dec(const int64_t N, DOUBLE *phi, DOUBLE *theta)
{
  int fix_ra  = 0;
  int fix_dec = 0;
  //Check input dec + ra
  for(int64_t i=0;i<N;i++) {
	if(phi[i] < 0.0) {
	  fix_ra = 1;
	  break;
	}
  }

  if(fix_ra == 1) {
	fprintf(stderr,"DDtheta> Out of range values found for ra. Expected ra to be in the range [0.0,360.0]. Found ra values in [-180,180] -- fixing that\n");
	for(int64_t i=0;i<N;i++) {
	  phi[i] += 180.0;
	}
  }

  for(int64_t i=0;i<N;i++) {
	if(theta[i] > 90.0) {
	  fix_dec = 1;
	}
  }

  if(fix_dec == 1) {
	fprintf(stderr,"DDtheta> Out of range values found for dec. Expected dec to be in the range [-90.0,90.0]. Found dec values in [0,180] -- fixing that\n");
	for(int64_t i=0;i<N;i++) {
	  theta[i] -= 90.0;
	}
  }
}


results_countpairs_theta * countpairs_theta(const int64_t ND1, DOUBLE *phi1, DOUBLE *theta1,
											const int64_t ND2, DOUBLE *phi2, DOUBLE *theta2,
#ifdef USE_OMP
											const int numthreads,
#endif
											const int autocorr,
											const char *binfile)
{

  fprintf(stderr,"Running in `%s' mode \n", autocorr == 1 ? "auto-correlation":"cross-correlation");
#if defined(LINK_IN_RA) && !defined(LINK_IN_DEC)
#error LINK_IN_DEC Makefile option must be enabled before LINK_IN_RA is selected
#endif 

#ifdef LINK_IN_DEC
  int rbin_refine_factor=1;
#ifdef LINK_IN_RA
  int phi_bin_refine_factor=1;
#endif
#endif

  //check the ra-dec inputs
  check_ra_dec(ND1, phi1, theta1);
  if(autocorr==0) {
	check_ra_dec(ND2, phi2, theta2);
  }


  double *theta_upp;
  int nthetabin;
  double thetamin,thetamax;
  setup_bins(binfile,&thetamin,&thetamax,&nthetabin,&theta_upp);
  assert(thetamin > 0.0 && thetamax > 0.0 && thetamin < thetamax && thetamax < 180.0 &&  "[thetamin, thetamax] are valid inputs");
  assert(nthetabin > 0 && "Number of theta bins is valid");
  uint64_t npairs[nthetabin];
  DOUBLE costheta_upp[nthetabin];
  for(int i=0;i<nthetabin;i++) {
	costheta_upp[i] = COSD(theta_upp[i]);
  }
  const DOUBLE costhetamin=costheta_upp[0];
  const DOUBLE costhetamax=costheta_upp[nthetabin-1];
#ifdef OUTPUT_THETAAVG  
  DOUBLE thetaavg[nthetabin];
  for(int i=0;i<nthetabin;i++) {
	thetaavg[i] = 0.0;
  }
#endif

#ifdef USE_AVX
  AVX_FLOATS m_costheta_upp[nthetabin] ;
  for(int i=0;i<nthetabin;i++) {  
    /* fprintf(stderr," i = %d theta_upp[i-1] = %lf cos(theta_upp[i-1] = %lf cos(theta_upp[i]) = %lf \n",i, theta_upp[i-1],COSD(theta_upp[i-1]),COSD(theta_upp[i])); */
    m_costheta_upp[i] = AVX_SET_FLOAT(costheta_upp[i]);
  }
/*   const AVX_FLOATS m_costhetamin=AVX_SET_FLOAT(costhetamin); */
  const AVX_FLOATS m_costhetamax=AVX_SET_FLOAT(costhetamax);
#ifdef OUTPUT_THETAVG
  AVX_FLOATS m_kbin[nthetabin];
  for(int i=0;i<nthetabin;i++) {
    m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
  }
#endif
#endif

#ifdef LINK_IN_DEC
  double dec_min=90.0,dec_max=-90.0;
#endif
  
#ifdef LINK_IN_RA
  double ra_min=360.0,ra_max=0.0;
#endif  

  for(int i=0;i<=nthetabin;i++) {
    npairs[i] = 0;
  }

  /*---Prepare-Data2--------------------------------*/
  DOUBLE *x2,*y2,*z2 ;
  x2=my_malloc(sizeof(*x2),ND2);
  y2=my_malloc(sizeof(*y2),ND2);
  z2=my_malloc(sizeof(*z2),ND2);
  
  for(int i=0;i<ND2;i++) {
    x2[i] = COSD(theta2[i])*COSD(phi2[i]) ;
    y2[i] = COSD(theta2[i])*SIND(phi2[i]) ;
    z2[i] = SIND(theta2[i]);

#ifdef LINK_IN_DEC
    if(theta2[i] < dec_min)
      dec_min = theta2[i];
    if(theta2[i] > dec_max)
      dec_max = theta2[i];
#endif

#ifdef LINK_IN_RA
    if(phi2[i] < ra_min)
      ra_min = phi2[i];
    if(phi2[i] > ra_max)
      ra_max = phi2[i];
#endif    
    
  }


#ifdef LINK_IN_DEC
  double dec_diff = dec_max-dec_min;
  double inv_dec_diff=1.0/dec_diff,decpos;
  int ngrid_dec=0,max_n=0;
#ifndef LINK_IN_RA

  cellarray_mocks *lattice2 = gridlink1D_theta(ND2, 
											  dec_min, dec_max, thetamax, 
											  x2, y2, z2,theta2,
											  &ngrid_dec,
											  &max_n,
											  rbin_refine_factor);

#else
  int *ngrid_ra=NULL;
  double inv_ra_diff=1.0/(ra_max-ra_min);
  fprintf(stderr,"ra_max = %lf ra_min = %lf\n",ra_max,ra_min);
 cellarray_mocks **lattice2 = gridlink2D_theta(ND2, dec_min, dec_max, thetamax,
										 x2, y2, z2,
										 theta2,
										 &ngrid_dec,
										 phi2,ra_min,ra_max,
										 &ngrid_ra,
										 &max_n,
										 rbin_refine_factor,
										 phi_bin_refine_factor);

#endif
  free(x2);free(y2);free(z2);
#endif

  int interrupted=0;
  int Nloop = 0;
  init_my_progressbar(ND1,&interrupted);

  /*---Loop-over-Data1-particles--------------------*/
  for(int i=0;i<ND1;i++) {
    my_progressbar(i,&interrupted);
    const DOUBLE x1 = COSD(theta1[i])*COSD(phi1[i]) ;
    const DOUBLE y1 = COSD(theta1[i])*SIND(phi1[i]) ;
    const DOUBLE z1 = SIND(theta1[i]) ;
    
#ifdef LINK_IN_DEC
    decpos = theta1[i];
    int dec_iz = (int)(ngrid_dec*(decpos-dec_min)*inv_dec_diff);
    if (dec_iz >= ngrid_dec) dec_iz--;
    assert(dec_iz >= 0 && dec_iz < ngrid_dec && "Declination is within lattice2 bounds");
    /* int dec_limits = (int) (ceil(thetamax*inv_dec_diff*ngrid_dec)); */
    int dec_limits = rbin_refine_factor;
    int min_dec = (dec_iz - dec_limits) <= 0 ? 0:dec_iz - dec_limits;
    int max_dec = (dec_iz + dec_limits) >= (ngrid_dec-1) ? (ngrid_dec-1):dec_iz + dec_limits;
    cellarray_mocks *cellstruct=NULL;    
    for(int idec=min_dec;idec<=max_dec;idec++) {
#ifndef LINK_IN_RA
      cellstruct = &(lattice2[idec]);
#else
      double rapos = phi1[i];
      int ra_iz = (int)(ngrid_ra[idec]*(rapos-ra_min)*inv_ra_diff);
      if (ra_iz >= ngrid_ra[idec]) ra_iz--;
      assert(ra_iz >= 0 && ra_iz < ngrid_ra[idec] && "RA position is within bounds");
      for(int ira_step=-phi_bin_refine_factor;ira_step<=phi_bin_refine_factor;ira_step++) {
		int ira = (ra_iz + ira_step + ngrid_ra[idec]) % ngrid_ra[idec];
		cellstruct = &(lattice2[idec][ira]);
#endif 
		x2 = cellstruct->x;
		y2 = cellstruct->y;
		z2 = cellstruct->z;
		Nloop = cellstruct->nelements;
#else //No linking in RA or DEC
		Nloop = ND2;
#endif
		
		
#ifdef USE_AVX
		const AVX_FLOATS m_x1 = AVX_SET_FLOAT(x1);
		const AVX_FLOATS m_y1 = AVX_SET_FLOAT(y1);
		const AVX_FLOATS m_z1 = AVX_SET_FLOAT(z1);
#endif    
		
	/*---Loop-over-Data2-particles--------------------*/
		int64_t j;
		for(j=0;j <=(Nloop-NVEC);j+=NVEC) {
#ifndef USE_AVX
		  DOUBLE costheta[NVEC];
		  int thetabin[NVEC];
#ifdef OUTPUT_THETAAVG
		  DOUBLE theta[NVEC];
#endif
		  int num_bad=0;
		  for(int k=0;k<NVEC;k++) {
			const DOUBLE tmp = x1*x2[j+k] + y1*y2[j+k] + z1*z2[j+k];
			costheta[k] = (tmp > 1.0) ? 1:(tmp < -1.0 ? -1.0:tmp);
#ifdef OUTPUT_THETAAVG
			theta[k]    =  INV_PI_OVER_180*ACOS(costheta[k]) ;
#endif
			if(costheta[k] > costhetamin || costheta[k] <= costhetamax) {
			  thetabin[k] = 0;
			  num_bad++;
			} else {
			  thetabin[k] = 1;//fill get filled in later
			}
		  }

		  //No pairs will be added just continue with the next iteration
		  if(num_bad == NVEC) {
			continue;
		  }
		  
		  //Now find the bins
		  for(int k=0;k<NVEC;k++) {
			if(thetabin[k]==0) continue;
			const DOUBLE this_cos_theta = costheta[k];
			for(int ibin=nthetabin-1;ibin>=1;ibin--) {
			  if(this_cos_theta <= costheta_upp[ibin-1]) {
				npairs[ibin]++;
				thetabin[k] = ibin;
				break;
			  }
			}
		  }

#ifdef OUTPUT_THETAAVG
#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif	  
	  for(int k=0;k<NVEC;k++) {
		thetaavg[thetabin[k]]+=theta[k];
	  }
#endif//OUTPUT_THETAAVG

#else

#ifdef OUTPUT_THETAAVG
	  union int8 {
		AVX_INTS m_ibin;
		int ibin[NVEC];
	  };
	  union int8 union_rpbin;
                                          
	  union float8{
		AVX_FLOATS m_Dperp;
		DOUBLE Dperp[NVEC];
	  };
	  union float8 union_mDperp;
#endif

	  //USE AVX intrinsics
	  const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(&x2[j]);
	  const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(&y2[j]);
	  const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(&z2[j]);
	  const AVX_FLOATS m_tmp1 = AVX_MULTIPLY_FLOATS(m_x2,m_x1);
	  const AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_y2,m_y1);
	  const AVX_FLOATS m_tmp3 = AVX_MULTIPLY_FLOATS(m_z2,m_z1);
	  const AVX_FLOATS m_costheta = AVX_ADD_FLOATS(m_tmp1,AVX_ADD_FLOATS(m_tmp2,m_tmp3));

	  AVX_FLOATS m_mask_left = AVX_COMPARE_FLOATS(m_costheta,m_costhetamax,_CMP_GT_OS);	  
/* 	  { */
/*      Does not seem to speed up the code */
/* 		AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_costheta,m_costhetamin,_CMP_LE_OS); */
/* 		AVX_FLOATS m_mask = AVX_BITWISE_AND(m1,m_mask_left); */
/* 		if(AVX_TEST_COMPARISON(m_mask) == 0) { */
/* 		  continue; */
/* 		} */
/* 	  } */

#ifdef OUTPUT_THETAAVG
	  //first do the acos to get the actual angles
	  const AVX_FLOATS m_inv_pi_over_180 = AVX_SET_FLOAT(INV_PI_OVER_180);
	  const AVX_FLOATS m_theta = AVX_ARC_COSINE(m_costheta);
	  union_mDperp.m_Dperp = AVX_MULTIPLY_FLOATS(m_theta,m_inv_pi_over_180);
	  AVX_FLOATS m_thetabin = AVX_SET_FLOAT((DOUBLE) 0.0);
#endif                                    


	  for(int kbin=nthetabin-1;kbin>=1;kbin--) {
	    const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_costheta,m_costheta_upp[kbin-1],_CMP_LE_OS);
	    const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
	    const int test = AVX_TEST_COMPARISON(m_bin_mask);
#ifdef OUTPUT_RPAVG
		m_thetabin = AVX_BLEND_FLOATS_WITH_MASK(m_thetabin,m_kbin[kbin], m_bin_mask);
#endif
	    npairs[kbin] += AVX_BIT_COUNT_INT(test);
		m_mask_left = AVX_COMPARE_FLOATS(m_costheta,m_costheta_upp[kbin-1],_CMP_GT_OS);
	    if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
	      break;
	    }
	  }

#ifdef OUTPUT_THETAAVG	  
	  union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_thetabin);
#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
	  for(int jj=0;jj<NVEC;jj++) {
		const int kbin = union_rpbin.ibin[jj];
		const DOUBLE theta = union_mDperp.Dperp[jj];
		thetaavg[kbin] += theta;
	  }
#endif


#endif//end of AVX section
	}//loop over particles in second data in chunks of NVEC
	
	
	//Take care of the remainder
	for(;j<Nloop;j++) {
	  const DOUBLE costheta = x1*x2[j] + y1*y2[j] + z1*z2[j] ;
	  if(costheta > costhetamin || costheta <= costhetamax) {
		continue;
	  }
	  
#ifdef OUTPUT_THETAAVG
	  const DOUBLE theta =  INV_PI_OVER_180*ACOS(costheta) ;
#endif	
	  for(int ibin=nthetabin-1;ibin>=1;ibin--) {
		if(costheta <= costheta_upp[ibin-1]) {
		  npairs[ibin]++;
#ifdef OUTPUT_THETAAVG
		  thetaavg[ibin] += theta;
#endif
		  break;
		}
	  }
	}//end of remainder loop


#ifdef LINK_IN_DEC
#ifdef LINK_IN_RA
      }
#endif    
      
    }//finish the loop over dec-cells in cellstruct
#endif    
  }//loop over ND1
  finish_myprogressbar(&interrupted);


#ifndef LINK_IN_DEC  
  free(x2);free(y2);free(z2);
#else
#ifndef LINK_IN_RA
  for(int i=0;i<ngrid_dec;i++) {
    free(lattice2[i].x);
    free(lattice2[i].y);
    free(lattice2[i].z);
  }
  free(lattice2);
#else
  for(int i=0;i<ngrid_dec;i++) {
    for(int j=0;j<ngrid_ra[i];j++) {
      free(lattice2[i][j].x);
      free(lattice2[i][j].y);
      free(lattice2[i][j].z);
    }
  }
  free(ngrid_ra);
  matrix_free((void **) lattice2,ngrid_dec);
#endif//LINK_IN_RA  
#endif//LINK_IN_DEC



  //prepare results
  //Pack in the results
  results_countpairs_theta *results = my_malloc(sizeof(*results), 1);
  results->nbin = nthetabin;
  results->npairs = my_malloc(sizeof(uint64_t), nthetabin);
  results->theta_upp   = my_malloc(sizeof(DOUBLE)  , nthetabin);
  results->theta_avg  = my_malloc(sizeof(DOUBLE)  , nthetabin);

  for(int i=0;i<nthetabin;i++) {
	results->npairs[i] = npairs[i];
	results->theta_upp[i] = theta_upp[i];
#ifdef OUTPUT_THETAAVG
	results->theta_avg[i] = thetaavg[i];
#else
	results->theta_avg[i] = 0.0;
#endif
  }
  free(theta_upp);

  return results;

}
