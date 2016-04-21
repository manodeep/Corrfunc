/* File: countpairs_rp_pi_kernels.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

#include "defs.h"
#include "utils.h"
#include "function_precision.h"


#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"

static inline void countpairs_rp_pi_avx_intrinsics(DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                                                   DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell, 
                                                   const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const int npibin, const DOUBLE *rupp_sqr, const DOUBLE pimax
#ifdef PERIODIC                                             
                                                   ,const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#endif                                             
                                                   ,const AVX_FLOATS *m_rupp_sqr

                                                   ,const AVX_FLOATS *m_kbin
#ifdef OUTPUT_RPAVG                                                   
                                                   ,DOUBLE *src_rpavg
#endif                         
                                                   ,uint64_t *src_npairs)
{
  const int64_t totnbins = (npibin+1)*(nbin+1);
  uint64_t *npairs = my_calloc(sizeof(*npairs), totnbins);
  const DOUBLE dpi = pimax/npibin;
  const DOUBLE inv_dpi = 1.0/dpi;

#ifdef OUTPUT_RPAVG
  DOUBLE rpavg[totnbins];
  for(int i=0;i<totnbins;i++) {
    rpavg[i] = ZERO;
  }
#endif
    
  for(int64_t i=0;i<N0;i++) {

#ifdef PERIODIC        
    const DOUBLE xpos = *x0++ + off_xwrap;
    const DOUBLE ypos = *y0++ + off_ywrap;
    const DOUBLE zpos = *z0++ + off_zwrap;
#else
    const DOUBLE xpos = *x0++;
    const DOUBLE ypos = *y0++;
    const DOUBLE zpos = *z0++;
#endif        
    DOUBLE *localz1 = ((DOUBLE *) z1);

    int64_t j= 0;
    if(same_cell == 1) {
      j = i+1;
      localz1 += j;
    } else {
      while(j < N1){
        const DOUBLE dz = *localz1++ - zpos;
        if(dz > -pimax) break;
        j++;
      }
      localz1--;
    }
    DOUBLE *localx1 = ((DOUBLE *) x1) + j;
    DOUBLE *localy1 = ((DOUBLE *) y1) + j;

    for(;j<=(N1 - AVX_NVEC);j+=AVX_NVEC) {
      const AVX_FLOATS m_xpos    = AVX_SET_FLOAT(xpos);
      const AVX_FLOATS m_ypos    = AVX_SET_FLOAT(ypos);
      const AVX_FLOATS m_zpos    = AVX_SET_FLOAT(zpos);
            
      union int8 {
        AVX_INTS m_ibin;
        int ibin[NVEC];
      };
      union int8 union_finalbin;

#ifdef OUTPUT_RPAVG
      union float8{
        AVX_FLOATS m_Dperp;
        DOUBLE Dperp[NVEC];
      };
      union float8 union_mDperp;
#endif
            
      const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(localx1);
      const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(localy1);
      const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(localz1);
            
      localx1 += AVX_NVEC;//this might actually exceed the allocated range but we will never dereference that
      localy1 += AVX_NVEC;
      localz1 += AVX_NVEC;

      const AVX_FLOATS m_pimax = AVX_SET_FLOAT((DOUBLE) pimax);
      const AVX_FLOATS m_sqr_rpmax = m_rupp_sqr[nbin-1];
      const AVX_FLOATS m_sqr_rpmin = m_rupp_sqr[0];
      const AVX_FLOATS m_inv_dpi    = AVX_SET_FLOAT(inv_dpi);
      const AVX_FLOATS m_zero = AVX_SET_FLOAT(ZERO);
      const AVX_FLOATS m_npibin     = AVX_SET_FLOAT((DOUBLE) npibin);
      const AVX_FLOATS m_one    = AVX_SET_FLOAT((DOUBLE) 1);

      AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_zpos, m_z1);
      const AVX_FLOATS m_sqr_xdiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_xpos,m_x1));//(x0 - x[j])^2
      const AVX_FLOATS m_sqr_ydiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_ypos,m_y1));//(y0 - y[j])^2
      /* const AVX_FLOATS m_sqr_zdiff = AVX_SQUARE_FLOAT(m_zdiff); */
      AVX_FLOATS r2  = AVX_ADD_FLOATS(m_sqr_xdiff, m_sqr_ydiff);
      m_zdiff = AVX_MAX_FLOATS(m_zdiff,AVX_SUBTRACT_FLOATS(m_zero,m_zdiff));//dz = fabs(dz) => dz = max(dz, -dz);
            
      AVX_FLOATS m_mask_left;
            
      //Do all the distance cuts using masks here in new scope
      {
        //the z2 arrays are sorted in increasing order. which means
        //the z2 value will increase in any future iteration of j.
        //that implies the zdiff values are also monotonically increasing
        //Therefore, if none of the zdiff values are less than pimax, then
        //no future iteration in j can produce a zdiff value less than pimax.
        AVX_FLOATS m_mask_pimax = AVX_COMPARE_FLOATS(m_zdiff,m_pimax,_CMP_LT_OS);
        if(AVX_TEST_COMPARISON(m_mask_pimax) == 0) {
          j=N1;
          break;
        }
                
        const AVX_FLOATS m_rpmax_mask = AVX_COMPARE_FLOATS(r2, m_sqr_rpmax, _CMP_LT_OS);
        const AVX_FLOATS m_rpmin_mask = AVX_COMPARE_FLOATS(r2, m_sqr_rpmin, _CMP_GE_OS);
        const AVX_FLOATS m_rp_mask = AVX_BITWISE_AND(m_rpmax_mask,m_rpmin_mask);
                
        //Create a combined mask by bitwise and of m1 and m_mask_left.
        //This gives us the mask for all sqr_rpmin <= r2 < sqr_rpmax
        m_mask_left = AVX_BITWISE_AND(m_mask_pimax,m_rp_mask);
                
        //If not, continue with the next iteration of j-loop
        if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
          continue;
        }
                
        //There is some r2 that satisfies sqr_rpmin <= r2 < sqr_rpmax && 0.0 <= dz^2 < pimax^2.
        r2 = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2, m_mask_left);
        m_zdiff = AVX_BLEND_FLOATS_WITH_MASK(m_pimax, m_zdiff, m_mask_left);
      }
            
#ifdef OUTPUT_RPAVG
      union_mDperp.m_Dperp = AVX_SQRT_FLOAT(r2);
#endif
            
      const AVX_FLOATS m_pibin = AVX_MULTIPLY_FLOATS(m_zdiff,m_inv_dpi);
      AVX_FLOATS m_rpbin     = AVX_SET_FLOAT((DOUBLE) 0);
      //AVX_FLOATS m_all_ones  = AVX_CAST_INT_TO_FLOAT(AVX_SET_INT(-1));
      for(int kbin=nbin-1;kbin>=1;kbin--) {
        const AVX_FLOATS m_mask_low = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_GE_OS);
        const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m_mask_low,m_mask_left);
        m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
        m_mask_left = AVX_COMPARE_FLOATS(r2, m_rupp_sqr[kbin-1],_CMP_LT_OS);
        //m_mask_left = AVX_XOR_FLOATS(m_mask_low, m_all_ones);//XOR with 0xFFFF... gives the bins that are smaller than m_rupp_sqr[kbin] (and is faster than cmp_p(s/d) in theory)
        const int test = AVX_TEST_COMPARISON(m_mask_left);
        if(test==0) {
          break;
        }
      }
      const AVX_FLOATS m_npibin_p1 = AVX_ADD_FLOATS(m_npibin,m_one);
      const AVX_FLOATS m_binproduct = AVX_ADD_FLOATS(AVX_MULTIPLY_FLOATS(m_rpbin,m_npibin_p1),m_pibin);
      union_finalbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_binproduct);

      //update the histograms
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma unroll(NVEC)
#endif
      for(int jj=0;jj<NVEC;jj++) {
        int ibin = union_finalbin.ibin[jj];
        npairs[ibin]++;
#ifdef OUTPUT_RPAVG
        rpavg [ibin] += union_mDperp.Dperp[jj];
#endif
      }
    }

            
    //remainder loop 
    for(;j<N1;j++){
      const DOUBLE dz = FABS(*localz1++ - zpos);
      const DOUBLE dx = *localx1++ - xpos;
      const DOUBLE dy = *localy1++ - ypos;

      if(dz >= pimax) {
        break;
      }
            
      const DOUBLE r2 = dx*dx + dy*dy;
      if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
        continue;
      }
            
#ifdef OUTPUT_RPAVG
      const DOUBLE r = SQRT(r2);
#endif

      int pibin = (int) (dz*inv_dpi);
      pibin = pibin > npibin ? npibin:pibin;
      for(int kbin=nbin-1;kbin>=1;kbin--) {
        if(r2 >= rupp_sqr[kbin-1]) {
          int ibin = kbin*(npibin+1) + pibin;
          npairs[ibin]++;
#ifdef OUTPUT_RPAVG
          rpavg [ibin] += r;
#endif
          break;
        }
      }

    }//remainder loop over second set of particles
  }//loop over first set of particles

	for(int i=0;i<totnbins;i++) {
		src_npairs[i] += npairs[i];
#ifdef OUTPUT_RPAVG
    src_rpavg[i] += rpavg[i];
#endif        
  }
  free(npairs);
    
}

#endif //__AVX__



#if defined (__SSE4_2__)
#include "sse_calls.h"

static inline void countpairs_rp_pi_sse_intrinsics(DOUBLE *x0, DOUBLE *y0, DOUBLE *z0, const int64_t N0,
                                                   DOUBLE *x1, DOUBLE *y1, DOUBLE *z1, const int64_t N1, const int same_cell,
                                                   const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const int npibin, const DOUBLE rupp_sqr[] , const DOUBLE pimax
#ifdef PERIODIC                                             
                                                   ,const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#endif                                             
                                                   ,const SSE_FLOATS *m_rupp_sqr
                                                   ,const SSE_FLOATS *m_kbin
#ifdef OUTPUT_RPAVG
                                                   ,DOUBLE *src_rpavg
#endif                         
                                                   ,uint64_t *src_npairs)
{

  const int64_t totnbins = (npibin+1) * (nbin+1);
  uint64_t *npairs = my_calloc(sizeof(*npairs), totnbins);
  const DOUBLE dpi = pimax/npibin;
  const DOUBLE inv_dpi = 1.0/dpi;

#ifdef OUTPUT_RPAVG
  DOUBLE rpavg[totnbins];
  for(int i=0;i<totnbins;i++) {
    rpavg[i] = ZERO;
  }
#endif

  for(int64_t i=0;i<N0;i++) {

#ifdef PERIODIC        
    const DOUBLE xpos = *x0++ + off_xwrap;
    const DOUBLE ypos = *y0++ + off_ywrap;
    const DOUBLE zpos = *z0++ + off_zwrap;
#else
    const DOUBLE xpos = *x0++;
    const DOUBLE ypos = *y0++;
    const DOUBLE zpos = *z0++;
#endif        

    DOUBLE *localz1 = ((DOUBLE *) z1) ;
    int64_t j= 0;
    if(same_cell == 1) {
      j = i+1;
      localz1 += j;
    } else {
      while(j < N1){
        const DOUBLE dz = *localz1++ - zpos;
        if(dz > -pimax) break;
        j++;
      }
      localz1--;
    }
    DOUBLE *localx1 = ((DOUBLE *) x1) + j;
    DOUBLE *localy1 = ((DOUBLE *) y1) + j;
        
		for(;j<=(N1 - SSE_NVEC);j+=SSE_NVEC){

      union int4{
        SSE_INTS m_ibin;
        int ibin[SSE_NVEC];
      };
      union int4 union_finalbin;

#ifdef OUTPUT_RPAVG            
      union float4{
        SSE_FLOATS m_Dperp;
        DOUBLE Dperp[SSE_NVEC];
      };
      union float4 union_mDperp;
#endif
      const SSE_FLOATS m_xpos = SSE_SET_FLOAT(xpos);
      const SSE_FLOATS m_ypos = SSE_SET_FLOAT(ypos);
      const SSE_FLOATS m_zpos = SSE_SET_FLOAT(zpos);

      const SSE_FLOATS m_x1 = SSE_LOAD_FLOATS_UNALIGNED(localx1);
      const SSE_FLOATS m_y1 = SSE_LOAD_FLOATS_UNALIGNED(localy1);
      const SSE_FLOATS m_z1 = SSE_LOAD_FLOATS_UNALIGNED(localz1);
            
      localx1 += SSE_NVEC;
      localy1 += SSE_NVEC;
      localz1 += SSE_NVEC;
            
      const SSE_FLOATS m_pimax = SSE_SET_FLOAT((DOUBLE) pimax);
      const SSE_FLOATS m_sqr_rpmax = m_rupp_sqr[nbin-1];
      const SSE_FLOATS m_sqr_rpmin = m_rupp_sqr[0];
      const SSE_FLOATS m_inv_dpi    = SSE_SET_FLOAT(inv_dpi);
      const SSE_FLOATS m_zero = SSE_SET_FLOAT(ZERO);
      const SSE_FLOATS m_npibin     = SSE_SET_FLOAT((DOUBLE) npibin);
      const SSE_FLOATS m_one    = SSE_SET_FLOAT((DOUBLE) 1);

      SSE_FLOATS m_zdiff = SSE_SUBTRACT_FLOATS(m_zpos, m_z1);
      const SSE_FLOATS m_sqr_xdiff = SSE_SQUARE_FLOAT(SSE_SUBTRACT_FLOATS(m_xpos,m_x1));//(x0 - x[j])^2
      const SSE_FLOATS m_sqr_ydiff = SSE_SQUARE_FLOAT(SSE_SUBTRACT_FLOATS(m_ypos,m_y1));//(y0 - y[j])^2
      SSE_FLOATS r2  = SSE_ADD_FLOATS(m_sqr_xdiff, m_sqr_ydiff);
      m_zdiff = SSE_MAX_FLOATS(m_zdiff,SSE_SUBTRACT_FLOATS(m_zero,m_zdiff));//dz = fabs(dz) => dz = max(dz, -dz);
            
      SSE_FLOATS m_mask_left;
            
      //Do all the distance cuts using masks here in new scope
      {
        //the z2 arrays are sorted in increasing order. which means
        //the z2 value will increase in any future iteration of j.
        //that implies the zdiff values are also monotonically increasing
        //Therefore, if none of the zdiff values are less than pimax, then
        //no future iteration in j can produce a zdiff value less than pimax.
        SSE_FLOATS m_mask_pimax = SSE_COMPARE_FLOATS_LT(m_zdiff,m_pimax);
        if(SSE_TEST_COMPARISON(m_mask_pimax) == 0) {
          j=N1;
          break;
        }
                
        const SSE_FLOATS m_rpmax_mask = SSE_COMPARE_FLOATS_LT(r2, m_sqr_rpmax);
        const SSE_FLOATS m_rpmin_mask = SSE_COMPARE_FLOATS_GE(r2, m_sqr_rpmin);
        const SSE_FLOATS m_rp_mask = SSE_BITWISE_AND(m_rpmax_mask,m_rpmin_mask);
                
        //Create a combined mask by bitwise and of m1 and m_mask_left.
        //This gives us the mask for all sqr_rpmin <= r2 < sqr_rpmax
        m_mask_left = SSE_BITWISE_AND(m_mask_pimax,m_rp_mask);
                
        //If not, continue with the next iteration of j-loop
        if(SSE_TEST_COMPARISON(m_mask_left) == 0) {
          continue;
        }
                
        //There is some r2 that satisfies sqr_rpmin <= r2 < sqr_rpmax && 0.0 <= dz^2 < pimax^2.
        r2 = SSE_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2, m_mask_left);
        m_zdiff = SSE_BLEND_FLOATS_WITH_MASK(m_pimax, m_zdiff, m_mask_left);
      }
            
#ifdef OUTPUT_RPAVG
      union_mDperp.m_Dperp = SSE_SQRT_FLOAT(r2);
#endif
            
      const SSE_FLOATS m_pibin = SSE_MULTIPLY_FLOATS(m_zdiff,m_inv_dpi);
      SSE_FLOATS m_rpbin     = SSE_SET_FLOAT((DOUBLE) 0);
      //SSE_FLOATS m_all_ones  = SSE_CAST_INT_TO_FLOAT(SSE_SET_INT(-1));
      for(int kbin=nbin-1;kbin>=1;kbin--) {
        const SSE_FLOATS m_mask_low = SSE_COMPARE_FLOATS_GE(r2,m_rupp_sqr[kbin-1]);
        const SSE_FLOATS m_bin_mask = SSE_BITWISE_AND(m_mask_low,m_mask_left);
        m_rpbin = SSE_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
        m_mask_left = SSE_COMPARE_FLOATS_LT(r2, m_rupp_sqr[kbin-1]);
        //m_mask_left = SSE_XOR_FLOATS(m_mask_low, m_all_ones);//XOR with 0xFFFF... gives the bins that are smaller than m_rupp_sqr[kbin] (and is faster than cmp_p(s/d) in theory)
        const int test = SSE_TEST_COMPARISON(m_mask_left);
        if(test==0) {
          break;
        }
      }
      const SSE_FLOATS m_npibin_p1 = SSE_ADD_FLOATS(m_npibin,m_one);
      const SSE_FLOATS m_binproduct = SSE_ADD_FLOATS(SSE_MULTIPLY_FLOATS(m_rpbin,m_npibin_p1),m_pibin);
      union_finalbin.m_ibin = SSE_TRUNCATE_FLOAT_TO_INT(m_binproduct);

      //update the histograms
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma unroll(NVEC)
#endif
      for(int jj=0;jj<SSE_NVEC;jj++) {
        int ibin = union_finalbin.ibin[jj];
        npairs[ibin]++;
#ifdef OUTPUT_RPAVG
        rpavg [ibin] += union_mDperp.Dperp[jj];
#endif
      }
    }


		for(;j<N1;j++) {
			const DOUBLE dx = *localx1++ - xpos;
			const DOUBLE dy = *localy1++ - ypos;
			const DOUBLE dz = FABS(*localz1++ - zpos);
      if(dz >= pimax) break;
            
			const DOUBLE r2 = dx*dx + dy*dy;
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
#ifdef OUTPUT_RPAVG
      const DOUBLE r = SQRT(r2);
#endif            

      int pibin = (int) (dz*inv_dpi);
      pibin = pibin > npibin ? npibin:pibin;
      for(int kbin=nbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
          int ibin = kbin*(npibin+1) + pibin;
					npairs[ibin]++;
#ifdef OUTPUT_RPAVG
          rpavg[ibin] += r;
#endif                    
					break;
				}
			}//searching for kbin loop
		}
  }
	for(int i=0;i<totnbins;i++) {
		src_npairs[i] += npairs[i];
#ifdef OUTPUT_RPAVG
    src_rpavg[i] += rpavg[i];
#endif        
        
	}
	free(npairs);
}

#endif //__SSE4_2__
