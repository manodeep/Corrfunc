/* File: countpairs_wp.c */
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
#include "function_precision.h"
#include "cellarray.h" //definition of struct cellarray_index
#include "utils.h" //all of the utilities
#include "gridlink.h"//function proto-type for gridlink
#include "countpairs_wp.h" //function proto-type
#include "progressbar.h" //for the progressbar

#include "sglib.h"

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif

#ifdef USE_OMP
#include <omp.h>
#endif


#ifndef PERIODIC
#warning "wp is only valid for PERIODIC boundary conditions. Ignoring the Makefile (non)-definition of PERIODIC"
#endif

void free_results_wp(results_countpairs_wp **results)
{
    if(results == NULL)
        return;

    if(*results == NULL)
        return;

    results_countpairs_wp *tmp = *results;
    free(tmp->npairs);
    free(tmp->rupp);
    free(tmp->wp);
    free(tmp->rpavg);
    free(tmp);
    tmp = NULL;
}

static inline void same_cell_wp_kernel(const cellarray_index * first,
                                       DOUBLE * restrict const X,DOUBLE * restrict const Y, DOUBLE * restrict const Z, 
                                       const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE rupp_sqr[] , const DOUBLE pimax
#ifdef USE_AVX
                                            ,const AVX_FLOATS m_rupp_sqr[] 
#ifdef OUTPUT_RPAVG
                                            ,const AVX_FLOATS m_kbin[]
#endif                                            
#endif                                            


#ifdef OUTPUT_RPAVG
                                       ,DOUBLE src_rpavg[]
#endif                         
                                       ,uint64_t *src_npairs)
{
    uint64_t npair[nbin];
    for(int i=0;i<nbin;i++) {
        npair[i] = 0;
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nbin];
    for(int i=0;i<nbin;i++) {
        rpavg[i] = 0.0;
    }
#endif    

    DOUBLE *x1 = X + first->start;
    DOUBLE *y1 = Y + first->start;
    DOUBLE *z1 = Z + first->start;
    
    for(int64_t i=0;i<first->nelements;i++) {
        const DOUBLE x1pos = *x1++;
        const DOUBLE y1pos = *y1++;
        const DOUBLE z1pos = *z1++;
        
        DOUBLE *x2 = x1;
        DOUBLE *y2 = y1;
        DOUBLE *z2 = z1;

        int64_t j=i+1;
        for(;j<=(first->nelements-NVEC);j+=NVEC) {
#if !(defined(USE_AVX) && defined(__AVX__))
            for(int jj=0;jj<NVEC;jj++) {
                const DOUBLE dx = *x2++ - x1pos;
                const DOUBLE dy = *y2++ - y1pos;
                const DOUBLE dz = *z2++ - z1pos;
                
                if(dz >= pimax) {
                    j = first->nelements;
                    break;
                }

                
                const DOUBLE r2 = dx*dx + dy*dy;
                if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
                    continue;
                }
                
#ifdef OUTPUT_RPAVG
                const DOUBLE r = SQRT(r2);
#endif
                for(int kbin=nbin-1;kbin>=1;kbin--) {
                    if(r2 >= rupp_sqr[kbin-1]) {
                        npair[kbin]++;
#ifdef OUTPUT_RPAVG
                        rpavg[kbin]+=r;
#endif
                        break;
                    }
                }
            }//end of jj loop
#else
            const AVX_FLOATS m_xpos    = AVX_SET_FLOAT(x1pos);
            const AVX_FLOATS m_ypos    = AVX_SET_FLOAT(y1pos);
            const AVX_FLOATS m_zpos    = AVX_SET_FLOAT(z1pos);
            
#ifdef OUTPUT_RPAVG
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
            const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(x2);
            const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(y2);
            const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(z2);
            
            x2 += NVEC;//this might actually exceed the allocated range but we will never dereference that
            y2 += NVEC;
            z2 += NVEC;
            
            const AVX_FLOATS m_pimax = AVX_SET_FLOAT(pimax);
            const AVX_FLOATS m_zero  = AVX_SET_FLOAT(ZERO);
            
            const AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_z2,m_zpos);//z2[j:j+NVEC-1] - z1
            const AVX_FLOATS m_sqr_xdiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_xpos,m_x2));//(x0 - x[j])^2
            const AVX_FLOATS m_sqr_ydiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_ypos,m_y2));//(y0 - y[j])^2
            AVX_FLOATS r2  = AVX_ADD_FLOATS(m_sqr_xdiff,m_sqr_ydiff);
            
            AVX_FLOATS m_mask_left;
            
            //Do all the distance cuts using masks here in new scope
            {
                //the z2 arrays are sorted in increasing order. which means
                //the z2 value will increase in any future iteration of j.
                //that implies the zdiff values are also monotonically increasing
                //Therefore, if none of the zdiff values are less than pimax, then
                //no future iteration in j can produce a zdiff value less than pimax.
                //The code terminates the j-loop early in that case (and also sets
                //j equal to first->nelements to ensure that the remainder loop
                //does not run either.
                AVX_FLOATS m_mask_pimax = AVX_COMPARE_FLOATS(m_zdiff,m_pimax,_CMP_LT_OS);
                if(AVX_TEST_COMPARISON(m_mask_pimax) == 0) {
                    //If the execution reaches here -> then none of the NVEC zdiff values
                    //are smaller than pimax. We can terminate the j-loop now.
                    
                    //set j so that the remainder loop does not run
                    j = first->nelements;
                    //break out of the j-loop
                    break;
                }
                
                /* //Create a mask with true bits when  0 <= zdiff < pimax. */
                m_mask_pimax = AVX_BITWISE_AND(AVX_COMPARE_FLOATS(m_zdiff,m_zero,_CMP_GE_OS),m_mask_pimax);
                const AVX_FLOATS m_rpmax_mask = AVX_COMPARE_FLOATS(r2, m_rupp_sqr[nbin-1], _CMP_LT_OS);
                const AVX_FLOATS m_rpmin_mask = AVX_COMPARE_FLOATS(r2, m_rupp_sqr[0], _CMP_GE_OS);
                const AVX_FLOATS m_rp_mask = AVX_BITWISE_AND(m_rpmax_mask,m_rpmin_mask);
                
                //Create a combined mask by bitwise and of m1 and m_mask_left.
                //This gives us the mask for all sqr_rpmin <= r2 < sqr_rpmax
                m_mask_left = AVX_BITWISE_AND(m_mask_pimax,m_rp_mask);

                
                //If not, continue with the next iteration of j-loop
                if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
                    continue;
                }
                
                //There is some r2 that satisfies sqr_rpmin <= r2 < sqr_rpmax && dz < pimax.
                r2 = AVX_BLEND_FLOATS_WITH_MASK(m_rupp_sqr[nbin-1], r2, m_mask_left);
            }
            
#ifdef OUTPUT_RPAVG
            union_mDperp.m_Dperp = AVX_SQRT_FLOAT(r2);
            AVX_FLOATS m_rpbin = AVX_SET_FLOAT((DOUBLE) 0.0);
#endif
            
            //Loop backwards through nbins. m_mask_left contains all the points that are less than rpmax
            for(int kbin=nbin-1;kbin>=1;kbin--) {
                const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_GE_OS);
                const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
                const int test2  = AVX_TEST_COMPARISON(m_bin_mask);
                npair[kbin] += AVX_BIT_COUNT_INT(test2);
#ifdef OUTPUT_RPAVG
                m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
#endif
                m_mask_left = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_LT_OS);
                const int test3 = AVX_TEST_COMPARISON(m_mask_left);
                if(test3 == 0) {
                    break;
                }
            }
            
#ifdef OUTPUT_RPAVG
            union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
            //protect the unroll pragma in case compiler is not icc.
#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
            for(int jj=0;jj<NVEC;jj++) {
                const int kbin = union_rpbin.ibin[jj];
                const DOUBLE r = union_mDperp.Dperp[jj];
                rpavg[kbin] += r;
            }
#endif//OUTPUT_RPAVG
#endif//end of AVX section               
        }//end of j-loop

        //remainder loop 
        for(;j<first->nelements;j++){
            const DOUBLE dz = *z2++ - z1pos;
            const DOUBLE dx = *x2++ - x1pos;
            const DOUBLE dy = *y2++ - y1pos;
            
            if(dz >= pimax ) {
                break;
            }
            
            const DOUBLE r2 = dx*dx + dy*dy;
            if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
                continue;
            }
            
#ifdef OUTPUT_RPAVG
            const DOUBLE r = SQRT(r2);
#endif

            for(int kbin=nbin-1;kbin>=1;kbin--) {
                if(r2 >= rupp_sqr[kbin-1]) {
                    npair[kbin]++;
#ifdef OUTPUT_RPAVG
                    rpavg[kbin] += r;
#endif
                    break;
                }
            }
        }//remainder loop over cellstruct first
    }//outer loop over cellstruct first

    for(int i=0;i<nbin;i++) {
        src_npairs[i] += npair[i];
#ifdef OUTPUT_RPAVG
        src_rpavg[i]  += rpavg[i];
#endif        
    }
}    

static inline void different_cell_wp_kernel(const cellarray_index * first, const cellarray_index *second,
                                            DOUBLE * restrict const X,DOUBLE * restrict const Y,DOUBLE * restrict const Z,
                                            const DOUBLE sqr_rpmax, const DOUBLE sqr_rpmin, const int nbin, const DOUBLE rupp_sqr[], const DOUBLE pimax,
                                            const DOUBLE off_xwrap, const DOUBLE off_ywrap, const DOUBLE off_zwrap
#ifdef USE_AVX
                                            ,const AVX_FLOATS m_rupp_sqr[] 
#ifdef OUTPUT_RPAVG
                                            ,const AVX_FLOATS m_kbin[]
#endif                                            
#endif                                            

                                            
#ifdef OUTPUT_RPAVG
                                            ,DOUBLE src_rpavg[]
#endif                         
                                            ,uint64_t *src_npairs)
{
    uint64_t npair[nbin];
    for(int i=0;i<nbin;i++) {
        npair[i] = 0;
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nbin];
    for(int i=0;i<nbin;i++) {
        rpavg[i] = ZERO;
    }
#endif    

    DOUBLE *x1 = X + first->start;
    DOUBLE *y1 = Y + first->start;
    DOUBLE *z1 = Z + first->start;
    /* int64_t first_offset = first->start; */
    
    for(int64_t i=0;i<first->nelements;i++) {
        const DOUBLE x1pos = *x1++ + off_xwrap;
        const DOUBLE y1pos = *y1++ + off_ywrap;
        const DOUBLE z1pos = *z1++ + off_zwrap;
        /* first_offset++; */
        
        DOUBLE *x2 = X + second->start;
        DOUBLE *y2 = Y + second->start;
        DOUBLE *z2 = Z + second->start;

        
        int64_t j=0;
        while(j < second->nelements){
            const DOUBLE dz = *z2++-z1pos;
            if(dz > -pimax) break;
            j++;
        }
        z2 -= 1;
        x2 += j;
        y2 += j;
        /* int64_t second_offset = second->start + j; */
        
        for(;j<=(second->nelements-NVEC);j+=NVEC) {
            
#if !(defined(USE_AVX) && defined(__AVX__))
            for(int jj=0;jj<NVEC;jj++) {
                const DOUBLE dx = *x2++ - x1pos;
                const DOUBLE dy = *y2++ - y1pos;
                const DOUBLE dz = *z2++ - z1pos;
                /* second_offset++; */
                
                if(dz >= pimax) {
                    j = second->nelements;
                    break;
                } /* else if((dz == ZERO) && (first_offset < second_offset)) { */
                /*     continue; */
                /* } */
                
                const DOUBLE r2 = dx*dx + dy*dy;
                if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
                    continue;
                }
                
#ifdef OUTPUT_RPAVG
                const DOUBLE r = SQRT(r2);
#endif
                for(int kbin=nbin-1;kbin>=1;kbin--) {
                    if(r2 >= rupp_sqr[kbin-1]) {
                        npair[kbin]++;
#ifdef OUTPUT_RPAVG
                        rpavg[kbin]+=r;
#endif
                        break;
                    }
                }
            }//end of jj-loop
#else //AVX section
            const AVX_FLOATS m_xpos    = AVX_SET_FLOAT(x1pos);
            const AVX_FLOATS m_ypos    = AVX_SET_FLOAT(y1pos);
            const AVX_FLOATS m_zpos    = AVX_SET_FLOAT(z1pos);
            /* const AVX_FLOATS m_first_offset = AVX_SET_FLOAT((DOUBLE) first_offset); */
            /* union off8{ */
            /*     AVX_FLOATS m_offsets; */
            /*     DOUBLE new_offsets[NVEC];     */
            /* }; */
            /* union off8 union_m_offsets; */
            
            /* for(int ii=0;ii<NVEC;ii++) { */
            /*     union_m_offsets.new_offsets[ii] = (DOUBLE) ii; */
            /* } */
            /* const AVX_FLOATS m_second_off = AVX_SET_FLOAT((DOUBLE) second_offset); */
            /* const AVX_FLOATS m_second_offset = AVX_ADD_FLOATS(m_second_off, union_m_offsets.m_offsets); */
            
#ifdef OUTPUT_RPAVG
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
            const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(x2);
            const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(y2);
            const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(z2);
            
            x2 += NVEC;//this might actually exceed the allocated range but we will never dereference that
            y2 += NVEC;
            z2 += NVEC;

            const AVX_FLOATS m_pimax = AVX_SET_FLOAT(pimax);
            /* const AVX_FLOATS m_sqr_pimax = AVX_SET_FLOAT(sqr_pimax); */
            /* const AVX_FLOATS m_zero  = AVX_SET_FLOAT(ZERO); */
            
            
            const AVX_FLOATS m_zdiff = AVX_SUBTRACT_FLOATS(m_z2,m_zpos);//z2[j:j+NVEC-1] - z1
            /* const AVX_FLOATS m_sqr_zdiff = AVX_SQUARE_FLOAT(m_zdiff); */
            const AVX_FLOATS m_sqr_xdiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_xpos,m_x2));//(x0 - x[j])^2
            const AVX_FLOATS m_sqr_ydiff = AVX_SQUARE_FLOAT(AVX_SUBTRACT_FLOATS(m_ypos,m_y2));//(y0 - y[j])^2
            AVX_FLOATS r2  = AVX_ADD_FLOATS(m_sqr_xdiff,m_sqr_ydiff);
            
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
                    //None of the dz^2 values satisfies dz^2 < pimax^2
                    // => no pairs can be added -> continue and process the next NVEC
                    j=second->nelements;
                    break;
                }
                
                /* Create a mask with true bits when  0 <= zdiff < pimax. 
                   Special case of dz == 0.0, then only consider those pairs half the times.
                   Otherwise, double-counting will occur.
                */
                /* This is corresponding condition from the non-avx part of the code 
                   if((sqr_dz == 0.0) && first_offset < second_offset) continue; */
                /* const AVX_FLOATS m_zeros_mask = AVX_COMPARE_FLOATS(m_zdiff, m_zero, _CMP_EQ_OS); */
                /* const AVX_FLOATS m_idx_mask   = AVX_COMPARE_FLOATS(m_first_offset, m_second_offset, _CMP_LT_OS); */
                /* const AVX_FLOATS m_no_task = AVX_BITWISE_AND(m_zeros_mask, m_idx_mask); */

                // NOT m_no_task AND m_mask_pimax. the intrinsic should really be called NOT_AND
                /* m_mask_pimax = AVX_AND_NOT(m_no_task,m_mask_pimax); */
                
                const AVX_FLOATS m_rpmax_mask = AVX_COMPARE_FLOATS(r2, m_rupp_sqr[nbin-1], _CMP_LT_OS);
                const AVX_FLOATS m_rpmin_mask = AVX_COMPARE_FLOATS(r2, m_rupp_sqr[0], _CMP_GE_OS);
                const AVX_FLOATS m_rp_mask = AVX_BITWISE_AND(m_rpmax_mask,m_rpmin_mask);
                
                //Create a combined mask by bitwise and of m1 and m_mask_left.
                //This gives us the mask for all sqr_rpmin <= r2 < sqr_rpmax
                m_mask_left = AVX_BITWISE_AND(m_mask_pimax,m_rp_mask);
                
                //If not, continue with the next iteration of j-loop
                if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
                    continue;
                }
                
                //There is some r2 that satisfies sqr_rpmin <= r2 < sqr_rpmax && 0.0 <= dz^2 < pimax^2.
                r2 = AVX_BLEND_FLOATS_WITH_MASK(m_rupp_sqr[nbin-1], r2, m_mask_left);
            }
            
#ifdef OUTPUT_RPAVG
            union_mDperp.m_Dperp = AVX_SQRT_FLOAT(r2);
            AVX_FLOATS m_rpbin = AVX_SET_FLOAT(ZERO);
#endif
            
            //Loop backwards through nbins. m_mask_left contains all the points that are less than rpmax
            for(int kbin=nbin-1;kbin>=1;kbin--) {
                const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_GE_OS);
                const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
                const int test2  = AVX_TEST_COMPARISON(m_bin_mask);
                npair[kbin] += AVX_BIT_COUNT_INT(test2);
#ifdef OUTPUT_RPAVG
                m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
#endif
                m_mask_left = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_LT_OS);
                const int test3 = AVX_TEST_COMPARISON(m_mask_left);
                if(test3 == 0) {
                    break;
                }
            }
            
#ifdef OUTPUT_RPAVG
            union_rpbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_rpbin);
            //protect the unroll pragma in case compiler is not icc.
#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
            for(int jj=0;jj<NVEC;jj++) {
                const int kbin = union_rpbin.ibin[jj];
                const DOUBLE r = union_mDperp.Dperp[jj];
                rpavg[kbin] += r;
            }
#endif//OUTPUT_RPAVG
            
#endif//end of AVX                                
        }//end of j-loop

        //remainder loop 
        for(;j<second->nelements;j++){
            const DOUBLE dz = *z2++ - z1pos;
            const DOUBLE dx = *x2++ - x1pos;
            const DOUBLE dy = *y2++ - y1pos;
            /* second_offset++; */

            if(dz >= pimax) {
                j = second->nelements;
                break;
            } /* else if((dz == ZERO) && (first_offset < second_offset)) { */
            /*     //different particles might have same z. Only count them once */
            /*     continue; */
            /* } */

            
            const DOUBLE r2 = dx*dx + dy*dy;
            if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
                continue;
            }
            
#ifdef OUTPUT_RPAVG
            const DOUBLE r = SQRT(r2);
#endif

            
            for(int kbin=nbin-1;kbin>=1;kbin--) {
                if(r2 >= rupp_sqr[kbin-1]) {
                    npair[kbin]++;
#ifdef OUTPUT_RPAVG
                    rpavg[kbin] += r;
#endif
                    break;
                }
            }
        }//remainder loop over cellstruct second
    }//loop over cellstruct first

    for(int i=0;i<nbin;i++) {
        src_npairs[i] += npair[i];
#ifdef OUTPUT_RPAVG
        src_rpavg[i]  += rpavg[i];
#endif        
    }
}



results_countpairs_wp *countpairs_wp(const int64_t ND, DOUBLE * restrict X, DOUBLE * restrict Y, DOUBLE * restrict Z,
                                     const double boxsize,
#ifdef USE_OMP
                                     const int numthreads,
#endif
                                     const char *binfile,
                                     const double pimax)
{

    int bin_refine_factor=2,zbin_refine_factor=1;
    int nmesh_x, nmesh_y, nmesh_z;


    /* #ifdef USE_OMP */
    /*  if(numthreads == 1) { */
    /*      bin_refine_factor=2; */
    /*      zbin_refine_factor=1; */
    /*  } else { */
    /*      //seems redundant - but I can not discount the possibility that some */
    /*      //combination of these refine factors will be faster on a different architecture. */
    /*      bin_refine_factor=2; */
    /*      zbin_refine_factor=1; */
    /*  } */
    /* #endif */

    /***********************
     *initializing the  bins
     ************************/
    double *rupp;
    double rpmin,rpmax;
    int nrpbins;
    setup_bins(binfile,&rpmin,&rpmax,&nrpbins,&rupp);
    assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
    assert(nrpbins > 0 && "Number of rp bins is valid");

    const DOUBLE xmin = 0.0, xmax=boxsize;
    const DOUBLE ymin = 0.0, ymax=boxsize;
    const DOUBLE zmin = 0.0, zmax=boxsize;
    DOUBLE rupp_sqr[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        rupp_sqr[i] = rupp[i]*rupp[i];
    }

    const DOUBLE sqr_rpmin = rupp_sqr[0];
    const DOUBLE sqr_rpmax = rupp_sqr[nrpbins-1];
    /* const DOUBLE sqr_pimax = (DOUBLE) pimax*pimax; */

    //set up the 3-d grid structure. Each element of the structure contains a
    //pointer to the cellarray structure that itself contains all the points
    cellarray_index *lattice = gridlink_index(ND, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    if(nmesh_x <= 10 && nmesh_y <= 10 && nmesh_z <= 10) {
        fprintf(stderr,"countpairs_wp> gridlink seems inefficient - boosting bin refine factor - should lead to better performance\n");
        bin_refine_factor *=2;
        zbin_refine_factor *=2;
        const int64_t totncells = nmesh_x*nmesh_y*(int64_t) nmesh_z;
        free_cellarray_index(lattice, totncells);
        lattice = gridlink_index(ND, X, Y, Z, xmin, xmax, ymin, ymax, zmin, zmax, rpmax, rpmax, pimax, bin_refine_factor, bin_refine_factor, zbin_refine_factor, &nmesh_x, &nmesh_y, &nmesh_z);
    }
    const int64_t totncells = nmesh_x*nmesh_y*(int64_t) nmesh_z;


#ifdef USE_OMP
    omp_set_num_threads(numthreads);
#pragma omp parallel for schedule(dynamic)
#endif
    for(int64_t icell=0;icell<totncells;icell++) {
        const cellarray_index *first=&(lattice[icell]);
        if(first->nelements > 0) {
          const int64_t start = first->start;
          DOUBLE *x = X + start;
          DOUBLE *y = Y + start;
          DOUBLE *z = Z + start;
          
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,x,i,j); \
            SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,y,i,j);               \
            SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,z,i,j) }
          
          SGLIB_ARRAY_QUICK_SORT(DOUBLE, z, first->nelements, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER
        }
    }
#if defined(USE_AVX) && defined(__AVX__)
    AVX_FLOATS m_rupp_sqr[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
    }
#ifdef OUTPUT_RPAVG
    AVX_FLOATS m_kbin[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
    }
#endif//RPAVG
#endif


#ifdef USE_OMP
    omp_set_num_threads(numthreads);
    uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, nrpbins);
#ifdef OUTPUT_RPAVG
    DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE), numthreads, nrpbins);
#endif//OUTPUT_RPAVG

#else//USE_OMP
    uint64_t npair[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        npair[i]=0;
    }
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nrpbins];
    for(int i=0;i<nrpbins;i++) {
        rpavg[i]=0;
    }
#endif//OUTPUT_RPAVG
#endif// USE_OMP

    int interrupted=0;
    int64_t numdone=0;
    init_my_progressbar(totncells,&interrupted);

#ifdef USE_OMP
#pragma omp parallel shared(numdone)
    {
        const int tid = omp_get_thread_num();
        uint64_t npair[nrpbins];
        for(int i=0;i<nrpbins;i++) {
            npair[i]=0;
        }
#ifdef OUTPUT_RPAVG
        DOUBLE rpavg[nrpbins];
        for(int i=0;i<nrpbins;i++) {
            rpavg[i]=0.0;
        }
#endif


#pragma omp for schedule(dynamic) nowait
#endif//USE_OMP
        for(int index1=0;index1<totncells;index1++) {

#ifdef USE_OMP
            if (omp_get_thread_num() == 0)
#endif
                my_progressbar(numdone,&interrupted);


#ifdef USE_OMP
#pragma omp atomic
#endif
            numdone++;

            /* First do the same-cell calculations */
            const cellarray_index *first  = &(lattice[index1]);
            if(first->nelements == 0) {
                continue;
            }
            
            same_cell_wp_kernel(first, X, Y, Z, sqr_rpmax, sqr_rpmin, nrpbins, rupp_sqr, pimax
#ifdef USE_AVX
                                ,m_rupp_sqr
#ifdef OUTPUT_RPAVG
                                ,m_kbin
#endif
#endif
                                
#ifdef OUTPUT_RPAVG
                                ,rpavg
#endif
                                ,npair);

            for(int64_t ngb=0;ngb<first->num_ngb;ngb++){
                cellarray_index *second = first->ngb_cells[ngb];
                const DOUBLE off_xwrap = first->xwrap[ngb];
                const DOUBLE off_ywrap = first->ywrap[ngb];
                const DOUBLE off_zwrap = first->zwrap[ngb];
                different_cell_wp_kernel(first, second,X, Y, Z,sqr_rpmax, sqr_rpmin, nrpbins, rupp_sqr, pimax, off_xwrap, off_ywrap, off_zwrap
#ifdef USE_AVX
                                         ,m_rupp_sqr
#ifdef OUTPUT_RPAVG
                                         ,m_kbin
#endif                                                 
#endif                                                 
                                         
#ifdef OUTPUT_RPAVG
                                         ,rpavg
#endif
                                         ,npair);
                
            }//ngb loop
        }//index1 loop
#ifdef USE_OMP
        for(int j=0;j<nrpbins;j++) {
            all_npairs[tid][j] = npair[j];
#ifdef OUTPUT_RPAVG
            all_rpavg[tid][j] = rpavg[j];
#endif
        }
    }//omp parallel
#endif
    finish_myprogressbar(&interrupted);


#ifdef USE_OMP
    uint64_t npair[nrpbins];
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nrpbins];
#endif
    for(int i=0;i<nrpbins;i++) {
        npair[i] = 0;
#ifdef OUTPUT_RPAVG
        rpavg[i] = 0.0;
#endif
    }

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<nrpbins;j++) {
            npair[j] += all_npairs[i][j];
#ifdef OUTPUT_RPAVG
            rpavg[j] += all_rpavg[i][j];
#endif
        }
    }
    matrix_free((void **) all_npairs,numthreads);
#ifdef OUTPUT_RPAVG
    matrix_free((void **) all_rpavg, numthreads);
#endif //OUTPUT_RPAVG
#endif//USE_OMP

#ifdef OUTPUT_RPAVG
    for(int i=0;i<nrpbins;i++) {
        if(npair[i] > 0) {
            rpavg[i] /= (DOUBLE) npair[i];
        }
    }
#endif

    free_cellarray_index(lattice, totncells);

    //Pack in the results
    results_countpairs_wp *results = my_malloc(sizeof(*results), 1);
    results->nbin  = nrpbins;
    results->pimax = pimax;
    results->npairs = my_malloc(sizeof(uint64_t), nrpbins);
    results->wp = my_malloc(sizeof(DOUBLE), nrpbins);
    results->rupp   = my_malloc(sizeof(DOUBLE), nrpbins);
    results->rpavg  = my_malloc(sizeof(DOUBLE), nrpbins);


    const DOUBLE avgweight2 = 1.0, avgweight1 = 1.0;
    const DOUBLE density=0.5*avgweight2*ND/(boxsize*boxsize*boxsize);//pairs are not double-counted
    DOUBLE rlow=0.0 ;
    DOUBLE prefac_density_DD=avgweight1*ND*density;
    DOUBLE twice_pimax = 2.0*pimax;

    for(int i=0;i<nrpbins;i++) {
        results->npairs[i] = npair[i];
        results->rupp[i] = rupp[i];
#ifdef OUTPUT_RPAVG
        results->rpavg[i] = rpavg[i];
#else
        results->rpavg[i] = 0.0;
#endif
        const DOUBLE weight0 = (DOUBLE) results->npairs[i];
        /* compute xi, dividing summed weight by that expected for a random set */
        const DOUBLE vol=M_PI*(results->rupp[i]*results->rupp[i]-rlow*rlow)*twice_pimax;
        const DOUBLE weightrandom = prefac_density_DD*vol;
        results->wp[i] = (weight0/weightrandom-1)*twice_pimax;
        rlow=results->rupp[i];
    }
    free(rupp);
    return results;
}
