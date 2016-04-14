/* File: countpairs_theta_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

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

#include "defs.h"
#include "gridlink_mocks.h"//function proto-type for gridlink
#include "countpairs_theta_mocks.h" //function proto-type
#include "cellarray_mocks.h" //definition of struct cellarray_mocks
#include "utils.h" //all of the utilities

#ifndef SILENT
#include "progressbar.h" //for the progressbar
#endif

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif

#if defined(USE_OMP) && defined(_OPENMP)
#include <omp.h>
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
        }
        assert(theta[i] <= 180.0 && "Declination should not be more than 180. Did you swap ra and dec?");
        if(theta[i] > 90.0) {
            fix_dec = 1;
        }
    }

    if(fix_ra == 1 || fix_dec == 1) {
        if(fix_ra == 1) {
            fprintf(stderr,ANSI_COLOR_YELLOW "DDtheta> Out of range values found for ra. Expected ra to be in the range [0.0,360.0]. Found ra values in [-180,180] -- fixing that" ANSI_COLOR_RESET "\n");
        }
        if(fix_dec == 1) {
            fprintf(stderr,ANSI_COLOR_YELLOW "DDtheta> Out of range values found for dec. Expected dec to be in the range [-90.0,90.0]. Found dec values in [0,180] -- fixing that" ANSI_COLOR_RESET "\n");
        }

        for(int64_t i=0;i<N;i++) {
            if(fix_ra == 1) {
                phi[i] += 180.0;
            }
            if(fix_dec == 1) {
                theta[i] -= 90.0;
            }
        }
    }
}


results_countpairs_theta * countpairs_theta_mocks(const int64_t ND1, DOUBLE *phi1, DOUBLE *theta1,
                                                  const int64_t ND2, DOUBLE *phi2, DOUBLE *theta2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                  const int numthreads,
#endif
                                                  const int autocorr,
                                                  const char *binfile)
{

#if defined(LINK_IN_RA) && !defined(LINK_IN_DEC)
#error LINK_IN_DEC Makefile option must be enabled before LINK_IN_RA is selected
#endif

    //check the ra-dec inputs
    check_ra_dec(ND1, phi1, theta1);
    if(autocorr==0) {
        check_ra_dec(ND2, phi2, theta2);
    }

    /* #if !defined(LINK_IN_DEC) && !defined(LINK_IN_RA) */
    /*  { */
    /* #error The brute force code does not work */
    /*      results_countpairs_theta *results = countpairs_theta_brute_force(ND1,phi1,theta1, */
    /*                                                                                                                                       ND2,phi2,theta2, */
    /* #if defined(USE_OMP) && defined(_OPENMP) */
    /*                                                                                                                                       numthreads, */
    /* #endif */
    /*                                                                                                                                       autocorr, */
    /*                                                                                                                                       binfile) ; */
    /*      return results; */
    /*  } */
    /* #endif */



    double *theta_upp;
    int nthetabin;
    double thetamin,thetamax;
    setup_bins(binfile,&thetamin,&thetamax,&nthetabin,&theta_upp);
    assert(thetamin > 0.0 && thetamax > 0.0 && thetamin < thetamax && thetamax < 180.0 &&  "[thetamin, thetamax] are valid inputs");
    assert(nthetabin > 0 && "Number of theta bins is valid");

#ifdef LINK_IN_DEC
    int rbin_refine_factor=2;
#if defined(USE_OMP) && !defined(LINK_IN_RA)
    if(rbin_refine_factor < numthreads) {
        rbin_refine_factor=numthreads;
    }
#endif

#ifdef LINK_IN_RA
    //roughly occurs around
    int phi_bin_refine_factor=2;
    if(thetamax >= 25) {
        fprintf(stderr,ANSI_COLOR_YELLOW "LINK_IN_RA can produce incorrect answers if the angular limits are large. Please cross-check with the output where LINK_IN_RA is not defined. " ANSI_COLOR_RESET"\n");
        fprintf(stderr,ANSI_COLOR_YELLOW "If increasing phi_bin_refine_factor changes the answer -- that is a good indication that the calculation is incorrect" ANSI_COLOR_RESET "\n");
        phi_bin_refine_factor=1;//do not change this line. Increases the chance of incorrect answers.
    }

#endif
#endif




#if !(defined(USE_OMP) && defined(_OPENMP))
    uint64_t npairs[nthetabin];
    for(int i=0;i<nthetabin;i++) npairs[i]=0;
#else
    assert(numthreads >=1 && "Number of requested threads must be at least 1");
    omp_set_num_threads(numthreads);
    uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, nthetabin);
#endif

    DOUBLE costheta_upp[nthetabin];
    for(int i=0;i<nthetabin;i++) {
        costheta_upp[i] = COSD(theta_upp[i]);
    }
    const DOUBLE costhetamin=costheta_upp[0];
    const DOUBLE costhetamax=costheta_upp[nthetabin-1];

#ifdef OUTPUT_THETAAVG
#if !(defined(USE_OMP) && defined(_OPENMP))
    DOUBLE thetaavg[nthetabin];
    for(int i=0;i<nthetabin;i++) thetaavg[i] = 0.0;
#else
    DOUBLE **all_thetaavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE),numthreads,nthetabin);
#endif
#endif

#if defined(USE_AVX) && defined(__AVX__)
    AVX_FLOATS m_costheta_upp[nthetabin] ;
    for(int i=0;i<nthetabin;i++) {
        /* fprintf(stderr," i = %d theta_upp[i-1] = %lf cos(theta_upp[i-1] = %lf cos(theta_upp[i]) = %lf \n",i, theta_upp[i-1],COSD(theta_upp[i-1]),COSD(theta_upp[i])); */
        m_costheta_upp[i] = AVX_SET_FLOAT(costheta_upp[i]);
    }
    /* const AVX_FLOATS m_costhetamin=AVX_SET_FLOAT(costhetamin); */
    const AVX_FLOATS m_costhetamax=AVX_SET_FLOAT(costhetamax);
#ifdef OUTPUT_THETAAVG
    AVX_FLOATS m_kbin[nthetabin];
    for(int i=0;i<nthetabin;i++) {
        m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
    }
#endif
#endif

    /*---Prepare-Data2--------------------------------*/
#ifdef LINK_IN_DEC
    double dec_min=90.0,dec_max=-90.0;
#endif

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
    }

    DOUBLE *x1,*y1,*z1;

    if (autocorr==0) {
        x1 = my_malloc(sizeof(*x1),ND1);
        y1 = my_malloc(sizeof(*y1),ND1);
        z1 = my_malloc(sizeof(*z1),ND1);

        for(int i=0;i<ND1;i++) {
            x1[i] = COSD(theta1[i])*COSD(phi1[i]) ;
            y1[i] = COSD(theta1[i])*SIND(phi1[i]) ;
            z1[i] = SIND(theta1[i]);

#ifdef LINK_IN_DEC
            if(theta1[i] < dec_min)
                dec_min = theta1[i];
            if(theta1[i] > dec_max)
                dec_max = theta1[i];
#endif
        }
    } else {
        x1 = x2;
        y1 = y2;
        z1 = z2;
    }

#ifdef LINK_IN_DEC
    DOUBLE dec_diff = dec_max-dec_min;
    DOUBLE inv_dec_diff=1.0/dec_diff;
    int ngrid_dec=0,max_n=0;
#ifndef LINK_IN_RA

    cellarray *lattice2 = gridlink1D_theta(ND2,
                                           dec_min, dec_max, thetamax,
                                           x2, y2, z2,theta2,
                                           &ngrid_dec,
                                           &max_n,
                                           rbin_refine_factor);
    /* cellarray *lattice1; */
    /* int ngrid_dec1; */
    /* if(autocorr==0) { */
    /*  lattice1 = gridlink1D_theta(ND1, */
    /*                                                          dec_min, dec_max, thetamax, */
    /*                                                          x1, y1, z1,theta1, */
    /*                                                          &ngrid_dec1, */
    /*                                                          &max_n, */
    /*                                                          rbin_refine_factor); */
    /*  assert(ngrid_dec1 == ngrid_dec && "The two lattices should have identical dec-bins"); */
    /* } else { */
    /*  lattice1 = lattice2; */
    /* } */

#else
    int *ngrid_ra=NULL;
    const DOUBLE ra_min=0.0,ra_max=360.0;
    const DOUBLE inv_ra_diff=1.0/(ra_max-ra_min);
    cellarray **lattice2 = gridlink2D_theta(ND2, dec_min, dec_max, thetamax,
                                            x2, y2, z2,
                                            theta2,
                                            &ngrid_dec,
                                            phi2,ra_min,ra_max,
                                            &ngrid_ra,
                                            &max_n,
                                            rbin_refine_factor,
                                            phi_bin_refine_factor);

    /* cellarray **lattice1; */
    /* int *ngrid_ra1 = NULL; */
    /* int ngrid_dec1 = 0; */
    /* if(autocorr==0) { */
    /*  lattice1 = gridlink2D_theta(ND1, dec_min, dec_max, thetamax, */
    /*                                                          x1, y1, z1, */
    /*                                                          theta1, */
    /*                                                          &ngrid_dec1, */
    /*                                                          phi1,ra_min,ra_max, */
    /*                                                          &ngrid_ra1, */
    /*                                                          &max_n, */
    /*                                                          rbin_refine_factor, */
    /*                                                          phi_bin_refine_factor); */
    /*  assert(ngrid_dec1 == ngrid_dec && "The two lattices should have identical dec-bins"); */
    /*  for(int i=0;i<ngrid_dec1;i++) { */
    /*      assert(ngrid_ra1[i] == ngrid_ra[i] && "The two lattices should have identical ra-bins"); */
    /*  } */
    /*  free(ngrid_ra1); */
    /* } else { */
    /*  lattice1 = lattice2; */
    /* } */


#endif
    if(autocorr == 0) {
        free(x2);free(y2);free(z2);
    }

    /* if(autocorr==0) { */
    /*  free(x1);free(y1);free(z1); */
    /* } */
#endif


#ifndef SILENT    
    int interrupted=0;
    int numdone=0;
    init_my_progressbar(ND1,&interrupted);
#endif    

    /*---Loop-over-Data1-particles--------------------*/
#if defined(USE_OMP) && defined(_OPENMP)
#if !defined(LINK_IN_DEC ) && !defined(LINK_IN_RA)
#pragma omp parallel shared(x2,y2,z2)
#else
#pragma omp parallel private(x2,y2,z2)
#endif
    {
        int tid = omp_get_thread_num();
        uint64_t npairs[nthetabin] __attribute__ ((aligned (ALIGNMENT)));
        for(int i=0;i<nthetabin;i++) npairs[i] = 0;
#ifdef OUTPUT_THETAAVG
        DOUBLE thetaavg[nthetabin] __attribute__ ((aligned (ALIGNMENT)));
        for(int i=0; i < nthetabin;i++) thetaavg[i] = 0.0;
#endif


#pragma omp for  schedule(dynamic)
#endif
        for(int i=0;i<ND1;i++) {

#ifndef SILENT          
#if defined(USE_OMP) && defined(_OPENMP)
            if (omp_get_thread_num() == 0)
#endif
                my_progressbar(numdone,&interrupted);

#if defined(USE_OMP) && defined(_OPENMP)
#pragma omp atomic
#endif
            numdone++;
#endif//SILENT            


            const DOUBLE x1pos = x1[i];
            const DOUBLE y1pos = y1[i];
            const DOUBLE z1pos = z1[i];

#ifdef LINK_IN_DEC
            const DOUBLE decpos = theta1[i];
            int dec_iz = (int)(ngrid_dec*(decpos-dec_min)*inv_dec_diff);
            if (dec_iz >= ngrid_dec) dec_iz--;
            assert(dec_iz >= 0 && dec_iz < ngrid_dec && "Declination is within lattice2 bounds");
            /* int dec_limits = (int) (ceil(thetamax*inv_dec_diff*ngrid_dec)); */
            /* const int dec_limits = rbin_refine_factor; */
            /* const int min_dec = (dec_iz - dec_limits) <= 0 ? 0:dec_iz - dec_limits; */
            /* const int max_dec = (dec_iz + dec_limits) >= (ngrid_dec-1) ? (ngrid_dec-1):dec_iz + dec_limits; */
            for(int iidec=-rbin_refine_factor;iidec<=rbin_refine_factor;iidec++) {
                const int idec = dec_iz + iidec;
                if(idec < 0 || idec >= ngrid_dec) continue;

#ifndef LINK_IN_RA
                const cellarray *second = &(lattice2[idec]);
#else
                const DOUBLE rapos = phi1[i];
                int ra_iz = (int)(ngrid_ra[idec]*(rapos-ra_min)*inv_ra_diff);
                if (ra_iz >= ngrid_ra[idec]) ra_iz--;
                assert(ra_iz >= 0 && ra_iz < ngrid_ra[idec] && "RA position is within bounds");
                for(int ira_step=-phi_bin_refine_factor;ira_step<=phi_bin_refine_factor;ira_step++) {
                    const int ira = (ra_iz + ira_step + ngrid_ra[idec]) % ngrid_ra[idec];
                    const cellarray *second = &(lattice2[idec][ira]);
#endif //LINK_IN_DEC

                    x2 = second->pos ;
                    y2 = second->pos + NVEC;
                    z2 = second->pos + 2*NVEC;
                    const int Nloop = second->nelements;
#else //No linking in RA or DEC
                    const int Nloop = ND2;
#endif


                    DOUBLE *localx2 = x2;
                    DOUBLE *localy2 = y2;
                    DOUBLE *localz2 = z2;

                    /*---Loop-over-Data2-particles--------------------*/
                    int j;
                    for(j=0;j <=(Nloop-NVEC);j+=NVEC) {
#if !(defined(USE_AVX) && defined(__AVX__))
                        DOUBLE costheta[NVEC];
                        int thetabin[NVEC];
#ifdef OUTPUT_THETAAVG
                        DOUBLE theta[NVEC];
#endif
                        int num_bad=0;
                        for(int k=0;k<NVEC;k++) {
                            const DOUBLE tmp = x1pos*localx2[k] + y1pos*localy2[k] + z1pos*localz2[k];
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
#ifdef LINK_IN_DEC
                        localx2 += 3*NVEC;
                        localy2 += 3*NVEC;
                        localz2 += 3*NVEC;
#else
                        localx2 += NVEC;
                        localy2 += NVEC;
                        localz2 += NVEC;
#endif


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

#else //USE_AVX
                        const AVX_FLOATS m_x1 = AVX_SET_FLOAT(x1pos);
                        const AVX_FLOATS m_y1 = AVX_SET_FLOAT(y1pos);
                        const AVX_FLOATS m_z1 = AVX_SET_FLOAT(z1pos);


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
                        const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(localx2);
                        const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(localy2);
                        const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(localz2);
#ifdef LINK_IN_DEC
                        localx2 += 3*NVEC;
                        localy2 += 3*NVEC;
                        localz2 += 3*NVEC;
#else
                        localx2 += NVEC;
                        localy2 += NVEC;
                        localz2 += NVEC;
#endif

                        const AVX_FLOATS m_tmp1 = AVX_MULTIPLY_FLOATS(m_x2,m_x1);
                        const AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_y2,m_y1);
                        const AVX_FLOATS m_tmp3 = AVX_MULTIPLY_FLOATS(m_z2,m_z1);
                        const AVX_FLOATS m_costheta = AVX_ADD_FLOATS(m_tmp1,AVX_ADD_FLOATS(m_tmp2,m_tmp3));

                        AVX_FLOATS m_mask_left = AVX_COMPARE_FLOATS(m_costheta,m_costhetamax,_CMP_GT_OS);
#ifdef DOUBLE_PREC
                        {
                            //Only seems to help if double precision is enabled -> wrapping it within the double-prec flags
                            const AVX_FLOATS m_costhetamin = AVX_SET_FLOAT(costhetamin);
                            AVX_FLOATS m1 = AVX_COMPARE_FLOATS(m_costheta,m_costhetamin,_CMP_LE_OS);
                            AVX_FLOATS m_mask = AVX_BITWISE_AND(m1,m_mask_left);
                            if(AVX_TEST_COMPARISON(m_mask) == 0) {
                                continue;
                            }
                        }
#endif

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
#ifdef OUTPUT_THETAAVG
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
                    for(int ipos=0;j<Nloop;j++,ipos++) {
                        const DOUBLE costheta = x1pos*localx2[ipos] + y1pos*localy2[ipos] + z1pos*localz2[ipos] ;
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
                }//finish the loop over ra-cells in second
#endif
            }//finish the loop over dec-cells in second
#endif//LINK_IN_DEC

        }//loop over i
#if defined(USE_OMP) && defined(_OPENMP)
        for(int j=0;j<nthetabin;j++) {
            all_npairs[tid][j] = npairs[j];
        }
#ifdef OUTPUT_THETAAVG
        for(int j=0;j<nthetabin;j++) {
            all_thetaavg[tid][j] = thetaavg[j];
        }
#endif

    }//close the omp parallel region
#endif

#ifndef SILENT    
    finish_myprogressbar(&interrupted);
#endif    

#if defined(USE_OMP) && defined(_OPENMP)
    uint64_t npairs[nthetabin];
    for(int i=0;i<nthetabin;i++) npairs[i] = 0;

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<nthetabin;j++) {
            npairs[j] += all_npairs[i][j];
        }
    }
#ifdef OUTPUT_THETAAVG
    DOUBLE thetaavg[nthetabin];
    for(int i=0;i<nthetabin;i++) thetaavg[i] = 0.0;

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<nthetabin;j++) {
            thetaavg[j] += all_thetaavg[i][j];
        }
    }
#endif//OUTPUT_THETAAVG
#endif//USE_OMP


#ifndef LINK_IN_DEC
    free(x2);free(y2);free(z2);
    if(autocorr==0) {
        free(x1);free(y1);free(z1);
    }
#else
#ifndef LINK_IN_RA
    for(int i=0;i<ngrid_dec;i++) {
        free(lattice2[i].pos);
        /* if(autocorr==0) { */
        /*  free(lattice1[i].pos); */
        /* } */
    }
    free(lattice2);
    /* if(autocorr==0) { */
    /*  free(lattice1); */
    /* } */
#else
    for(int i=0;i<ngrid_dec;i++) {
        for(int j=0;j<ngrid_ra[i];j++) {
            free(lattice2[i][j].pos);
            /* if(autocorr==0) { */
            /*  free(lattice1[i][j].pos); */
            /* } */
        }
    }
    free(ngrid_ra);
    matrix_free((void **) lattice2,ngrid_dec);
    /* if(autocorr==0) { */
    /*  matrix_free((void **) lattice1,ngrid_dec); */
    /* } */
#endif//LINK_IN_RA
#endif//LINK_IN_DEC


#if defined(USE_OMP) && defined(_OPENMP)
    matrix_free((void **) all_npairs, numthreads);
#ifdef OUTPUT_THETAVG
    matrix_free((void **) all_thetaavg, numthreads);
#endif
#endif


#ifdef OUTPUT_THETAAVG
    for(int i=1;i<nthetabin;i++) {
        if(npairs[i] > 0) {
            thetaavg[i] /= (DOUBLE) npairs[i];
        }
    }
#endif

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
