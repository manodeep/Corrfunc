/* File: countpairs_rp_pi_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_interp.h>

#include "defs.h"
/* #include "sglib.h" */
#include "utils.h"
#include "cellarray_mocks.h"
#include "gridlink_mocks.h"
#include "countpairs_rp_pi_mocks.h"
#include "cosmology_params.h"
#include "set_cosmo_dist.h"

#ifndef SILENT
#include "progressbar.h"
#endif

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif

#if defined(USE_OMP) && defined(_OPENMP)
#include <omp.h>
#endif


void free_results_mocks(results_countpairs_mocks **results)
{
    if(results==NULL)
        return;

    if(*results==NULL)
        return;

    results_countpairs_mocks *tmp = *results;

    free(tmp->npairs);
    free(tmp->rupp);
    free(tmp->rpavg);
    free(tmp);
    tmp = NULL;
}



void check_ra_dec_cz(const int64_t N, DOUBLE *phi, DOUBLE *theta, DOUBLE *cz)
{

    assert(N > 0 && "Number of data-points must be non-zero");
    assert(phi != NULL && theta != NULL && cz != NULL && "Input arrays can not be NULL");
    /* fprintf(stderr,"N = %"PRId64" phi = %p theta = %p cz = %p\n",N,phi,theta,cz); */
    int fix_cz  = 0;
    int fix_ra  = 0;
    int fix_dec = 0;

    const DOUBLE max_cz_threshold = 10.0;//if I find that max cz is smaller than this threshold, then I will assume z has been supplied rather than cz
    DOUBLE max_cz = 0.0;
    //Check input cz -> ensure that cz contains cz and not z
    for(int64_t i=0;i<N;i++) {
        if(cz[i] > max_cz) max_cz = cz[i];
        if(phi[i] < 0.0) {
            fix_ra = 1;
        }
        if(theta[i] > 90.0) {
            fix_dec = 1;
        }
        assert(theta[i] <= 180 && "Declination can not be more than 180 deg. Did you switch the ra/dec variables?");
    }
    if(max_cz < max_cz_threshold) fix_cz = 1;

    //Only run the loop if something needs to be fixed
    if(fix_cz==1 || fix_ra == 1 || fix_dec == 1) {
        if(fix_ra == 1) {
            fprintf(stderr,"countpairs_mocks> Out of range values found for ra. Expected ra to be in the range [0.0,360.0]. Found ra values in [-180,180] -- fixing that\n");
        }
        if(fix_dec == 1) {
            fprintf(stderr,"countpairs_mocks> Out of range values found for dec. Expected dec to be in the range [-90.0,90.0]. Found dec values in [0,180] -- fixing that\n");
        }
        if(fix_cz == 1)  {
            fprintf(stderr,"countpairs_mocks> Out of range values found for cz. Expected input to be `cz' but found `z' instead. max_cz (found in input) = %"DOUBLE_FORMAT" threshold = %"DOUBLE_FORMAT"\n",max_cz,max_cz_threshold);
        }

        for(int64_t i=0;i<N;i++) {
            if(fix_ra==1) {
                phi[i] += 180.0;
            }
            if(fix_dec==1) {
                theta[i] -= 90.0;
            }
            if(fix_cz == 1) {
                cz[i] *= SPEED_OF_LIGHT;//input was z -> convert to cz
            }
        }
    }
}


results_countpairs_mocks * countpairs_mocks(const int64_t ND1, DOUBLE *phi1, DOUBLE *theta1, DOUBLE *czD1,
                                            const int64_t ND2, DOUBLE *phi2, DOUBLE *theta2, DOUBLE *czD2,
#if defined(USE_OMP) && defined(_OPENMP)
                                            const int numthreads,
#endif
                                            const int autocorr,
                                            const char *binfile,
                                            const DOUBLE pimax,
                                            const int cosmology)
{
    /* DOUBLE logrpmax,logrpmin,dlogrp,inv_dlogrp; */
    DOUBLE dpi,inv_dpi;

    int zbin_refine_factor=2;
#if defined(USE_OMP) && defined(_OPENMP)
    if(zbin_refine_factor < numthreads/2)
        zbin_refine_factor=numthreads/2;
#endif

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

    //Check inputs
    assert(ND1 > 0 && "Number of data-points must be non-zero");
    assert(phi1 != NULL && theta1 != NULL && czD1 != NULL && "Input arrays can not be NULL");
    if(autocorr==0) {
        assert(ND2 > 0 && "Number of data-points in data-set 2 must be non-zero");
        assert(phi2 != NULL && theta2 != NULL && czD2 != NULL && "Input arrays for data-set 2 can not be NULL");
    }

    //Try to initialize cosmology - code will exit if comoslogy is not implemented.
    init_cosmology(cosmology);

    //Check inputs
    check_ra_dec_cz(ND1, phi1, theta1, czD1);
    if(autocorr==0) {
        check_ra_dec_cz(ND2, phi2, theta2, czD2);
    }


    /***********************
     *initializing the  bins
     ************************/
    double *rupp;
    int nrpbin ;
    double rpmin,rpmax;
    setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
    assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
    assert(nrpbin > 0 && "Number of rp bins is valid");

    const DOUBLE sqr_max_sep = rpmax*rpmax + pimax*pimax;

    //Change cz into co-moving distance
    DOUBLE *d1     = my_malloc(sizeof(*d1),ND1);
    DOUBLE *d2 = NULL;
    if(autocorr==0) {
        d2     = my_malloc(sizeof(*d2),ND2);
    } else {
        d2 = d1;
    }

    {
        //Setup variables to do the cz->comoving distance
        int Nzdc;
        double *zz,*ddc;
        zz=my_calloc(sizeof(*zz),COSMO_DIST_SIZE);
        ddc=my_calloc(sizeof(*ddc),COSMO_DIST_SIZE);
        Nzdc = set_cosmo_dist(MAX_REDSHIFT_FOR_COSMO_DIST, COSMO_DIST_SIZE, zz, ddc, cosmology);

        gsl_interp *interpolation;
        gsl_interp_accel *accelerator;
        const DOUBLE inv_speed_of_light = 1.0/SPEED_OF_LIGHT;
        accelerator =  gsl_interp_accel_alloc();
        interpolation = gsl_interp_alloc (gsl_interp_linear,Nzdc);
        gsl_interp_init(interpolation, zz, ddc, Nzdc);
        for(int64_t i=0;i<ND1;i++) {
            d1[i] = gsl_interp_eval(interpolation, zz, ddc, czD1[i]*inv_speed_of_light, accelerator);
        }

        if(autocorr==0) {
            for(int64_t i=0;i<ND2;i++) {
                d2[i] = gsl_interp_eval(interpolation, zz, ddc, czD2[i]*inv_speed_of_light, accelerator);
            }
        }
        free(zz);free(ddc);
        gsl_interp_free(interpolation);
        gsl_interp_accel_free(accelerator);
    }

    //Now I am going to sort the d1 array into increasing order.
    /* #define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,theta1,i,j); \ */
    /*      SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,d1,i,j);                                          \ */
    /*      SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,phi1,i,j) } */

    /*  SGLIB_ARRAY_QUICK_SORT(DOUBLE, d1, ND1, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER); */


    /*---Gridlink-variables----------------*/
    int ngrid;
    const int totnbins = (nrpbin+1)*(npibin+1);
#if !(defined(USE_OMP) && defined(_OPENMP))
    uint64_t npairs[totnbins];
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[totnbins];
#endif
    for(int i=0; i <totnbins;i++) {
        npairs[i] = 0;
#ifdef OUTPUT_RPAVG
        rpavg[i] = 0.0;
#endif
    }
#else //USE_OMP
    omp_set_num_threads(numthreads);
    uint64_t **all_npairs = (uint64_t **) matrix_calloc(sizeof(uint64_t), numthreads, totnbins);
#ifdef OUTPUT_RPAVG
    DOUBLE **all_rpavg = (DOUBLE **) matrix_calloc(sizeof(DOUBLE),numthreads,totnbins);
#endif
#endif //USE_OMP

    DOUBLE sqr_rpmin = rpmin*rpmin;
    DOUBLE sqr_rpmax = rpmax*rpmax;
    DOUBLE sqr_pimax = pimax*pimax;

    dpi = pimax/(DOUBLE)npibin ;
    inv_dpi = 1.0/dpi;
    DOUBLE rupp_sqr[nrpbin];
    rupp_sqr[0] = sqr_rpmin;
    for(int i=0;i<nrpbin;i++) {
        rupp_sqr[i] = rupp[i]*rupp[i];
    }

#if defined(USE_AVX) && defined(__AVX__)
    AVX_FLOATS m_rupp_sqr[nrpbin];
    AVX_FLOATS m_kbin[nrpbin];
    for(int i=0;i<nrpbin;i++) {
        m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
        m_kbin[i] = AVX_SET_FLOAT((DOUBLE) i);
    }
#endif


    /*---Prepare-Data2--------------------------------*/
#ifdef LINK_IN_DEC
    DOUBLE dec_min=90.0,dec_max=-90.0;
#endif


    DOUBLE d2min=1000. ;
    DOUBLE d2max=0. ;

    get_max_min_data(ND1, d1, &d2min, &d2max
#ifdef LINK_IN_DEC
                     ,theta1,&dec_min,&dec_max
#endif
                     );

    if(autocorr==0) {
        get_max_min_data(ND2, d2, &d2min, &d2max
#ifdef LINK_IN_DEC
                         ,theta2,&dec_min,&dec_max
#endif
                         );
    }

    ngrid=0;
    int max_n;

#ifdef LINK_IN_DEC
    int *ngrid_dec;
    const DOUBLE dec_diff = dec_max - dec_min;
    const DOUBLE inv_dec_diff=1.0/dec_diff;
    /* fprintf(stderr,"dec_min = %lf dec_max = %lf\n",dec_min,dec_max); */
#ifndef LINK_IN_RA

    cellarray_mocks **lattice2 = gridlink2D(ND2,d2min,d2max,pimax,
                                            dec_min,dec_max,rpmax,
                                            d2, theta2, phi2,
                                            &ngrid, &ngrid_dec, &max_n,
                                            rbin_refine_factor,
                                            zbin_refine_factor);

#else
    //Linking in cz, Dec, RA
    const DOUBLE ra_max=360.0,ra_min=0.0;
    /* const DOUBLE ra_max = 267.0, ra_min = 109.0; */
    const DOUBLE inv_ra_diff=1.0/(ra_max-ra_min);
    int **ngrid_ra=NULL;
    /* fprintf(stderr,"ra_min = %lf ra_max = %lf\n",ra_min,ra_max); */
    cellarray_mocks ***lattice2 = gridlink3D(ND2,d2min,d2max,pimax,
                                             dec_min,dec_max,rpmax,
                                             d2, theta2, phi2,
                                             &ngrid,
                                             &ngrid_dec,
                                             ra_min,ra_max,
                                             &ngrid_ra,
                                             &max_n,
                                             phibin_refine_factor,
                                             rbin_refine_factor,
                                             zbin_refine_factor);

#endif
    //Need cz_binsize for LINK_IN_DEC option
    /* const DOUBLE cz_binsize=(d2max-d2min)/ngrid; */
#else
    //Only linking in cz
    cellarray_mocks *lattice2 = gridlink1D(ND2, d2min, d2max, pimax, theta2, phi2, d2, &ngrid, &max_n,zbin_refine_factor);
#endif

    /* const DOUBLE cz_binsize=(d2max-d2min)/ngrid; */
    const DOUBLE inv_cz_binsize=ngrid/(d2max-d2min);

#ifndef SILENT    
    int interrupted=0,numdone=0;
    init_my_progressbar(ND1,&interrupted);
#endif    

#if defined(USE_OMP) && defined(_OPENMP)
#ifndef SILENT    
#pragma omp parallel shared(numdone)
#else
#pragma omp parallel
#endif//SILENT    
    {
        const int tid = omp_get_thread_num();
        uint64_t npairs[totnbins] __attribute__ ((aligned(ALIGNMENT)));
        for(int i=0;i<totnbins;i++) npairs[i] = 0;
#ifdef OUTPUT_RPAVG
        DOUBLE rpavg[totnbins] __attribute__ ((aligned(ALIGNMENT)));
        for(int i=0;i<totnbins;i++) rpavg[i] = 0.0;
#endif

#pragma omp for  schedule(dynamic)
#endif//USE_OMP
        /*---Loop-over-Data1-particles--------------------*/
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
            
            const DOUBLE x1 = d1[i]*COSD(theta1[i])*COSD(phi1[i]) ;
            const DOUBLE y1 = d1[i]*COSD(theta1[i])*SIND(phi1[i]) ;
            const DOUBLE z1 = d1[i]*SIND(theta1[i]) ;

            /*---Deterpmine-central-grid-cell-of-search--------*/
            int icen = (int)((d1[i]-d2min)*inv_cz_binsize) ;
            if(icen<0) icen = 0 ;
            if(icen>=ngrid) icen = ngrid-1 ;

            const int min_iz = (icen - zbin_refine_factor) <= 0 ? 0:icen - zbin_refine_factor;
            const int max_iz = (icen + zbin_refine_factor) >= (ngrid-1) ? (ngrid-1):icen + zbin_refine_factor;

            /*---Loop-over-surrounding-cells----------------*/
            for(int icell=min_iz;icell<=max_iz;icell++) {
#ifdef LINK_IN_DEC

                DOUBLE decpos = theta1[i];
                /* DOUBLE dmin_iz = (icell < icen) ? icell:icen; */
                /* dmin_iz *= (d2max-d2min)/ngrid; */
                /* dmin_iz += d2min; */
                /* DOUBLE dmin_iz = d2min + ((icen + icell)*cz_binsize)*0.5; */
                /* double theta = rpmax/(2.0*dmin_iz); */
                /* DOUBLE max_dec_sep=ASIN(theta)*2.0*INV_PI_OVER_180; */
                /* int dec_limits = (int) (ceil(max_dec_sep*inv_dec_diff*ngrid_dec[icell])); */
                /* fprintf(stderr,"icen = %d ngrid_dec[%d] = %d rbin_refine_factor=%d dec_limits=%d\n", */
                /*        icen,icell,ngrid_dec[icell],RBIN_REFINE_FACTOR,dec_limits); */
                int dec_limits = rbin_refine_factor;
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
                    int ra_limits = phibin_refine_factor;
                    for(int ira_step=-ra_limits;ira_step<=ra_limits;ira_step++) {
                        int ira = (ra_iz + ira_step + ngrid_ra[icell][idec]) % ngrid_ra[icell][idec];
                        //Linked in CZ, DEC and RA
                        const cellarray_mocks *cellstruct = &(lattice2[icell][idec][ira]);
#else
                        //Linked in CZ + DEC
                        const cellarray_mocks *cellstruct = &(lattice2[icell][idec]);
#endif

#else
                        //LINKED only in CZ
                        const cellarray_mocks *cellstruct = &(lattice2[icell]);
#endif
                        DOUBLE *x2  = cellstruct->pos;
                        DOUBLE *y2  = cellstruct->pos + NVEC;
                        DOUBLE *z2  = cellstruct->pos + 2*NVEC;
                        /*                  DOUBLE *cz2 = cellstruct->pos + 3*NVEC; */
                        /*                  const DOUBLE TWO=2.0; */
                        /*                  const DOUBLE sqr_d1 = d1[i]*d1[i]; */

#if !(defined(USE_AVX) && defined(__AVX__))

                        DOUBLE *localx2  = x2;
                        DOUBLE *localy2  = y2;
                        DOUBLE *localz2  = z2;
                        /*                  DOUBLE *localcz2 = cz2; */

                        for(int j=0;j<cellstruct->nelements;j+=NVEC) {
                            int block_size2=cellstruct->nelements - j;
                            if(block_size2 > NVEC) block_size2=NVEC;
                            for(int jj=0;jj<block_size2;jj++) {
                                const DOUBLE parx = x1 + localx2[jj];
                                const DOUBLE pary = y1 + localy2[jj];
                                const DOUBLE parz = z1 + localz2[jj];

                                const DOUBLE perpx = x1 - localx2[jj];
                                const DOUBLE perpy = y1 - localy2[jj];
                                const DOUBLE perpz = z1 - localz2[jj];

                                const DOUBLE sqr_s = perpx*perpx + perpy*perpy + perpz*perpz;
                                if(sqr_s >= sqr_max_sep) continue;

                                const DOUBLE tmp  = (parx*perpx+pary*perpy+parz*perpz);
                                const DOUBLE tmp1 = (parx*parx+pary*pary+parz*parz);
                                const DOUBLE sqr_Dpar = (tmp*tmp)/tmp1;
                                if(sqr_Dpar >= sqr_pimax) continue;

                                const int pibin  = (sqr_Dpar >= sqr_pimax) ? npibin:(int) (SQRT(sqr_Dpar)*inv_dpi);
                                const DOUBLE sqr_Dperp  = sqr_s - sqr_Dpar;
                                if(sqr_Dperp >= sqr_rpmax || sqr_Dperp < sqr_rpmin) continue;
#ifdef OUTPUT_RPAVG
                                const DOUBLE rp = SQRT(sqr_Dperp);
#endif

                                for(int kbin=nrpbin-1;kbin>=1;kbin--) {
                                    if(sqr_Dperp >= rupp_sqr[kbin-1]) {
                                        const int ibin = kbin*(npibin+1) + pibin;
                                        npairs[ibin]++;
#ifdef OUTPUT_RPAVG
                                        rpavg[ibin]+=rp;
#endif
                                        break;
                                    }
                                }
                            }//end of jj-loop

                            //increment localx2/localy2 etc
                            localx2   += 3*NVEC;//this might actually exceed the allocated range but we will never dereference that
                            localy2   += 3*NVEC;
                            localz2   += 3*NVEC;
                            /*                      localcz2  += 4*NVEC; */
                        }

#else //Use AVX intrinsics
                        AVX_FLOATS m_xpos    = AVX_SET_FLOAT(x1);
                        AVX_FLOATS m_ypos    = AVX_SET_FLOAT(y1);
                        AVX_FLOATS m_zpos    = AVX_SET_FLOAT(z1);
                        /*                  AVX_FLOATS m_sqr_d1  = AVX_SET_FLOAT(sqr_d1); */
                        union int8 {
                            AVX_INTS m_ibin;
                            int ibin[NVEC];
                        };


#ifdef OUTPUT_RPAVG
                        union float8{
                            AVX_FLOATS m_Dperp;
                            DOUBLE Dperp[NVEC];
                        };

#endif


                        DOUBLE *localx2  = x2;
                        DOUBLE *localy2  = y2;
                        DOUBLE *localz2  = z2;
                        /*                  DOUBLE *localcz2 = cz2; */
                        const AVX_FLOATS m_sqr_pimax  = AVX_SET_FLOAT(sqr_pimax);
                        const AVX_FLOATS m_sqr_rpmax  = AVX_SET_FLOAT(sqr_rpmax);
                        const AVX_FLOATS m_max_sep = AVX_SET_FLOAT(sqr_max_sep);
                        const AVX_FLOATS m_inv_dpi    = AVX_SET_FLOAT(inv_dpi);
                        const AVX_FLOATS m_sqr_rpmin  = AVX_SET_FLOAT(sqr_rpmin);
                        const AVX_FLOATS m_npibin     = AVX_SET_FLOAT((DOUBLE) npibin);
                        const AVX_FLOATS m_zero       = AVX_SET_FLOAT((DOUBLE) 0.0);
                        const AVX_FLOATS m_one    = AVX_SET_FLOAT((DOUBLE) 1);

                        int j;
                        for(j=0;j<=(cellstruct->nelements-NVEC);j+=NVEC){
                            const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(localx2);
                            const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(localy2);
                            const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(localz2);
                            /*                      const AVX_FLOATS m_cz2 = AVX_LOAD_FLOATS_UNALIGNED(localcz2); */

                            localx2  += 3*NVEC;
                            localy2  += 3*NVEC;
                            localz2  += 3*NVEC;
                            /*                      localcz2 += 4*NVEC; */

                            AVX_FLOATS m_sqr_Dpar, m_sqr_Dperp;
                            {
                                const AVX_FLOATS m_parx = AVX_ADD_FLOATS(m_x2, m_xpos);
                                const AVX_FLOATS m_pary = AVX_ADD_FLOATS(m_y2, m_ypos);
                                const AVX_FLOATS m_parz = AVX_ADD_FLOATS(m_z2, m_zpos);

                                const AVX_FLOATS m_perpx = AVX_SUBTRACT_FLOATS(m_xpos, m_x2);
                                const AVX_FLOATS m_perpy = AVX_SUBTRACT_FLOATS(m_ypos, m_y2);
                                const AVX_FLOATS m_perpz = AVX_SUBTRACT_FLOATS(m_zpos, m_z2);

                                const AVX_FLOATS m_term1 = AVX_MULTIPLY_FLOATS(m_parx, m_perpx);
                                const AVX_FLOATS m_term2 = AVX_MULTIPLY_FLOATS(m_pary, m_perpy);
                                const AVX_FLOATS m_term3 = AVX_MULTIPLY_FLOATS(m_parz, m_perpz);
                                const AVX_FLOATS m_numerator = AVX_SQUARE_FLOAT(AVX_ADD_FLOATS(m_term1, AVX_ADD_FLOATS(m_term2, m_term3)));

                                const AVX_FLOATS m_sqr_perpx = AVX_SQUARE_FLOAT(m_perpx);
                                const AVX_FLOATS m_sqr_perpy = AVX_SQUARE_FLOAT(m_perpy);
                                const AVX_FLOATS m_sqr_perpz = AVX_SQUARE_FLOAT(m_perpz);
                                const AVX_FLOATS m_sqr_sep = AVX_ADD_FLOATS(m_sqr_perpx, AVX_ADD_FLOATS(m_sqr_perpy, m_sqr_perpz));//3-d separation
                                //The 3-d separation (| s.s |)^2 *must* be less than (pimax^2 + rpmax^2). If not, one of the
                                //constraints for counting the pair (i.e., rp < rpmax, \pi < pimax) must be violated and
                                //we would discard the pair.
                                const AVX_FLOATS m_mask_3d_sep = AVX_COMPARE_FLOATS(m_sqr_sep, m_max_sep, _CMP_LT_OQ);

                                const AVX_FLOATS m_sqr_norm_l = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_parx), AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_pary), AVX_SQUARE_FLOAT(m_parz)));

                                //\pi ^2 = |s.l| ^2 / |l|^2
                                //However, division is slow -> so we will check if \pimax^2 * |l| ^2 < |s.l|^2. If not, then the
                                //value of \pi (after division) *must* be larger than \pimax -> in which case we would
                                //not count that pair anway.
                                const AVX_FLOATS m_sqr_pimax_times_l = AVX_MULTIPLY_FLOATS(m_sqr_pimax, m_sqr_norm_l);
                                const AVX_FLOATS m_mask_pimax_sep = AVX_COMPARE_FLOATS(m_numerator, m_sqr_pimax_times_l, _CMP_LT_OQ);// is pi < pimax ?
                                //If the bits are all 0, then *none* of the pairs satisfy the pimax + rpmax constraints.
                                const AVX_FLOATS m_mask = AVX_BITWISE_AND(m_mask_3d_sep, m_mask_pimax_sep);
                                if(AVX_TEST_COMPARISON(m_mask)==0) {
                                    continue;
                                }

#ifndef FAST_DIVIDE
                                m_sqr_Dpar = AVX_DIVIDE_FLOATS(m_numerator,m_sqr_norm_l);
                                //The divide is the actual operation we need
                                // but divides are about 10x slower than multiplies. So, I am replacing it
                                //with a approximate reciprocal in floating point
                                // + 2 iterations of newton-raphson in case of DOUBLE
#else //following blocks do an approximate reciprocal followed by two iterations of Newton-Raphson
#ifndef DOUBLE_PREC
                                //Taken from Intel's site: https://software.intel.com/en-us/articles/wiener-filtering-using-intel-advanced-vector-extensions
                                // (which has bugs in it, just FYI). Plus, https://techblog.lankes.org/2014/06/16/avx-isnt-always-faster-then-see/
                                __m256 rc  = _mm256_rcp_ps(m_sqr_norm_l);
                                __m256 two = AVX_SET_FLOAT(2.0f);
                                __m256 rc1 = _mm256_mul_ps(rc,
                                                           _mm256_sub_ps(two,
                                                                         _mm256_mul_ps(m_sqr_norm_l,rc)));

                                __m256 rc2 = _mm256_mul_ps(rc1,
                                                           _mm256_sub_ps(two,
                                                                         _mm256_mul_ps(m_sqr_norm_l,rc1)));

                                m_sqr_Dpar = _mm256_mul_ps ( m_numerator , rc2 );

#else
                                //we have to do this for doubles now.
                                //if the vrcpps instruction is not generated, there will
                                //be a ~70 cycle performance hit from switching between
                                //AVX and SSE modes.
                                __m128 float_tmp1 =  _mm256_cvtpd_ps(m_sqr_norm_l);
                                __m128 float_inv_tmp1 = _mm_rcp_ps(float_tmp1);
                                AVX_FLOATS rc = _mm256_cvtps_pd(float_inv_tmp1);

                                //We have the double->float->approx. reciprocal->double process done.
                                //Now improve the accuracy of the divide with newton-raphson.

                                //Ist iteration of NewtonRaphson
                                AVX_FLOATS two = AVX_SET_FLOAT(2.0);
                                AVX_FLOATS rc1 = _mm256_mul_pd(rc,
                                                               _mm256_sub_pd(two,
                                                                             _mm256_mul_pd(m_sqr_norm_l,rc)));
                                //2nd iteration of NewtonRaphson
                                AVX_FLOATS rc2 = _mm256_mul_pd(rc1,
                                                               _mm256_sub_pd(two,
                                                                             _mm256_mul_pd(m_sqr_norm_l,rc1)));
                                m_sqr_Dpar = AVX_MULTIPLY_FLOATS(m_numerator,rc2);
#endif//DOUBLE_PREC


#endif//FAST_DIVIDE

                                m_sqr_Dperp = AVX_SUBTRACT_FLOATS(m_sqr_sep,m_sqr_Dpar);
                            }


                            const AVX_FLOATS m_Dpar = AVX_SQRT_FLOAT(m_sqr_Dpar);

                            AVX_FLOATS m_mask_left;
                            //Do the mask filters in a separate scope
                            {
                                const AVX_FLOATS m_mask_pimax = AVX_COMPARE_FLOATS(m_sqr_Dpar,m_sqr_pimax,_CMP_LT_OQ);
                                const AVX_FLOATS m_rpmax_mask = AVX_COMPARE_FLOATS(m_sqr_Dperp, m_sqr_rpmax, _CMP_LT_OQ);
                                const AVX_FLOATS m_rpmin_mask = AVX_COMPARE_FLOATS(m_sqr_Dperp, m_sqr_rpmin, _CMP_GE_OQ);
                                const AVX_FLOATS m_rp_mask = AVX_BITWISE_AND(m_rpmax_mask,m_rpmin_mask);

                                m_mask_left = AVX_BITWISE_AND(m_mask_pimax, m_rp_mask);
                                if(AVX_TEST_COMPARISON(m_mask_left)==0) {
                                    continue;
                                }

                                m_sqr_Dperp = AVX_BLEND_FLOATS_WITH_MASK(m_zero,m_sqr_Dperp,m_mask_left);
                                m_sqr_Dpar  = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_pimax,m_sqr_Dpar,m_mask_left);
                            }
#ifdef OUTPUT_RPAVG
                            union float8 union_mDperp;
                            union_mDperp.m_Dperp = AVX_BLEND_FLOATS_WITH_MASK(m_zero,AVX_SQRT_FLOAT(m_sqr_Dperp),m_mask_left);
#endif
                            const AVX_FLOATS m_mask = m_mask_left;
                            AVX_FLOATS m_rpbin = AVX_SET_FLOAT((DOUBLE) 0);
                            for(int kbin=nrpbin-1;kbin>=1;kbin--) {
                                const AVX_FLOATS m_mask_low = AVX_COMPARE_FLOATS(m_sqr_Dperp,m_rupp_sqr[kbin-1],_CMP_GE_OQ);
                                const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m_mask_low,m_mask_left);
                                m_rpbin = AVX_BLEND_FLOATS_WITH_MASK(m_rpbin,m_kbin[kbin], m_bin_mask);
                                m_mask_left = AVX_COMPARE_FLOATS(m_sqr_Dperp, m_rupp_sqr[kbin-1],_CMP_LT_OQ);
                                if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
                                    break;
                                }
                            }

                            /* Compute the 1-D index to the [rpbin, pibin] := rpbin*(npibin+1) + pibin */
                            /*                      const AVX_FLOATS m_Dpar = AVX_SQRT_FLOAT(m_sqr_Dpar); */
                            const AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_Dpar,m_inv_dpi);
                            const AVX_FLOATS m_pibin = AVX_BLEND_FLOATS_WITH_MASK(m_npibin, m_tmp2, m_mask);
                            const AVX_FLOATS m_npibin_p1 = AVX_ADD_FLOATS(m_npibin,m_one);
                            const AVX_FLOATS m_binproduct = AVX_ADD_FLOATS(AVX_MULTIPLY_FLOATS(m_rpbin,m_npibin_p1),m_pibin);
                            union int8 union_finalbin;
                            union_finalbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_binproduct);

#if  __INTEL_COMPILER
#pragma unroll(NVEC)
#endif
                            for(int jj=0;jj<NVEC;jj++) {
                                const int ibin=union_finalbin.ibin[jj];

                                npairs[ibin]++;
#ifdef OUTPUT_RPAVG
                                rpavg [ibin] += union_mDperp.Dperp[jj];
#endif
                                /* fprintf(stderr,"i=%d j=%d union_rpbin.ibin[jj] = %d union_pibin.ibin[jj] = %d\n",i,j,union_rpbin.ibin[jj],union_pibin.ibin[jj]); */
                            }
                        }

                        //Take care of the remainder
                        for(int ipos=0;j<cellstruct->nelements;ipos++,j++) {
                            const DOUBLE parx = x1 + localx2[ipos];
                            const DOUBLE pary = y1 + localy2[ipos];
                            const DOUBLE parz = z1 + localz2[ipos];

                            const DOUBLE perpx = x1 - localx2[ipos];
                            const DOUBLE perpy = y1 - localy2[ipos];
                            const DOUBLE perpz = z1 - localz2[ipos];

                            const DOUBLE sqr_s = perpx*perpx + perpy*perpy + perpz*perpz;
                            if(sqr_s >= sqr_max_sep) continue;

                            const DOUBLE tmp  = (parx*perpx+pary*perpy+parz*perpz);
                            const DOUBLE tmp1 = (parx*parx+pary*pary+parz*parz);
                            const DOUBLE sqr_Dpar = (tmp*tmp)/tmp1;
                            if(sqr_Dpar >= sqr_pimax) continue;

                            const int pibin  = (sqr_Dpar >= sqr_pimax) ? npibin:(int) (SQRT(sqr_Dpar)*inv_dpi);
                            const DOUBLE sqr_Dperp  = sqr_s - sqr_Dpar;
                            if(sqr_Dperp >= sqr_rpmax || sqr_Dperp < sqr_rpmin) continue;
#ifdef OUTPUT_RPAVG
                            const DOUBLE rp = SQRT(sqr_Dperp);
#endif
                            for(int kbin=nrpbin-1;kbin>=1;kbin--) {
                                if(sqr_Dperp >= rupp_sqr[kbin-1]) {
                                    const int ibin = kbin*(npibin+1) + pibin;
                                    npairs[ibin]++;
#ifdef OUTPUT_RPAVG
                                    rpavg[ibin]+=rp;
#endif
                                    break;
                                }
                            }
                        }
#endif  //END of the AVX/NO-AVX section

#ifdef LINK_IN_DEC
#ifdef LINK_IN_RA
                    }//loop over ra bins
#endif
                }//loop over dec bins
#endif
            }//icell loop over cz cells
        }//i loop over ND1 particles
#if defined(USE_OMP) && defined(_OPENMP)
        for(int i=0;i<totnbins;i++) {
            all_npairs[tid][i] = npairs[i];
#ifdef OUTPUT_RPAVG
            all_rpavg[tid][i] = rpavg[i];
#endif
        }

    }//close the omp parallel region
#endif//USE_OMP

#ifndef SILENT    
    finish_myprogressbar(&interrupted);
#endif    

#if defined(USE_OMP) && defined(_OPENMP)
    uint64_t npairs[totnbins];
#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[totnbins];
#endif
    for(int i=0;i<totnbins;i++) {
        npairs[i] = 0;
#ifdef OUTPUT_RPAVG
        rpavg[i] = 0.0;
#endif
    }

    for(int i=0;i<numthreads;i++) {
        for(int j=0;j<totnbins;j++) {
            npairs[j] += all_npairs[i][j];
#ifdef OUTPUT_RPAVG
            rpavg[j] += all_rpavg[i][j];
#endif
        }
    }
    matrix_free((void **) all_npairs, numthreads);
#ifdef OUTPUT_RPAVG
    matrix_free((void **) all_rpavg, numthreads);
#endif
#endif //USE_OMP



#ifdef OUTPUT_RPAVG
    for(int i=0;i<totnbins;i++){
        if(npairs[i] > 0) {
            rpavg[i] /= ((DOUBLE) npairs[i] );
        }
    }
#endif


#ifndef LINK_IN_DEC
    for(int i=0;i < ngrid;i++) {
        free(lattice2[i].pos);
        /* free(lattice2[i].y); */
        /* free(lattice2[i].z); */
        /* free(lattice2[i].cz); */
    }
    free(lattice2);
#else
#ifndef LINK_IN_RA
    for(int i=0; i < ngrid; i++) {
        for(int j=0;j<ngrid_dec[i];j++) {
            free(lattice2[i][j].pos);
            /* free(lattice2[i][j].y); */
            /* free(lattice2[i][j].z); */
            /* free(lattice2[i][j].cz); */
        }
        free(lattice2[i]);
    }
    free(lattice2);
#else
    //LINK_IN_RA
    int max_nmesh_dec=0;
    for(int i=0; i < ngrid; i++) {
        for(int j=0;j< ngrid_dec[i];j++) {
            if(ngrid_dec[i] > max_nmesh_dec) max_nmesh_dec = ngrid_dec[i];
            for(int k=0;k<ngrid_ra[i][j];k++){
                free(lattice2[i][j][k].pos);
            }
        }
    }

    volume_free((void ***) lattice2, ngrid, max_nmesh_dec);
    matrix_free((void **) ngrid_ra, ngrid);
    //LINK_IN_RA
#endif
    free(ngrid_dec);
#endif

    free(d1);
    if(autocorr==0) free(d2);

    //Pack in the results
    results_countpairs_mocks *results = my_malloc(sizeof(*results), 1);
    results->nbin   = nrpbin;
    results->npibin = npibin;
    results->pimax  = pimax;
    results->npairs = my_malloc(sizeof(uint64_t), totnbins);
    results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
    results->rpavg  = my_malloc(sizeof(DOUBLE)  , totnbins);

    for(int i=0;i<nrpbin;i++) {
        results->rupp[i] = rupp[i];
        for(int j=0;j<npibin;j++) {
            int index = i*(npibin+1) + j;
            assert(index < totnbins && "index must be in range");
            results->npairs[index] = npairs[index];
#ifdef OUTPUT_RPAVG
            results->rpavg[index] = rpavg[index];
#else
            results->rpavg[index] = 0.0;
#endif
        }
    }
    free(rupp);

    return results;
}
