/* File: test_periodic.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include "tests_common.h"
#include "io.h"

#include "../DD/countpairs.h"
#include "../DDrppi/countpairs_rp_pi.h"
#include "../DDsmu/countpairs_s_mu.h"
#include "../wp/countpairs_wp.h"
#include "../xi/countpairs_xi.h"
#include "../vpf/countspheres.h"

char tmpoutputfile[]="./test_periodic_output.txt";

int test_periodic_DD(const char *correct_outputfile);
int test_periodic_DDrppi(const char *correct_outputfile);
int test_periodic_DDsmu(const char *correct_outputfile);
int test_wp(const char *correct_outputfile);
int test_vpf(const char *correct_outputfile);
int test_xi(const char *correct_outputfile);

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,
                               const char *secondfilename, const char *secondformat);


//Global variables
int64_t ND1;
double *X1=NULL,*Y1=NULL,*Z1=NULL,*weights1=NULL;

int64_t ND2;
double *X2=NULL,*Y2=NULL,*Z2=NULL,*weights2=NULL;

char current_file1[MAXLEN+1],current_file2[MAXLEN+1];

struct config_options options;
//end global variables

int test_periodic_DD(const char *correct_outputfile)
{
    int autocorr = (X1==X2) ? 1:0;
    results_countpairs results;
    int ret = EXIT_FAILURE;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    BEGIN_INTEGRATION_TEST_SECTION

        //Do the straight-up DD counts
        int status = countpairs(ND1,X1,Y1,Z1,
                                ND2,X2,Y2,Z2,
                                nthreads,
                                autocorr,
                                binfile,
                                &results,
                                &options,
                                &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        FILE *fp = my_fopen(correct_outputfile,"r");
        for(int i=1;i<results.nbin;i++) {
            uint64_t npairs;
            double rpavg, weightavg;
            ret = EXIT_FAILURE;
            int nitems = fscanf(fp,"%"SCNu64" %lf %*f %*f %lf", &npairs, &rpavg, &weightavg);
            if(nitems != 3) {
                ret = EXIT_FAILURE;//not required but showing intent
                break;
            }
            int floats_equal = AlmostEqualRelativeAndAbs_double(rpavg, results.rpavg[i], maxdiff, maxreldiff);
            int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[i], maxdiff, maxreldiff);

            //Check for exact equality of npairs and float "equality" for rpavg
            if(npairs == results.npairs[i] && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
                ret = EXIT_SUCCESS;
            } else {
                ret = EXIT_FAILURE;//not required but showing intent
                fprintf(stderr,"Failed. True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[i]);
                fprintf(stderr,"Failed. True rpavg = %e Computed rpavg = %e. floats_equal = %d\n", rpavg, results.rpavg[i], floats_equal);
                fprintf(stderr,"Failed. True weightavg = %e Computed weightavg = %e. weights_equal = %d\n", weightavg, results.weightavg[i], weights_equal);
                break;
            }

        }
        fclose(fp);

    END_INTEGRATION_TEST_SECTION(free_results(&results));

    if(ret != EXIT_SUCCESS && results.nbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results(&results);
            return EXIT_FAILURE;
        }
        double rlow=results.rupp[0];
        for(int i=1;i<results.nbin;i++) {
            fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf \n",results.npairs[i],results.rpavg[i],rlow,results.rupp[i],results.weightavg[i]);
            rlow = results.rupp[i];
        }
        fclose(fp);
    }

    free_results(&results);
    return ret;
}


int test_periodic_DDrppi(const char *correct_outputfile)
{
    results_countpairs_rp_pi results;
    int ret = EXIT_FAILURE;
    int autocorr = (X1==X2) ? 1:0;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    BEGIN_INTEGRATION_TEST_SECTION
        int status = countpairs_rp_pi(ND1,X1,Y1,Z1,
                                      ND2,X2,Y2,Z2,
                                      nthreads,
                                      autocorr,
                                      binfile,
                                      pimax,
                                      &results,
                                      &options,
                                      &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        const int npibin = results.npibin;
        FILE *fp = my_fopen(correct_outputfile, "r");
        for(int i=1;i<results.nbin;i++) {
            for(int j=0;j<npibin;j++) {
                int index = i*(npibin+1) + j;
                uint64_t npairs;
                double rpavg, weightavg;
                ret = EXIT_FAILURE;
                int nitems = fscanf(fp,"%"SCNu64" %lf %*f %*f %lf%*[^\n]", &npairs, &rpavg, &weightavg);
                if(nitems != 3) {
                    ret = EXIT_FAILURE;//not required but showing intent
                    i = results.nbin;
                    break;
                }
                int floats_equal = AlmostEqualRelativeAndAbs_double(rpavg, results.rpavg[index], maxdiff, maxreldiff);
                int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[index], maxdiff, maxreldiff);

                //Check for exact equality of npairs and float "equality" for rpavg
                if(npairs == results.npairs[index] && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
                    ret = EXIT_SUCCESS;
                } else {
                    fprintf(stderr,"True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[index]);
                    fprintf(stderr,"True rpavg  = %20.12e Computed rpavg = %20.12e. floats_equal = %d\n", rpavg, results.rpavg[index], floats_equal);
                    fprintf(stderr,"True weightavg = %e Computed weightavg = %e. weights_equal = %d\n", weightavg, results.weightavg[index], weights_equal);

                    ret = EXIT_FAILURE;//not required but showing intent
                    i = results.nbin;
                    break;
                }
            }
        }
        fclose(fp);
    END_INTEGRATION_TEST_SECTION(free_results_rp_pi(&results));

    if(ret != EXIT_SUCCESS && results.nbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results_rp_pi(&results);
            return EXIT_FAILURE;
        }

        const int npibin = results.npibin;
        const double dpi = pimax/(double)results.npibin ;
        for(int i=1;i<results.nbin;i++) {
            const double logrp = log10(results.rupp[i]);
            for(int j=0;j<npibin;j++) {
                int index = i*(npibin+1) + j;
                fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf\n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi, results.weightavg[index]);
            }
        }
        fclose(fp);
    }

    //free the result structure
    free_results_rp_pi(&results);
    return ret;
}

int test_periodic_DDsmu(const char *correct_outputfile)
{
    int autocorr = (X1==X2) ? 1:0;
    int ret = EXIT_FAILURE;
    results_countpairs_s_mu results;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    BEGIN_INTEGRATION_TEST_SECTION
        int status = countpairs_s_mu(ND1,X1,Y1,Z1,
                                     ND2,X2,Y2,Z2,
                                     nthreads,
                                     autocorr,
                                     binfile,
                                     theory_mu_max,
                                     nmu_bins,
                                     &results,
                                     &options,
                                     &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        ret = EXIT_FAILURE;
        const int nmubin = results.nmu_bins;
        FILE *fp = my_fopen(correct_outputfile, "r");
        for(int i=1;i<results.nsbin;i++) {
            for(int j=0;j<nmubin;j++) {
                int index = i*(nmubin+1) + j;
                uint64_t npairs;
                double savg, weightavg;
                ret = EXIT_FAILURE;
                int nitems = fscanf(fp,"%"SCNu64" %lf %*f %*f %lf%*[^\n]", &npairs, &savg, &weightavg);
                if(nitems != 3) {
                    ret = EXIT_FAILURE;//not required but showing intent
                    i = results.nsbin;
                    break;
                }
                int floats_equal = AlmostEqualRelativeAndAbs_double(savg, results.savg[index], maxdiff, maxreldiff);
                int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[index], maxdiff, maxreldiff);

                //Check for exact equality of npairs and float "equality" for savg
                if(npairs == results.npairs[index] && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
                    ret = EXIT_SUCCESS;
                } else {
                    fprintf(stderr,"True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[index]);
                    fprintf(stderr,"True savg  = %20.12e Computed savg = %20.12e. floats_equal = %d\n", savg, results.savg[index], floats_equal);
                    fprintf(stderr,"True weightavg = %e Computed weightavg = %e. weights_equal = %d\n", weightavg, results.weightavg[index], weights_equal);

                    ret = EXIT_FAILURE;//not required but showing intent
                    i = results.nsbin;
                    break;
                }
            }
        }
        fclose(fp);
    END_INTEGRATION_TEST_SECTION(free_results_s_mu(&results));

    if(ret != EXIT_SUCCESS && results.nsbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results_s_mu(&results);
            return EXIT_FAILURE;
        }
        const int nmubin = results.nmu_bins;
        const double dmu = theory_mu_max/(double)results.nmu_bins;
        for(int i=1;i<results.nsbin;i++) {
            const double logs = log10(results.supp[i]);
            for(int j=0;j<nmubin;j++) {
                int index = i*(nmubin+1) + j;
                fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf\n",results.npairs[index],results.savg[index],logs,(j+1)*dmu, results.weightavg[index]);
            }
        }
        fclose(fp);
    }

    //free the result structure
    free_results_s_mu(&results);
    return ret;
}


int test_wp(const char *correct_outputfile)
{
    results_countpairs_wp results;
    int ret = EXIT_FAILURE;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;

    BEGIN_INTEGRATION_TEST_SECTION
        int status = countpairs_wp(ND1,X1,Y1,Z1,
                                   boxsize,
                                   nthreads,
                                   binfile,
                                   pimax,
                                   &results,
                                   &options,
                                   &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        FILE *fp=my_fopen(correct_outputfile,"r");
        for(int i=1;i<results.nbin;i++) {
            uint64_t npairs;
            double rpavg,wp,weightavg;
            ret = EXIT_FAILURE;
            int nitems = fscanf(fp,"%lf %lf %*f %*f %"SCNu64" %lf%*[^\n]", &wp, &rpavg, &npairs, &weightavg);//discard rlow and rupp
            if(nitems != 4) {
                ret = EXIT_FAILURE;//not required but showing intent
                break;
            }
            int rpavg_equal = AlmostEqualRelativeAndAbs_double(rpavg, results.rpavg[i], maxdiff, maxreldiff);
            int weightavg_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[i], maxdiff, maxreldiff);
            int wp_equal = AlmostEqualRelativeAndAbs_double(wp, results.wp[i], maxdiff, maxreldiff);

            //Check for exact equality of npairs and float "equality" for rpavg + wp
            if(npairs == results.npairs[i] && rpavg_equal == EXIT_SUCCESS && wp_equal == EXIT_SUCCESS && weightavg_equal == EXIT_SUCCESS) {
                ret = EXIT_SUCCESS;
            } else {
                ret = EXIT_FAILURE;//not required but showing intent
                fprintf(stderr,"Failed. True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[i]);
                fprintf(stderr,"Failed. True wp = %e Computed results wp = %e\n", wp, results.wp[i]);
                fprintf(stderr,"Failed. True rpavg = %e Computed rpavg = %e.\n",
                        rpavg, results.rpavg[i]);
                fprintf(stderr,"Failed. True weightavg = %e Computed weightavg = %e.\n",
                        weightavg, results.weightavg[i]);
                fprintf(stderr," wp_equal = %d rpavg_equal = %d weightavg_equal = %d\n", wp_equal, rpavg_equal, weightavg_equal);
                break;
            }
        }
        fclose(fp);
    END_INTEGRATION_TEST_SECTION(free_results_wp(&results));

    /* Test failed. Output the results into a temporary file */
    if(ret != EXIT_SUCCESS && results.nbin > 0 ) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results_wp(&results);
            return EXIT_FAILURE;
        }
        double rlow=results.rupp[0];
        for(int i=1;i<results.nbin;++i) {
            fprintf(fp,"%e\t%e\t%e\t%e\t%12"PRIu64"\t%e \n",results.wp[i],results.rpavg[i],rlow,results.rupp[i],results.npairs[i], results.weightavg[i]);
            rlow=results.rupp[i];
        }
        fclose(fp);
    }

    //free the result structure
    free_results_wp(&results);
    return ret;
}

int test_vpf(const char *correct_outputfile)
{
    const double rmax = 10.0;
    const int nbin = 10;
    const int nc = 10000;
    const int num_pN=6;
    const unsigned long seed=-1234;
    results_countspheres results;
    int ret = EXIT_FAILURE;

    BEGIN_INTEGRATION_TEST_SECTION
        int status = countspheres(ND1, X1, Y1, Z1,
                                  rmax, nbin, nc,
                                  num_pN,
                                  seed,
                                  &results,
                                  &options, NULL);

        if(status != EXIT_SUCCESS) {
            return status;
        }

        //Output the results
        FILE *fp=my_fopen(correct_outputfile,"r");
        if(fp == NULL) {
            free_results_countspheres(&results);
            return EXIT_FAILURE;
        }
        for(int ibin=0;ibin<results.nbin;ibin++) {
            double r;
            int nitems = fscanf(fp, "%lf", &r);
            if(nitems != 1) {
                return EXIT_FAILURE;
            }
            ret = EXIT_FAILURE;
            for(int i=0;i<num_pN;i++) {
                double pN;
                nitems = fscanf(fp, " %lf ", &pN);
                if(nitems != 1) {
                    return EXIT_FAILURE;
                }

                /* Not quite sure how this is working. The correct output columns only have 4 digits printed,
                   but I am comparing here with ~1e-9 in abs. diff. The only way the comparison should work is
                   if the conversion to 4 digits during printf, round-trips during scanf. But surely there must
                   be a lot more doubles that can be fit within those missing digits of precision.

                   I would have thought the comparison would require maxdiff ~ 1e-4. -- MS
                */
                int floats_equal = AlmostEqualRelativeAndAbs_double(pN, (results.pN)[ibin][i], maxdiff, maxreldiff);
                if(floats_equal != EXIT_SUCCESS) {
                    ibin=results.nbin;
                    ret=EXIT_FAILURE;
                    break;
                }
                ret = EXIT_SUCCESS;
            }
        }
        fclose(fp);
    END_INTEGRATION_TEST_SECTION(free_results_countspheres(&results));

    if(ret != EXIT_SUCCESS && results.nbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        const double rstep = rmax/(double)nbin ;
        for(int ibin=0;ibin<results.nbin;ibin++) {
            const double r=(ibin+1)*rstep;
            fprintf(fp,"%lf ", r);
            for(int i=0;i<num_pN;i++) {
                fprintf(fp," %10.4e", (results.pN)[ibin][i]);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }

    //free the result structure
    free_results_countspheres(&results);
    return ret;
}

int test_xi(const char *correct_outputfile)
{
    results_countpairs_xi results;
    int ret = EXIT_FAILURE;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;

    BEGIN_INTEGRATION_TEST_SECTION
        int status = countpairs_xi(ND1,X1,Y1,Z1,
                                   boxsize,
                                   nthreads,
                                   binfile,
                                   &results,
                                   &options,
                                   &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        FILE *fp=my_fopen(correct_outputfile,"r");
        for(int i=1;i<results.nbin;i++) {
            uint64_t npairs;
            double ravg,xi,weightavg;
            ret = EXIT_FAILURE;
            int nitems = fscanf(fp,"%lf %lf %*f %*f %"SCNu64" %lf[^\n]", &xi, &ravg, &npairs, &weightavg); //discard rlow and rupp
            if(nitems != 4) {
                ret = EXIT_FAILURE;//not required but showing intent
                break;
            }
            int ravg_equal = AlmostEqualRelativeAndAbs_double(ravg, results.ravg[i], maxdiff, maxreldiff);
            int weightavg_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[i], maxdiff, maxreldiff);
            int xi_equal = AlmostEqualRelativeAndAbs_double(xi, results.xi[i], maxdiff, maxreldiff);

            //Check for exact equality of npairs and float "equality" for ravg + xi
            if(npairs == results.npairs[i] && ravg_equal == EXIT_SUCCESS && xi_equal == EXIT_SUCCESS && weightavg_equal == EXIT_SUCCESS) {
                ret = EXIT_SUCCESS;
            } else {
                ret = EXIT_FAILURE;//not required but showing intent
                fprintf(stderr,"Failed. True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[i]);
                fprintf(stderr,"Failed. True xi = %e Computed results xi = %e\n", xi, results.xi[i]);
                fprintf(stderr,"Failed. True ravg = %e Computed ravg = %e.\n",
                        ravg, results.ravg[i]);
                fprintf(stderr,"Failed. True weightavg = %e Computed weightavg = %e.\n",
                        weightavg, results.weightavg[i]);
                fprintf(stderr," xi_equal = %d ravg_equal = %d weightavg_equal = %d\n", xi_equal, ravg_equal, weightavg_equal);
                break;
            }
        }
        fclose(fp);
    END_INTEGRATION_TEST_SECTION(free_results_xi(&results));

    if(ret != EXIT_SUCCESS && results.nbin > 0){
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results_xi(&results);
            return EXIT_FAILURE;
        }
        double rlow=results.rupp[0];
        for(int i=1;i<results.nbin;++i) {
            fprintf(fp,"%e\t%e\t%e\t%e\t%12"PRIu64"\t%e \n",results.xi[i],results.ravg[i],rlow,results.rupp[i],results.npairs[i], results.weightavg[i]);
            rlow=results.rupp[i];
        }
        fclose(fp);
    }

    //free the result structure
    free_results_xi(&results);
    return ret;
}

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,
                               const char *secondfilename, const char *secondformat)
{
    int free_X2=0;
    if(X2 != NULL && X2 != X1) {
        free_X2=1;
    }


    //Check to see if data has to be read for X1/Y1/Z1
    if (strncmp(current_file1,firstfilename,strlen(current_file1)) != 0) {
        //replace the data-set
        if(X1 != NULL) {
            free(X1);free(Y1);free(Z1);free(weights1);
        }

        //Since X2 was pointing to X1, need to signal that the memory location is no longer valid
        if(free_X2 == 0) {
            X2 = NULL;
            Y2 = NULL;
            Z2 = NULL;
            weights2 = NULL;
        }
        ND1 = read_positions(firstfilename,firstformat, sizeof(double), 4, &X1, &Y1, &Z1, &weights1);
        strncpy(current_file1,firstfilename,MAXLEN);
    }

    //first check if only one unique file is asked for
    if(strncmp(firstfilename,secondfilename,strlen(firstfilename))==0) {
        //But X2 might have read-in a different file->avoid a memory-leak
        if(free_X2 == 1) {
            free(X2);free(Y2);free(Z2);free(weights2);
            free_X2 = 0;//not essential since the code returns after this section
        }
        X2=X1;
        Y2=Y1;
        Z2=Z1;
        weights2=weights1;
        ND2=ND1;
        strncpy(current_file2,secondfilename,MAXLEN);
        return;
    }


    //Check to see if data has to be read for X2/Y2/Z2
    if (strncmp(current_file2,secondfilename,strlen(current_file2)) != 0 || X2 == NULL) {
        //replace the data-set
        if(free_X2 == 1) {
            free(X2);free(Y2);free(Z2);free(weights2);
        }
        ND2 = read_positions(secondfilename,secondformat, sizeof(double), 4, &X2, &Y2, &Z2, &weights2);
        strncpy(current_file2,secondfilename,MAXLEN);
    }
}

int main(int argc, char **argv)
{
    struct timeval tstart,t0,t1;
    options = get_config_options();
    options.need_avg_sep=1;
    options.verbose=0;
    options.periodic=1;
    options.fast_divide_and_NR_steps=0;
    options.copy_particles=1;
    options.float_type=sizeof(double);
    options.boxsize = boxsize;

    char file[]="../tests/data/gals_Mr19.ff";
    char fileformat[]="f";

    gettimeofday(&tstart,NULL);

    //set the globals
    ND1 = read_positions(file,fileformat, sizeof(double), 4, &X1, &Y1, &Z1, &weights1);

    ND2 = ND1;
    X2 = X1;
    Y2 = Y1;
    Z2 = Z1;
    weights2 = weights1;

    strncpy(current_file1,file,MAXLEN);
    strncpy(current_file2,file,MAXLEN);
    reset_bin_refine_factors(&options);

    int failed=0;
    int status;

    const char alltests_names[][MAXLEN] = {"Mr19 DDrppi (periodic)",
                                           "Mr19 DD (periodic)",
                                           "Mr19 wp (periodic)",
                                           "Mr19 vpf (periodic)",
                                           "Mr19 xi (periodic)",
                                           "Mr19 DDsmu (periodic)",
                                           "CMASS DDrppi DD (periodic)",
                                           "CMASS DDrppi DR (periodic)",
                                           "CMASS DDrppi RR (periodic)"};
    const int ntests = sizeof(alltests_names)/(sizeof(char)*MAXLEN);
    const int function_pointer_index[] = {1,0,2,3,4,5,1,1,1};//0->DD, 1->DDrppi,2->wp, 3->vpf, 4->xi, 5->DDsmu

    const char correct_outputfiles[][MAXLEN] = {"Mr19_DDrppi_periodic",
                                                "Mr19_DD_periodic",
                                                "Mr19_wp",
                                                "Mr19_vpf_periodic",
                                                "Mr19_xi",
                                                "Mr19_DDsmu_periodic",
                                                "cmass_DD_periodic",
                                                "cmass_DR_periodic",
                                                "cmass_RR_periodic"};
    const char firstfilename[][MAXLEN] = {"../tests/data/gals_Mr19.ff",
                                          "../tests/data/gals_Mr19.ff",
                                          "../tests/data/gals_Mr19.ff",
                                          "../tests/data/gals_Mr19.ff",
                                          "../tests/data/gals_Mr19.ff",
                                          "../tests/data/gals_Mr19.ff",
                                          "../tests/data/cmassmock_Zspace.ff",
                                          "../tests/data/cmassmock_Zspace.ff",
                                          "../tests/data/random_Zspace.ff"};
    const char firstfiletype[][MAXLEN] = {"f","f","f","f","f","f","f","f","f"};
    const char secondfilename[][MAXLEN] = {"../tests/data/gals_Mr19.ff",
                                           "../tests/data/gals_Mr19.ff",
                                           "../tests/data/gals_Mr19.ff",
                                           "../tests/data/gals_Mr19.ff",
                                           "../tests/data/gals_Mr19.ff",
                                           "../tests/data/gals_Mr19.ff",
                                           "../tests/data/cmassmock_Zspace.ff",
                                           "../tests/data/random_Zspace.ff",
                                           "../tests/data/random_Zspace.ff"};
    const char secondfiletype[][MAXLEN] = {"f","f","f","f","f","f","f","f","f"};
    const double allpimax[]             = {40.0,40.0,40.0,40.0,40.0,40.0,80.0,80.0,80.0};
    const double allboxsize[]             = {boxsize,boxsize,boxsize,boxsize,boxsize,boxsize,0.,0.,0.};

    int (*allfunctions[]) (const char *) = {test_periodic_DD,
                                            test_periodic_DDrppi,
                                            test_wp,
                                            test_vpf,
                                            test_xi,
                                            test_periodic_DDsmu};
    const int numfunctions=6;//6 functions total

    int total_tests=0,skipped=0;

    if(argc==1) {
        //nothing was passed at the command-line -> run all tests
        for(int i=0;i<ntests;i++) {
            int function_index = function_pointer_index[i];
            assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
            const char *testname = alltests_names[i];
            int skip_test=test_all_files_present(2,firstfilename[i],secondfilename[i]);
            if(skip_test != 0) {
                fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". File(s) not found\n", testname);
                skipped++;
                continue;
            }
            read_data_and_set_globals(firstfilename[i],firstfiletype[i],secondfilename[i],secondfiletype[i]);
            pimax=allpimax[i];
            options.boxsize = allboxsize[i];
            gettimeofday(&t0,NULL);
            status = (*allfunctions[function_index])(correct_outputfiles[i]);
            gettimeofday(&t1,NULL);
            double pair_time = ADD_DIFF_TIME(t0,t1);
            total_tests++;
            if(status==EXIT_SUCCESS) {
                fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                run_system_call(execstring);//ignoring status
            } else {
                fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                failed++;
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,i);
                run_system_call(execstring);//ignoring status
                fprintf(stderr, ANSI_COLOR_RED "Failed output copied to %s.%d correct output is in %s"ANSI_COLOR_RESET"\n", tmpoutputfile, i, correct_outputfiles[i]);
            }
        }
    } else {
        //run specific tests
        for(int i=1;i<argc;i++){
            int this_test_num = atoi(argv[i]);
            if(this_test_num >= 0 && this_test_num < ntests) {
                int function_index = function_pointer_index[this_test_num];
                assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
                const char *testname = alltests_names[this_test_num];
                int skip_test=test_all_files_present(2,firstfilename[this_test_num],secondfilename[this_test_num]);
                if(skip_test != 0) {
                    fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". File(s) not found\n", testname);
                    skipped++;
                    continue;
                }
                total_tests++;
                read_data_and_set_globals(firstfilename[this_test_num],firstfiletype[this_test_num],secondfilename[this_test_num],secondfiletype[this_test_num]);
                pimax=allpimax[this_test_num];
                options.boxsize = allboxsize[this_test_num];
                gettimeofday(&t0,NULL);
                status = (*allfunctions[function_index])(correct_outputfiles[this_test_num]);
                gettimeofday(&t1,NULL);
                double pair_time = ADD_DIFF_TIME(t0,t1);
                if(status==EXIT_SUCCESS) {
                    fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                    char execstring[MAXLEN];
                    my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                    run_system_call(execstring);//ignoring status
                } else {
                    fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                    failed++;
                    char execstring[MAXLEN];
                    my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,this_test_num);
                    run_system_call(execstring);//ignoring status
                    fprintf(stderr, ANSI_COLOR_RED "Failed output copied to %s.%d correct output is in %s"ANSI_COLOR_RESET"\n", tmpoutputfile, this_test_num, correct_outputfiles[this_test_num]);
                }
            }
        }
    }


    gettimeofday(&t1,NULL);
    double total_time = ADD_DIFF_TIME(tstart,t1);
    if(failed > 0) {
        fprintf(stderr,ANSI_COLOR_RED "FAILED %d out of %d tests. Total time = %8.2lf seconds " ANSI_COLOR_RESET "\n", failed, total_tests, total_time);
    } else {
        fprintf(stderr,ANSI_COLOR_GREEN "PASSED: ALL %d tests. Total time = %8.2lf seconds " ANSI_COLOR_RESET "\n", total_tests, total_time);
    }
    if(skipped > 0) {
        fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: %d tests" ANSI_COLOR_RESET "\n", skipped);
        fprintf(stderr,ANSI_COLOR_MAGENTA "Tests are skipped on the PyPI installed code-base. Please use the git repo if you want to run the entire suite of tests"ANSI_COLOR_RESET"\n\n");
    }

    if(X2 != X1) {
        free(X2);free(Y2);free(Z2);free(weights2);
    }
    free(X1);free(Y1);free(Z1);free(weights1);
    return failed;
}
