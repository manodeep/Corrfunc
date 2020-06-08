/* File: tests_mocks.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include "tests_common.h"
#include "io.h"
#include "utils.h"
#include "cosmology_params.h"

#if !(defined(__INTEL_COMPILER)) && defined(USE_AVX)
#warning Test suite for mocks are faster with Intel compiler, icc, AVX libraries.
#endif


#include "../DDrppi_mocks/countpairs_rp_pi_mocks.h"
#include "../DDsmu_mocks/countpairs_s_mu_mocks.h"
#include "../DDtheta_mocks/countpairs_theta_mocks.h"
#include "../vpf_mocks/countspheres_mocks.h"

char tmpoutputfile[]="../tests/tests_mocks_output.txt";

int test_DDrppi_mocks(const char *correct_outputfile);
int test_DDtheta_mocks(const char *correct_outputfile);
int test_vpf_mocks(const char *correct_outputfile);
int test_DDsmu_mocks(const char *correct_outputfile);

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat);

//Global variables
int64_t ND1;
double *RA1=NULL,*DEC1=NULL,*CZ1=NULL,*weights1=NULL;

int64_t ND2;
double *RA2=NULL,*DEC2=NULL,*CZ2=NULL,*weights2=NULL;

const int cosmology_flag=1;
char current_file1[MAXLEN+1],current_file2[MAXLEN+1];

struct config_options options;
//end of global variables

int test_DDrppi_mocks(const char *correct_outputfile)
{
    results_countpairs_mocks results;
    int ret = EXIT_FAILURE;
    assert(RA1 != NULL && DEC1 != NULL && CZ1 != NULL && "ERROR: In test suite for DDrppi ra/dec/cz can not be NULL pointers");
    int autocorr = (RA1==RA2) ? 1:0;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    //Do DD(rp,pi) counts
    BEGIN_INTEGRATION_TEST_SECTION
        int status = countpairs_mocks(ND1,RA1,DEC1,CZ1,
                                      ND2,RA2,DEC2,CZ2,
                                      nthreads,
                                      autocorr,
                                      binfile,
                                      pimax,
                                      cosmology_flag,
                                      &results,
                                      &options,
                                      &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        FILE *fp=my_fopen(correct_outputfile,"r");
        if(fp == NULL) {
            free_results_mocks(&results);
            return EXIT_FAILURE;
        }
        const int npibin = results.npibin;
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
    END_INTEGRATION_TEST_SECTION(free_results_mocks(&results));

    /* If the test failed, then write the output into a temporary file */
    if(ret != EXIT_SUCCESS && results.nbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results_mocks(&results);
            return EXIT_FAILURE;
        }
        const double dpi = pimax/(double)results.npibin ;
        const int npibin = results.npibin;
        for(int i=1;i<results.nbin;i++) {
            const double logrp = log10(results.rupp[i]);
            for(int j=0;j<npibin;j++) {
                int index = i*(npibin+1) + j;
                fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf %20.8lf \n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi, results.weightavg[index]);
            }
        }
        fclose(fp);
    }

    free_results_mocks(&results);
    return ret;
}

int test_DDsmu_mocks(const char *correct_outputfile)
{
    results_countpairs_mocks_s_mu results;
    int ret = EXIT_FAILURE;

    assert(RA1 != NULL && DEC1 != NULL && CZ1 != NULL && "ERROR: In test suite for DDsmu ra/dec/cz can not be NULL pointers");
    int autocorr = (RA1==RA2) ? 1:0;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    BEGIN_INTEGRATION_TEST_SECTION
        //Do DD(s,mu) counts
        int status = countpairs_mocks_s_mu(ND1,RA1,DEC1,CZ1,
                                           ND2,RA2,DEC2,CZ2,
                                           nthreads,
                                           autocorr,
                                           binfile,
                                           mocks_mu_max,
                                           nmu_bins,
                                           cosmology_flag,
                                           &results,
                                           &options,
                                           &extra);
        if(status != EXIT_SUCCESS) {
            return status;
        }

        FILE *fp=my_fopen(correct_outputfile,"r");
        if(fp == NULL) {
            free_results_mocks_s_mu(&results);
            return EXIT_FAILURE;
        }

        const int nmubin = results.nmu_bins;
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
    END_INTEGRATION_TEST_SECTION(free_results_mocks_s_mu(&results));

    /* If the test failed, then write the output into a temporary file */
    if(ret != EXIT_SUCCESS && results.nsbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        if(fp == NULL) {
            free_results_mocks_s_mu(&results);
            return EXIT_FAILURE;
        }
        const double dmu= 1.0/(double)results.nmu_bins;
        const int nmubin = results.nmu_bins;
        for(int i=1;i<results.nsbin;i++) {
            const double logrp = log10(results.supp[i]);
            for(int j=0;j<nmubin;j++) {
                int index = i*(nmubin+1) + j;
                fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf %20.8lf \n",results.npairs[index],results.savg[index],logrp,(j+1)*dmu, results.weightavg[index]);
            }
        }
        fclose(fp);
    }

    free_results_mocks_s_mu(&results);
    return ret;
}

int test_DDtheta_mocks(const char *correct_outputfile)
{
    int autocorr = (RA1==RA2) ? 1:0;
    int ret = EXIT_FAILURE;
    results_countpairs_theta results;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    BEGIN_DDTHETA_INTEGRATION_TEST_SECTION;
        int status = countpairs_theta_mocks(ND1,RA1,DEC1,
                                            ND2,RA2,DEC2,
                                            nthreads,
                                            autocorr,
                                            angular_binfile,
                                            &results,
                                            &options,
                                            &extra);

        if(status != EXIT_SUCCESS) {
            return status;
        }

        /*---Output-Pairs-------------------------------------*/
        FILE *fp=my_fopen(correct_outputfile,"r");
        if(fp == NULL) {
            free_results_countpairs_theta(&results);
            return EXIT_FAILURE;
        }
        for(int i=1;i<results.nbin;i++) {
            uint64_t npairs;
            double theta_avg, weightavg;
            ret = EXIT_FAILURE;
            int nitems = fscanf(fp,"%"SCNu64" %lf %*f %*f %lf%*[^\n]", &npairs, &theta_avg, &weightavg);
            if(nitems != 3) {
                ret = EXIT_FAILURE;//not required but showing intent
                break;
            }
            int floats_equal = AlmostEqualRelativeAndAbs_double(theta_avg, results.theta_avg[i], maxdiff, maxreldiff);
            int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[i], maxdiff, maxreldiff);

            //Check for exact equality of npairs and float "equality" for theta_avg
            if(npairs == results.npairs[i] && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
                ret = EXIT_SUCCESS;
            } else {
                ret = EXIT_FAILURE;//not required but showing intent
                if(npairs != results.npairs[i]) {
                    fprintf(stderr,"Failed (bin #%2d). True npairs = %"PRIu64 " Computed npairs = %"PRIu64"\n", i, npairs, results.npairs[i]);
                }
                if(floats_equal != EXIT_SUCCESS) {
                    fprintf(stderr,"Failed (bin #%2d). True thetaavg = %e Computed thetaavg = %e\n", i, theta_avg, results.theta_avg[i]);
                }
                if(weights_equal != EXIT_SUCCESS) {
                    fprintf(stderr,"Failed (bin #%2d). True weightavg = %e Computed weightavg = %e.\n",
                            i, weightavg, results.weightavg[i]);
                }
                /* break; */
            }
        }
        fclose(fp);
        END_DDTHETA_INTEGRATION_TEST_SECTION(free_results_countpairs_theta(&results));

    if(ret != EXIT_SUCCESS && results.nbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        double theta_low = results.theta_upp[0];
        for(int i=1;i<results.nbin;i++) {
            fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf\n",
                    results.npairs[i],results.theta_avg[i],theta_low,results.theta_upp[i], results.weightavg[i]);
            theta_low=results.theta_upp[i];
        }
        fclose(fp);
    }

    //free the result structure
    free_results_countpairs_theta(&results);
    return ret;
}

int test_vpf_mocks(const char *correct_outputfile)
{
    const double rmax=10.0;
    const int nbin=10;
    const int nc=10000;
    const int num_pN=6;
    const int64_t Nran=nc;//Hack. Need to set it to nc so that the loop runs
    double *xran=NULL,*yran=NULL,*zran=NULL;
    const int threshold_neighbors=1;
    const char centers_file[]="../tests/data/Mr19_centers_xyz_forVPF_rmax_10Mpc.txt";
    results_countspheres_mocks results;
    int ret = EXIT_FAILURE;

    BEGIN_INTEGRATION_TEST_SECTION
        int status = countspheres_mocks(ND1, RA1, DEC1, CZ1,
                                        Nran, xran, yran, zran,
                                        threshold_neighbors,
                                        rmax, nbin, nc,
                                        num_pN,
                                        centers_file,
                                        cosmology_flag,
                                        &results,
                                        &options, NULL);
        if(status != EXIT_SUCCESS) {
            return status;
        }


        //Output the results
        FILE *fp=my_fopen(correct_outputfile,"r");
        if(fp == NULL) {
            free_results_countspheres_mocks(&results);
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
    END_INTEGRATION_TEST_SECTION(free_results_countspheres_mocks(&results));

    /* If the test failed, then write the output into a temporary file */
    if(ret != EXIT_SUCCESS && results.nbin > 0) {
        FILE *fp=my_fopen(tmpoutputfile,"w");
        const double rstep = rmax/(double)nbin ;
        for(int ibin=0;ibin<results.nbin;ibin++) {
            const double r=(ibin+1)*rstep;
            fprintf(fp,"%10.2lf ", r);
            for(int i=0;i<num_pN;i++) {
                fprintf(fp," %10.4e", (results.pN)[ibin][i]);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }

    //free the result structure
    free_results_countspheres_mocks(&results);
    return ret;
}


void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat)
{
    int free_RA2=0;
    if(RA2 != NULL && RA2 != RA1) {
        free_RA2=1;
    }


    //Check to see if data has to be read for RA1/DEC1/CZ1
    if (strcmp(current_file1,firstfilename) != 0) {
        /* fprintf(stderr,"Freeing the first data-set and replacing with data from file `%s'\n",firstfilename); */
        //replace the data-set
        if(RA1 != NULL) {
            free(RA1);free(DEC1);free(CZ1);free(weights1);
            if(free_RA2 == 0) {
                RA2  = NULL;
                DEC2 = NULL;
                CZ2  = NULL;
                weights2 = NULL;
            }
        }
        ND1 = read_positions(firstfilename,firstformat, sizeof(double), 4, &RA1, &DEC1, &CZ1, &weights1);
    }

    //first check if only one unique file is asked for
    if(strncmp(firstfilename,secondfilename,strlen(firstfilename))==0) {
        //But RA2 might have read-in a different file->avoid a memory-leak
        if(free_RA2 == 1) {
            free(RA2);free(DEC2);free(CZ2);free(weights2);
            free_RA2 = 0;//not essential since the code returns after this section
        }
        /* fprintf(stderr,"Second data-set is the same as the first data-set. First file = `%s' and second file = `%s'\n",firstfilename,secondfilename); */
        RA2=RA1;
        DEC2=DEC1;
        CZ2=CZ1;
        weights2=weights1;
        ND2=ND1;
        strncpy(current_file2,secondfilename,MAXLEN);
        return;
    }


    //Check to see if data has to be read for RA2/DEC2/CZ2
    if (strncmp(current_file2,secondfilename,strlen(current_file2)) != 0 || RA2 == NULL) {
        //replace the data-set
        if(free_RA2 == 1) {
            free(RA2);free(DEC2);free(CZ2);free(weights2);
        }
        /* fprintf(stderr,"Second data-set is different -- reading in the new data-set from `%s'\n",secondfilename); */
        ND2 = read_positions(secondfilename,secondformat, sizeof(double), 4, &RA2, &DEC2, &CZ2, &weights2);
        strncpy(current_file2,secondfilename,MAXLEN);
    }
}


int main(int argc, char **argv)
{
    struct timeval tstart,t0,t1;
    char file[]="../tests/data/Mr19_mock_northonly.rdcz.dat";
    char fileformat[]="a";

    options = get_config_options();
    options.need_avg_sep=1;
    options.verbose=0;
    options.periodic=0;
    options.float_type=sizeof(double);
    options.fast_divide_and_NR_steps=0;
    options.fast_acos=0;
    options.link_in_ra = 1;
    options.link_in_dec = 1;
    options.copy_particles = 1;

    int status = init_cosmology(cosmology_flag);
    if(status != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }
    gettimeofday(&tstart,NULL);

    //set the globals.
    ND1 = read_positions(file,fileformat, sizeof(double), 4, &RA1, &DEC1, &CZ1, &weights1);

    ND2 = ND1;
    RA2 = RA1;
    DEC2 = DEC1;
    CZ2 = CZ1;
    weights2 = weights1;

    strncpy(current_file1,file,MAXLEN);
    strncpy(current_file2,file,MAXLEN);
    reset_bin_refine_factors(&options);

    int failed=0;

    const char alltests_names[][MAXLEN] = {"Mr19 mocks DDrppi (DD)",
                                           "Mr19 mocks wtheta (DD)",
                                           "Mr19 mocks vpf (data)",
                                           "Mr19 mocks DDrppi (DR)",
                                           "Mr19 mocks wtheta (DR)",
                                           "Mr19 mocks vpf (randoms)",
                                           "Mr19 mocks DDsmu (RR)",
                                           "Mr19 mocks DDsmu (DR)"};
    const int ntests = sizeof(alltests_names)/(sizeof(char)*MAXLEN);
    const int function_pointer_index[] = {0,1,2,0,1,2,3,3};//0->DDrppi, 1->wtheta, 2->vpf, 3->DDsmu
    assert(sizeof(function_pointer_index)/sizeof(int) == ntests && "Number of tests should equal the number of functions");

    const char correct_outputfiles[][MAXLEN] = {"../tests/Mr19_mock.DD", /* Test 0 Mr19 DD */
                                                "../tests/Mr19_mock_wtheta.DD", /* Test 1 Mr19 wtheta DD*/
                                                "../tests/Mr19_mock_vpf", /* Test 2 Mr19 mocks vpf */
                                                "../tests/Mr19_mock.DR", /* Test 3 Mr19 DR */
                                                "../tests/Mr19_mock_wtheta.DR", /* Test 4 Mr19 wtheta DR */
                                                "../tests/Mr19_randoms_vpf", /* Test 5 Mr19 randoms vpf */
                                                "../tests/Mr19_mock_DDsmu.RR", /* Test 6 Mr19 RR smu */
                                                "../tests/Mr19_mock_DDsmu.DR"}; /* Test 7 Mr19 DR smu */
    const char firstfilename[][MAXLEN] = {"../tests/data/Mr19_mock_northonly.rdcz.dat",
                                          "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                          "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff"};
    const char firstfiletype[][MAXLEN]  = {"a","a","a","f","f","f","f","f"};
    const char secondfilename[][MAXLEN] = {"../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                           "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat"};
    const char secondfiletype[][MAXLEN] = {"a","a","a","a","a","f","f","a"};

    const double allpimax[]             = {40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0};

    int (*allfunctions[]) (const char *) = {test_DDrppi_mocks,test_DDtheta_mocks,test_vpf_mocks,test_DDsmu_mocks};
    const int numfunctions=4;//4 functions total

    int total_tests=0,skipped=0;

    if(argc==1) {
        //nothing was passed at the command-line run all tests
        for(int i=0;i<ntests;i++) {
            int function_index = function_pointer_index[i];
            assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
            const char *testname = alltests_names[i];
            int skip_test=test_all_files_present(2,firstfilename[i],secondfilename[i]);
            if(skip_test != 0) {
                fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". Test data-file(s) (`%s',`%s') not found\n", testname, firstfilename[i], secondfilename[i]);
                skipped++;
                continue;
            }
            read_data_and_set_globals(firstfilename[i],firstfiletype[i],secondfilename[i],secondfiletype[i]);
            pimax=allpimax[i];
            gettimeofday(&t0,NULL);
            status = (*allfunctions[function_index])(correct_outputfiles[i]);
            gettimeofday(&t1,NULL);
            double pair_time=ADD_DIFF_TIME(t0,t1);
            total_tests++;
            if(status==EXIT_SUCCESS) {
                fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                run_system_call(execstring);//can ignore the status here

            } else {
                fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                failed++;
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,i);
                run_system_call(execstring);//can ignore the status here
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
                    fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". Test data-file(s) (`%s',`%s') not found\n", testname, firstfilename[this_test_num], secondfilename[this_test_num]);
                    skipped++;
                    continue;
                }

                total_tests++;
                read_data_and_set_globals(firstfilename[this_test_num],firstfiletype[this_test_num],secondfilename[this_test_num],secondfiletype[this_test_num]);
                pimax=allpimax[this_test_num];
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
            } else {
                fprintf(stderr,ANSI_COLOR_YELLOW "WARNING: Test = %d is not a valid test index. Valid test indices range between [0,%d] " ANSI_COLOR_RESET "\n",this_test_num,ntests-1);
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

    if(RA2 != RA1) {
        free(RA2);free(DEC2);free(CZ2);
    }
    free(RA1);free(DEC1);free(CZ1);
    return EXIT_SUCCESS;
}
