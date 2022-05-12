/* File: test_periodic.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/
//#define INTEGRATION_TESTS

#include "tests_common.h"
#include "io.h"

#include "../DD/countpairs.h"

char tmpoutputfile[]="./test_periodic_lin_output.txt";
char binfile_lin[]="../tests/bins_lin";

int test_periodic_DD(const char *correct_outputfile);
void read_data_and_set_globals(const char *firstfilename, const char *firstformat);
int write_bins_to_file(const double rmin, const double rmax, const double nbins, const char *linear_binfile);
int compare_two_results(const results_countpairs *results_reference, const results_countpairs *results_test);
int test_custom_and_linear_bins(void);

//Global variables
int64_t ND1;
double *X1=NULL,*Y1=NULL,*Z1=NULL,*weights1=NULL;

char current_file1[MAXLEN+1];


struct config_options options;
//end global variables

int test_periodic_DD(const char *correct_outputfile)
{
    int autocorr = 1;
    results_countpairs results;
    int ret = EXIT_FAILURE;

    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;

    BEGIN_INTEGRATION_TEST_SECTION

        //Do the straight-up DD counts
        int status = countpairs(ND1,X1,Y1,Z1,
                                ND1,X1,Y1,Z1,
                                nthreads,
                                autocorr,
                                binfile_lin,
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
            int npairs_equal = npairs == results.npairs[i] ? EXIT_SUCCESS : EXIT_FAILURE;
            int floats_equal = AlmostEqualRelativeAndAbs_double(rpavg, results.rpavg[i], maxdiff, maxreldiff);
            int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[i], maxdiff, maxreldiff);

            //Check for exact equality of npairs and float "equality" for rpavg
            if(npairs_equal == EXIT_SUCCESS && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
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

void read_data_and_set_globals(const char *firstfilename, const char *firstformat)
{
    //Check to see if data has to be read for X1/Y1/Z1
    if (strncmp(current_file1,firstfilename,strlen(current_file1)) != 0) {
        //replace the data-set
        if(X1 != NULL) {
            free(X1);free(Y1);free(Z1);free(weights1);
        }

        ND1 = read_positions(firstfilename,firstformat, sizeof(double), 4, &X1, &Y1, &Z1, &weights1);
        strncpy(current_file1,firstfilename,MAXLEN);
    }
}

int write_bins_to_file(const double rmin, const double rmax, const double nbins, const char *linear_binfile)
{
    FILE *fp = fopen(linear_binfile,"w");
    if(fp == NULL) {
        fprintf(stderr,"Failed to open file %s for writing\n",linear_binfile);
        return EXIT_FAILURE;
    }
    const double dr = (rmax-rmin)/nbins;
    for(int i=0;i<nbins;i++) {
        const double rp_low = rmin + i*dr;
        const double rp_upp = rmin + (i+1)*dr;
        fprintf(fp,"%a %a\n",rp_low,rp_upp);
    }
    fclose(fp);

    return EXIT_SUCCESS;
}

int compare_two_results(const results_countpairs *results_reference, const results_countpairs *results_test)
{
    int ret = 0;
    for(int i=1;i<results_reference->nbin;i++) {
        int npairs_equal = results_reference->npairs[i] == results_test->npairs[i] ? EXIT_SUCCESS : EXIT_FAILURE;
        int floats_equal = AlmostEqualRelativeAndAbs_double(results_reference->rpavg[i],
                                                            results_test->rpavg[i],  maxdiff, maxreldiff);
        int weights_equal = AlmostEqualRelativeAndAbs_double(results_reference->weightavg[i],
                                                             results_test->weightavg[i], maxdiff, maxreldiff);

        // Check for exact equality of npairs and float "equality" for rpavg
        if(npairs_equal == EXIT_SUCCESS
            && floats_equal == EXIT_SUCCESS
            && weights_equal == EXIT_SUCCESS) {
            continue;
        } else {
            ret++;
            if(npairs_equal != EXIT_SUCCESS){
                fprintf(stderr,"\n[nbin = %d] Failed. True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", i-1, results_reference->npairs[i], results_test->npairs[i]);
            }
            if(floats_equal != EXIT_SUCCESS) {
                fprintf(stderr,"[nbin = %d] Failed. True rpavg = %e Computed rpavg = %e. floats_equal = %d\n", i-1, results_reference->rpavg[i], results_test->rpavg[i], floats_equal);
            }
            if(weights_equal != EXIT_SUCCESS) {
                fprintf(stderr,"[nbin = %d] Failed. True weightavg = %e Computed weightavg = %e. weights_equal = %d\n", i-1, results_reference->weightavg[i], results_test->weightavg[i], weights_equal);
            }
            break;
        }
    }

    return ret;
}

int test_custom_and_linear_bins(void)
{
    int autocorr = 1;
    const double rmin = 0.1;
    const double rmax_lower = 1.0, rmax_upper=50.0, rmax_step=3.0;
    const double max_nbins = 20;
    const int nbins_step = 3;
    const char *linear_binfile = "../tests/test_custom_and_linear_bins";

    int ntests = 0, nfailed = 0;
    for(double rmax=rmax_lower;rmax<rmax_upper; rmax+=rmax_step){
        for(int nbins=1;nbins<max_nbins;nbins+=nbins_step) {
            int status = write_bins_to_file(rmin, rmax, nbins, linear_binfile);
            if (status != EXIT_SUCCESS) {
                return status;
            }
            for(int iset=num_instructions-1;iset>=0;iset--) {
                fprintf(stderr,"[rmin, rmax, nbins] = [%0.2f, %0.2f, %0d], isa = %10s ",rmin, rmax, nbins, isa_name[iset]);
                results_countpairs results_reference, results_linear;
                // Set up the weights pointers
                weight_method_t weight_method = PAIR_PRODUCT;
                struct extra_options extra = get_extra_options(weight_method);
                extra.weights0.weights[0] = weights1;
                options.instruction_set = valid_instruction_sets[iset];

                options.bin_type = BIN_CUSTOM;
                //Do the straight-up DD counts
                status = countpairs(ND1,X1,Y1,Z1,
                                    ND1,X1,Y1,Z1,
                                    nthreads,
                                    autocorr,
                                    linear_binfile,
                                    &results_reference,
                                    &options,
                                    &extra);
                if(status != EXIT_SUCCESS) {
                    return status;
                }

                options.bin_type = BIN_LIN;
                status = countpairs(ND1,X1,Y1,Z1,
                                    ND1,X1,Y1,Z1,
                                    nthreads,
                                    autocorr,
                                    linear_binfile,
                                    &results_linear,
                                    &options,
                                    &extra);
                if(status != EXIT_SUCCESS) {
                    return status;
                }
                ntests++;
                int ret = compare_two_results(&results_reference, &results_linear);
                if(ret != EXIT_SUCCESS) {
                    nfailed++;
                    fprintf(stderr,ANSI_COLOR_RED "FAILED" ANSI_COLOR_RESET"\n");
                } else {
                    fprintf(stderr,ANSI_COLOR_GREEN "PASSED" ANSI_COLOR_RESET"\n");
                }
            }
        }
    }
    fprintf(stderr,ANSI_COLOR_RESET"ntests run = %d nfailed = %d"ANSI_COLOR_RESET"\n", ntests, nfailed);

    return nfailed;
}


int main(int argc, char **argv)
{
    (void) argc;
    (void) argv;

    struct timeval tstart,t0,t1;
    options = get_config_options();
    options.need_avg_sep=1;
    options.verbose=0;
    options.periodic=1;
    options.fast_divide_and_NR_steps=0;
    options.copy_particles=1;
    options.float_type=sizeof(double);
    options.boxsize = boxsize;
    options.bin_type = BIN_LIN;

    char file[]="../tests/data/gals_Mr19.ff";
    char fileformat[]="f";

    gettimeofday(&tstart,NULL);

    //set the globals
    ND1 = read_positions(file,fileformat, sizeof(double), 4, &X1, &Y1, &Z1, &weights1);
    nthreads = get_nthreads_from_affinity();

    strncpy(current_file1,file,MAXLEN);
    reset_bin_refine_factors(&options);

    int failed=0;
    int status;

    //Test the linear bins by comparing the results obtained
    // with bin_type=BIN_CUSTOM and bin_type=BIN_LIN. This is
    // a comprehensive test running over a large set of
    // (rmax, nbins) for each instruction set.
#ifdef INTEGRATION_TESTS
    status = test_custom_and_linear_bins();
    if(status != EXIT_SUCCESS) {
        failed++;
        return failed;
    }
#endif

    const char alltests_names[][MAXLEN] = {"Mr19 DD (periodic)",
                                          };
    const int ntests = sizeof(alltests_names)/(sizeof(char)*MAXLEN);
    const int function_pointer_index[] = {0,};//0->DD, 1->DDrppi,2->wp, 3->vpf, 4->xi, 5->DDsmu

    const char correct_outputfiles[][MAXLEN] = {"Mr19_DD_periodic_lin",
                                                };
    const char firstfilename[][MAXLEN] = {"../tests/data/gals_Mr19.ff",
                                         };
    const char firstfiletype[][MAXLEN] = {"f",};
    const double allpimax[]             = {40.0,};
    const double allboxsize[]             = {boxsize,};

    int (*allfunctions[]) (const char *) = {test_periodic_DD,
                                            };
    const int numfunctions=1;

    int total_tests=0,skipped=0;

    for (int i = 0; i < ntests; i++){
        int function_index = function_pointer_index[i];
        assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
        const char *testname = alltests_names[i];
        int skip_test=test_all_files_present(1,firstfilename[i]);
        if(skip_test != 0) {
            fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". File(s) not found\n", testname);
            skipped++;
            continue;
        }
        read_data_and_set_globals(firstfilename[i],firstfiletype[i]);
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

    free(X1);free(Y1);free(Z1);free(weights1);
    return failed;
}
