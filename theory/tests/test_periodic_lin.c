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

char tmpoutputfile[]="./test_periodic_lin_output.txt";
char binfile_lin[]="../tests/bins_lin";

int test_periodic_DD(const char *correct_outputfile);
void read_data_and_set_globals(const char *firstfilename, const char *firstformat);


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

int main(int argc, char **argv)
{
    (void) argc;
    (void) argv;
    
    struct timeval tstart,t0,t1;
    options = get_config_options();
    options.need_avg_sep=1;
    options.verbose=1;
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

    strncpy(current_file1,file,MAXLEN);
    reset_bin_refine_factors(&options);

    int failed=0;
    int status;

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
    const int numfunctions=1;//6 functions total

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
