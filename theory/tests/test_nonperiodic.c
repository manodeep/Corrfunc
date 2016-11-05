/* File: test_nonperiodic.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <inttypes.h>

#ifndef MAXLEN
#define MAXLEN 500
#endif

#include "defs.h"
#include "utils.h"
#include "io.h"
#include "ftread.h"

#include "../xi_of_r/countpairs.h"
#include "../xi_rp_pi/countpairs_rp_pi.h"

char tmpoutputfile[]="./test_nonperiodic_output.txt";

int test_nonperiodic_DD(const char *correct_outputfile);
int test_nonperiodic_DDrppi(const char *correct_outputfile);
void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat,
                               const char *firstweightsfilename, const char *firstweightsfileformat,
                               const char *secondweightsfilename, const char *secondweightsfileformat);

//Global variables
int ND1;
double *X1=NULL,*Y1=NULL,*Z1=NULL,*weights1=NULL;

int ND2;
double *X2=NULL,*Y2=NULL,*Z2=NULL,*weights2=NULL;

char binfile[]="bins";
double pimax=40.0;
double boxsize=420.0;
#if defined(_OPENMP)
const int nthreads=4;
#else
const int nthreads=1;
#endif

char current_file1[MAXLEN],current_file2[MAXLEN];
struct config_options options;

const double maxdiff = 1e-9;
const double maxreldiff = 1e-6;
//end of global variables

int test_nonperiodic_DD(const char *correct_outputfile)
{
    int autocorr = (X1==X2) ? 1:0;
    
    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    //Do the straight-up DD counts
    results_countpairs results;
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

    int ret = EXIT_FAILURE;
    FILE *fp=my_fopen(correct_outputfile,"r");
    if(fp == NULL) {
        free_results(&results);
        return EXIT_FAILURE;
    }
    for(int i=1;i<results.nbin;i++) {
        uint64_t npairs;
        double rpavg, weightavg;
        ret = EXIT_FAILURE;
        int nitems = fscanf(fp,"%"SCNu64" %lf %*f %*f %lf%*[^\n]", &npairs, &rpavg, &weightavg);
        if(nitems != 3) {
            break;
        }
        int floats_equal = AlmostEqualRelativeAndAbs_double(rpavg, results.rpavg[i], maxdiff, maxreldiff);
        int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[i], maxdiff, maxreldiff);

        //Check for exact equality of npairs and float "equality" for rpavg
        if(npairs == results.npairs[i] && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
            ret = EXIT_SUCCESS;
        } else {
            ret = EXIT_FAILURE;//not required but showing intent 
            fprintf(stderr,"True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[i]);
            fprintf(stderr,"True rpavg  = %e Computed rpavg = %e. floats_equal = %d\n", rpavg, results.rpavg[i], floats_equal);
            fprintf(stderr,"True weightavg = %e Computed weightavg = %e. weights_equal = %d\n", weightavg, results.weightavg[i], weights_equal);
            break;
        }
    }
    fclose(fp);

    /* If the test failed, then write to temporary file, so a comparison can be made */
    if(ret != EXIT_SUCCESS) {
        fp=my_fopen(tmpoutputfile,"w"); 
        double rlow = results.rupp[0];
        for(int i=1;i<results.nbin;i++) {
            fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf \n",results.npairs[i],results.rpavg[i],rlow,results.rupp[i], results.weightavg[i]);
            rlow=results.rupp[i];
        }
        fclose(fp);
    }

    free_results(&results);
    return ret;
}

int test_nonperiodic_DDrppi(const char *correct_outputfile)
{
    int autocorr = (X1==X2) ? 1:0;
    
    // Set up the weights pointers
    weight_method_t weight_method = PAIR_PRODUCT;
    struct extra_options extra = get_extra_options(weight_method);
    extra.weights0.weights[0] = weights1;
    extra.weights1.weights[0] = weights2;

    results_countpairs_rp_pi results;
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

    int ret = EXIT_FAILURE;
    const int npibin = results.npibin;
    const double dpi = pimax/(double)results.npibin ;
    FILE *fp=my_fopen(correct_outputfile,"r");
    if(fp == NULL) {
        free_results_rp_pi(&results);
        return EXIT_FAILURE;
    }

    
    for(int i=1;i<results.nbin;i++) {
        for(int j=0;j<npibin;j++) {
            int index = i*(npibin+1) + j;
            uint64_t npairs;
            double rpavg, weightavg;
            ret = EXIT_FAILURE;
            int nitems = fscanf(fp,"%"SCNu64" %lf %*f %*f %lf%*[^\n]", &npairs, &rpavg, &weightavg);
            if(nitems != 3) {
                i = results.nbin;
                ret = EXIT_FAILURE;//not required but showing intent
                break;
            }
            int floats_equal = AlmostEqualRelativeAndAbs_double(rpavg, results.rpavg[index], maxdiff, maxreldiff);
            int weights_equal = AlmostEqualRelativeAndAbs_double(weightavg, results.weightavg[index], maxdiff, maxreldiff);
            
            //Check for exact equality of npairs and float "equality" for rpavg
            if(npairs == results.npairs[index] && floats_equal == EXIT_SUCCESS && weights_equal == EXIT_SUCCESS) {
                ret = EXIT_SUCCESS;
            } else {
                fprintf(stderr,"Failed. True npairs = %"PRIu64 " Computed results npairs = %"PRIu64"\n", npairs, results.npairs[index]);
                fprintf(stderr,"Failed. True rpavg = %e Computed rpavg = %e. floats_equal = %d\n", rpavg, results.rpavg[index], floats_equal);
                fprintf(stderr,"Failed. True weightavg = %e Computed weightavg = %e. weights_equal = %d\n", weightavg, results.weightavg[index], weights_equal);
                ret = EXIT_FAILURE;//not required but showing intent 
                i=results.nbin;
                break;
            }
        }
    }
    fclose(fp);

    if(ret != EXIT_SUCCESS) {
        fp = my_fopen(tmpoutputfile,"w");
        for(int i=1;i<results.nbin;i++) {
            const double logrp = log10(results.rupp[i]);
            for(int j=0;j<npibin;j++) {
                int index = i*(npibin+1) + j;
                fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf %20.8lf\n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi, results.weightavg[index]);
            }
        }
        fclose(fp);
    }

    //free the result structure
    free_results_rp_pi(&results);
    return ret;
}

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat,
                               const char *firstweightsfilename, const char *firstweightsfileformat,
                               const char *secondweightsfilename, const char *secondweightsfileformat)
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
        ND1 = read_positions(firstfilename,firstformat, sizeof(double), 3, &X1, &Y1, &Z1);
        strncpy(current_file1,firstfilename,MAXLEN);
        int64_t wND1 = read_columns_into_array(firstweightsfilename, firstweightsfileformat, sizeof(double), 1, (void **) &weights1);
        if(wND1 != ND1){
            fprintf(stderr, "Error: number of weights read from %s did not match number of positions read from %s\n", firstweightsfilename, firstfilename);
            return;
        }
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
        ND2 = read_positions(secondfilename,secondformat, sizeof(double), 3, &X2, &Y2, &Z2);
        strncpy(current_file2,secondfilename,MAXLEN);
        int64_t wND2 = read_columns_into_array(secondweightsfilename, secondweightsfileformat, sizeof(double), 1, (void **) &weights2);
        if(wND2 != ND2){
            fprintf(stderr,"Error: number of weights read from %s did not match number of positions read from %s\n", secondweightsfilename, secondfilename);
            return;
        }
    }
}


int main(int argc, char **argv)
{
    struct timeval tstart,t0,t1;
    char file[]="../tests/data/gals_Mr19.ff";
    char fileformat[]="f";
    
    char weights_file[] = "../tests/data/gals_Mr19_weights_random.csv";
    char weights_fileformat[] = "c";
    
    options = get_config_options();
    options.need_avg_sep=1;
    options.verbose=0;
    options.periodic=0;
    options.float_type=sizeof(double);
    //options.instruction_set = FALLBACK;

    gettimeofday(&tstart,NULL);

    //set the globals
    ND1 = read_positions(file,fileformat, sizeof(double), 3, &X1, &Y1, &Z1);
    int64_t wND1 = read_columns_into_array(weights_file, weights_fileformat, sizeof(double), 1, (void **) &weights1);
    XASSERT(wND1 == ND1, "Error: number of weights read from %s did not match number of positions read from %s\n", weights_file, file);
    
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

    const char alltests_names[][MAXLEN] = {"Mr19 DD (nonperiodic)","Mr19 DDrppi (nonperiodic)","CMASS DDrppi DR (nonperiodic)"};
    const int ntests = sizeof(alltests_names)/(sizeof(char)*MAXLEN);
    const int function_pointer_index[] = {0,1,1};//0->DD, 1->DDrppi

    const char correct_outputfiles[][MAXLEN] = {"Mr19_DD_nonperiodic","Mr19_DDrppi_nonperiodic","cmass_DR_nonperiodic"};
    const char firstfilename[][MAXLEN] = {"../tests/data/gals_Mr19.ff","../tests/data/gals_Mr19.ff","../tests/data/cmassmock_Zspace.ff"};
    const char firstfiletype[][MAXLEN] = {"f","f","f"};
    const char secondfilename[][MAXLEN] = {"../tests/data/gals_Mr19.ff","../tests/data/gals_Mr19.ff","../tests/data/random_Zspace.ff"};
    const char secondfiletype[][MAXLEN] = {"f","f","f"};
    
    const char firstweightsfilename[][MAXLEN] = {"../tests/data/gals_Mr19_weights_random.csv","../tests/data/gals_Mr19_weights_random.csv","../tests/data/cmassmock_Zspace_weights_random.csv"};
    const char firstweightsfiletype[][MAXLEN] = {"c","c","c"};
    const char secondweightsfilename[][MAXLEN] = {"../tests/data/gals_Mr19_weights_random.csv","../tests/data/gals_Mr19_weights_random.csv","../tests/data/random_Zspace_weights_random.csv"};
    const char secondweightsfiletype[][MAXLEN] = {"c","c","c"};

    const double allpimax[]             = {40.0,40.0,80.0};

    int (*allfunctions[]) (const char *) = {test_nonperiodic_DD,test_nonperiodic_DDrppi};
    const int numfunctions=2;//2 functions total

    int total_tests=0,skipped=0;

    if(argc==1) {
        //nothing was passed at the command-line run all tests
        for(int i=0;i<ntests;i++) {
            int function_index = function_pointer_index[i];
            assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
            const char *testname = alltests_names[i];
            int skip_test=test_all_files_present(4,firstfilename[i],secondfilename[i],firstweightsfilename[i],secondweightsfilename[i]);
            if(skip_test != 0) {
                fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". Test data-file(s) (`%s',`%s') not found\n", testname, firstfilename[i], secondfilename[i]);
                skipped++;
                continue;
            }

            read_data_and_set_globals(firstfilename[i],firstfiletype[i],secondfilename[i],secondfiletype[i],firstweightsfilename[i], firstweightsfiletype[i], secondweightsfilename[i], secondweightsfiletype[i]);
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
                int skip_test=test_all_files_present(4,firstfilename[this_test_num],secondfilename[this_test_num],firstweightsfilename[this_test_num],secondweightsfilename[this_test_num]);
                if(skip_test != 0) {
                    fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". Test data-file(s) (`%s',`%s') not found\n", testname, firstfilename[this_test_num], secondfilename[this_test_num]);
                    skipped++;
                    continue;
                }
                total_tests++;
                read_data_and_set_globals(firstfilename[this_test_num],firstfiletype[this_test_num],secondfilename[this_test_num],secondfiletype[this_test_num],firstweightsfilename[this_test_num], firstweightsfiletype[this_test_num], secondweightsfilename[this_test_num], secondweightsfiletype[this_test_num]);
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
