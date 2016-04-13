/* File: tests_mocks.c */
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

#ifndef MAXLEN
#define MAXLEN 500
#endif

#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif

#ifndef OUTPUT_RPAVG
#define OUTPUT_RPAVG
#endif

#ifndef OUTPUT_THETAAVG
#define OUTPUT_THETAAVG
#endif

#ifdef FAST_DIVIDE
#undef FAST_DIVIDE
#endif

#if !(defined(__INTEL_COMPILER)) && defined(USE_AVX)
#warning Test suite for mocks are faster with Intel compiler, icc, AVX libraries.
#endif

#ifndef SILENT
#define SILENT
#endif

#include "function_precision.h"
#include "io.h"
#include "defs.h"
#include "utils.h"


//Including the C files directly
#include "gridlink_mocks.c"
#include "io.c"
#include "ftread.c"
#include "../DDrppi/countpairs_rp_pi_mocks.c"
#include "../wtheta/countpairs_theta_mocks.c"
#include "../vpf/countspheres_mocks.c"


char tmpoutputfile[]="./tests_mocks_output.txt";

int test_DDrppi_mocks(const char *correct_outputfile);
int test_wtheta_mocks(const char *correct_outputfile);
int test_vpf_mocks(const char *correct_outputfile);

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat);

//Global variables
int ND1;
DOUBLE *RA1=NULL,*DEC1=NULL,*CZ1=NULL;

int ND2;
DOUBLE *RA2=NULL,*DEC2=NULL,*CZ2=NULL;

char binfile[]="../tests/bins";
char angular_binfile[]="../tests/angular_bins";
DOUBLE pimax=40.0;
double boxsize=420.0;
#if defined(USE_OMP) && defined(_OPENMP)
const int nthreads=4;
#endif
const int cosmology_flag=1;
char current_file1[MAXLEN],current_file2[MAXLEN];

//end of global variables


int test_DDrppi_mocks(const char *correct_outputfile)
{
    assert(RA1 != NULL && DEC1 != NULL && CZ1 != NULL && "ERROR: In test suite for DDrppi ra/dec/cz can not be NULL pointers");
    int autocorr = (RA1==RA2) ? 1:0;

    //Do DD(rp,pi) counts
    results_countpairs_mocks *results  = countpairs_mocks(ND1,RA1,DEC1,CZ1,
                                                          ND2,RA2,DEC2,CZ2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                          nthreads,
#endif
                                                          autocorr,
                                                          binfile,
                                                          pimax,
                                                          cosmology_flag);

    FILE *fp=my_fopen(tmpoutputfile,"w");
    const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
    const int npibin = results->npibin;
    for(int i=1;i<results->nbin;i++) {
        const double logrp = LOG10(results->rupp[i]);
        for(int j=0;j<npibin;j++) {
            int index = i*(npibin+1) + j;
            fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->rpavg[index],logrp,(j+1)*dpi);
        }
    }
    fclose(fp);

    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

    free_results_mocks(&results);
    return ret;
}

int test_wtheta_mocks(const char *correct_outputfile)
{
    int autocorr = (RA1==RA2) ? 1:0;

    results_countpairs_theta *results = countpairs_theta_mocks(ND1,RA1,DEC1,
                                                               ND2,RA2,DEC2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                               nthreads,
#endif
                                                               autocorr,
                                                               angular_binfile) ;

    /*---Output-Pairs-------------------------------------*/
    FILE *fp=my_fopen(tmpoutputfile,"w");
    DOUBLE theta_low = results->theta_upp[0];
    for(int i=1;i<results->nbin;i++) {
        fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->theta_avg[i],theta_low,results->theta_upp[i]);
        theta_low=results->theta_upp[i];
    }
    fclose(fp);

    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

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
    DOUBLE *xran=NULL,*yran=NULL,*zran=NULL;
    const int threshold_neighbors=1;
    const char centers_file[]="../tests/data/Mr19_centers_xyz_forVPF_rmax_10Mpc.txt";

    results_countspheres_mocks *results = countspheres_mocks(ND1, RA1, DEC1, CZ1,
                                                             Nran, xran, yran, zran,
                                                             threshold_neighbors,
                                                             rmax, nbin, nc,
                                                             num_pN,
                                                             centers_file,
                                                             cosmology_flag);



    //Output the results
    FILE *fp=my_fopen(tmpoutputfile,"w");
    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    for(int ibin=0;ibin<results->nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        fprintf(fp,"%10.2"DOUBLE_FORMAT" ", r);
        for(int i=0;i<num_pN;i++) {
            fprintf(fp," %10.4e", (results->pN)[ibin][i]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);
    assert(ret == EXIT_SUCCESS);
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
    if (strncmp(current_file1,firstfilename,strlen(current_file1)) != 0) {
        /* fprintf(stderr,"Freeing the first data-set and replacing with data from file `%s'\n",firstfilename); */
        //replace the data-set
        if(RA1 != NULL) {
            free(RA1);free(DEC1);free(CZ1);
            if(free_RA2 == 0) {
                RA2  = NULL;
                DEC2 = NULL;
                CZ2  = NULL;
            }
        }
        ND1 = read_positions(firstfilename,firstformat, sizeof(DOUBLE), 3, &RA1, &DEC1, &CZ1);
        strncpy(current_file1,firstfilename,MAXLEN);
    }

    //first check if only one unique file is asked for
    if(strncmp(firstfilename,secondfilename,strlen(firstfilename))==0) {
        //But RA2 might have read-in a different file->avoid a memory-leak
        if(free_RA2 == 1) {
            free(RA2);free(DEC2);free(CZ2);
            free_RA2 = 0;//not essential since the code returns after this section
        }
        /* fprintf(stderr,"Second data-set is the same as the first data-set. First file = `%s' and second file = `%s'\n",firstfilename,secondfilename); */
        RA2=RA1;
        DEC2=DEC1;
        CZ2=CZ1;
        ND2=ND1;
        strncpy(current_file2,secondfilename,MAXLEN);
        return;
    }


    //Check to see if data has to be read for RA2/DEC2/CZ2
    if (strncmp(current_file2,secondfilename,strlen(current_file2)) != 0 || RA2 == NULL) {
        //replace the data-set
        if(free_RA2 == 1) {
            free(RA2);free(DEC2);free(CZ2);
        }
        /* fprintf(stderr,"Second data-set is different -- reading in the new data-set from `%s'\n",secondfilename); */
        ND2 = read_positions(secondfilename,secondformat, sizeof(DOUBLE), 3, &RA2, &DEC2, &CZ2);
        strncpy(current_file2,secondfilename,MAXLEN);
    }
}


int main(int argc, char **argv)
{
    struct timeval tstart,t0,t1;
    char file[]="../tests/data/Mr19_mock_northonly.rdcz.dat";
    char fileformat[]="a";

#ifdef PERIODIC
#error PERIODIC must not be defined for running non-periodic tests
#endif

    init_cosmology(cosmology_flag);
    gettimeofday(&tstart,NULL);

    //set the globals.
    ND1 = read_positions(file,fileformat, sizeof(DOUBLE), 3, &RA1, &DEC1, &CZ1);
    ND2 = ND1;
    RA2 = RA1;
    DEC2 = DEC1;
    CZ2 = CZ1;

    strncpy(current_file1,file,MAXLEN);
    strncpy(current_file2,file,MAXLEN);

    int failed=0;
    int status;

    const char alltests_names[][MAXLEN] = {"Mr19 mocks DDrppi (DD)","Mr19 mocks wtheta (DD)","Mr19 mocks vpf (data)","Mr19 mocks DDrppi (DR)", "Mr19 mocks wtheta (DR)","Mr19 mocks vpf (randoms)"};
    const int ntests = sizeof(alltests_names)/(sizeof(char)*MAXLEN);
    const int function_pointer_index[] = {0,1,2,0,1,2};//0->DDrppi, 1->wtheta, 2->vpf
    assert(sizeof(function_pointer_index)/sizeof(int) == ntests && "Number of tests should equal the number of functions");

    const char correct_outoutfiles[][MAXLEN] = {"../tests/Mr19_mock.DD", /* Test 0 Mr19 DD */
                                                "../tests/Mr19_mock_wtheta.DD", /* Test 1 Mr19 wtheta DD*/
                                                "../tests/Mr19_mock_vpf", /* Test 2 Mr19 mocks vpf */
                                                "../tests/Mr19_mock.DR", /* Test 3 Mr19 DR */
                                                "../tests/Mr19_mock_wtheta.DR", /* Test 4 Mr19 wtheta DR */
                                                "../tests/Mr19_randoms_vpf"}; /* Test 5 Mr19 randoms vpf */
    const char firstfilename[][MAXLEN] = {"../tests/data/Mr19_mock_northonly.rdcz.dat",
                                          "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                          "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff",
                                          "../tests/data/Mr19_randoms_northonly.rdcz.ff"};
    const char firstfiletype[][MAXLEN]  = {"a","a","a","f","f","f"};
    const char secondfilename[][MAXLEN] = {"../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_mock_northonly.rdcz.dat",
                                           "../tests/data/Mr19_randoms_northonly.rdcz.ff"};
    const char secondfiletype[][MAXLEN] = {"a","a","a","a","a","f"};
    const DOUBLE allpimax[]             = {40.0,40.0,40.0,40.0,40.0,40.0};

    int (*allfunctions[]) (const char *) = {test_DDrppi_mocks,test_wtheta_mocks,test_vpf_mocks};
    const int numfunctions=3;//3 functions total

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
            status = (*allfunctions[function_index])(correct_outoutfiles[i]);
            gettimeofday(&t1,NULL);
            double pair_time=ADD_DIFF_TIME(t0,t1);
            total_tests++;
            if(status==EXIT_SUCCESS) {
                fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                run_system_call(execstring);
            } else {
                fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                failed++;
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,i);
                run_system_call(execstring);

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
                status = (*allfunctions[function_index])(correct_outoutfiles[this_test_num]);
                gettimeofday(&t1,NULL);
                double pair_time = ADD_DIFF_TIME(t0,t1);
                if(status==EXIT_SUCCESS) {
                    fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                    char execstring[MAXLEN];
                    my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                    run_system_call(execstring);
                } else {
                    fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                    failed++;
                    char execstring[MAXLEN];
                    my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,this_test_num);
                    run_system_call(execstring);
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
