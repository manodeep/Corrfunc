#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#ifndef MAXLEN
#define MAXLEN 500
#endif

#ifdef PERIODIC
#undef PERIODIC 
#endif

#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif

#ifndef OUTPUT_RPAVG
#define OUTPUT_RPAVG
#endif

#include "function_precision.h"
#include "io.h"
#include "defs.h"
#include "utils.h"


//Including the C files directly
#include "../utils/gridlink.c"
#include "../io/io.c"
#include "../io/ftread.c"
#include "../xi_of_r/countpairs.c"
#include "../xi_rp_pi/countpairs_rp_pi.c"
#include "../wp/countpairs_wp.c"

char tmpoutputfile[]="./tmp_output.txt";

int test_nonperiodic_DD(void);
int test_nonperiodic_DDrppi(void);


//Global variables
int ND1;
DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;

int ND2;
DOUBLE *x2=NULL,*y2=NULL,*z2=NULL;

char binfile[]="bins";
const DOUBLE pimax=40.0;
double boxsize=420.0;
const int autocorr=1;
#ifdef USE_OMP
const int nthreads=4;
#endif
//end of global variables


int test_periodic_DD(void)
{
  struct timeval t0,t1; 
  assert(ND1==ND2 && x1==x2 && y1==y2 && z1==z2 && "Running test_periodic_DD() - the pointers should be identical");

  //Do the straight-up DD counts                                                                                                                                                                                                             
  gettimeofday(&t0,NULL);
  results_countpairs *results = countpairs(ND1,x1,y1,z1,
										   ND2,x2,y2,z2,
#ifdef USE_OMP
										   nthreads,
#endif
										   autocorr,
										   binfile);
  gettimeofday(&t1,NULL);
  double pair_time = ADD_DIFF_TIME(t0,t1);
  DOUBLE rlow=results->rupp[0];
  FILE *fp=NULL;
  
  fp=my_fopen(tmpoutputfile,"w");
  for(int i=1;i<results->nbin;i++) {
	fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->rpavg[i],rlow,results->rupp[i]);
	rlow=results->rupp[i];
  }
  fclose(fp);
	
  char execstring[MAXLEN];
  my_snprintf(execstring,MAXLEN,"diff -q Mr19_DD_periodic %s",tmpoutputfile);
  int ret=system(execstring);
  if(ret==EXIT_SUCCESS) {
	fprintf(stderr,ANSI_COLOR_GREEN "PASSED: 3-D periodic DD calculation. Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", pair_time);
  } else {
	fprintf(stderr,ANSI_COLOR_RED "FAILED: 3-D periodic DD calculation. Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", pair_time);
  }

  free_results(&results);
  return ret;
}

int test_periodic_DDrppi(void)
{
  struct timeval t0,t1; 
  assert(ND1==ND2 && x1==x2 && y1==y2 && z1==z2 && "Running test_periodic_DD() - the pointers should be identical");
  gettimeofday(&t0,NULL);
  
  results_countpairs_rp_pi *results = countpairs_rp_pi(ND1,x1,y1,z1,
													   ND2,x2,y2,z2,
#ifdef USE_OMP
													   nthreads,
#endif
													   autocorr,
													   binfile,
													   pimax);
  
  gettimeofday(&t1,NULL);
  double pair_time = ADD_DIFF_TIME(t0,t1);
  const int npibin = results->npibin;
  const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
  FILE *fp=my_fopen(tmpoutputfile,"w");
  for(int i=1;i<results->nbin;i++) {
	const double logrp = LOG10(results->rupp[i]);
	for(int j=0;j<npibin;j++) {
	  int index = i*(npibin+1) + j;
	  fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->rpavg[index],logrp,(j+1)*dpi);
	}
  }
  fclose(fp);
  char execstring[MAXLEN];
  my_snprintf(execstring,MAXLEN,"diff -q Mr19_DDrppi_periodic %s",tmpoutputfile);
  int ret=system(execstring);
  if(ret==EXIT_SUCCESS) {
	fprintf(stderr,ANSI_COLOR_GREEN "PASSED: periodic DDrppi calculation. Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", pair_time);
  } else {
	fprintf(stderr,ANSI_COLOR_RED "FAILED: periodic DDrppi calculation. Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", pair_time);
  }
	

  //free the result structure
  free_results_rp_pi(&results);
  return ret;
}

int test_wp(void)
{
  struct timeval t0,t1; 
  assert(ND1==ND2 && x1==x2 && y1==y2 && z1==z2 && "Running test_periodic_DD() - the pointers should be identical");

  gettimeofday(&t0,NULL);
  results_countpairs_wp *results = countpairs_wp(ND1,x1,y1,z1,
												 boxsize,
#ifdef USE_OMP
												 nthreads,
#endif
												 binfile,
												 pimax);
  gettimeofday(&t1,NULL);
  double pair_time = ADD_DIFF_TIME(t0,t1);
  DOUBLE rlow=results->rupp[0];
  FILE *fp=my_fopen(tmpoutputfile,"w");
  for(int i=1;i<results->nbin;++i) {
	fprintf(fp,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",results->wp[i],results->rpavg[i],rlow,results->rupp[i],results->npairs[i]);
	rlow=results->rupp[i];
  }
  fclose(fp);
  char execstring[MAXLEN];
  my_snprintf(execstring,MAXLEN,"diff -q Mr19_wp %s",tmpoutputfile);
  int ret=system(execstring);
  if(ret==EXIT_SUCCESS) {
	fprintf(stderr,ANSI_COLOR_GREEN "PASSED: periodic wp calculation. Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", pair_time);
  } else {
	fprintf(stderr,ANSI_COLOR_RED "FAILED: periodic wp calculation. Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", pair_time);
  }
  
  //free the result structure
  free_results_wp(&results);
  return ret;
}

int main(int argc, char **argv)
{
  struct timeval t0,t1;
  char file[]="data/gals_Mr19.ff";
  char fileformat[]="f";
  char binfile[]="bins";
  DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;
  struct timeval t0,t1; 
#ifdef USE_OMP
  const int nthreads=4;
#endif
  const int64_t ND1 = read_positions(file,fileformat,(void **) &x1,(void **) &y1,(void **) &z1,sizeof(DOUBLE));
  int autocorr=1;
  DOUBLE *x2 = x1;
  DOUBLE *y2 = y1;
  DOUBLE *z2 = z1;
  int64_t ND2 = ND1;

  gettimeofday(&t0,NULL);
  
  int failed=0;
  int status;
  
  const char allfunction_names[][MAXLEN] = {"tests_periodic_DD","tests_periodic_DDrppi","tests_wp"};
  const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);
  int (*allfunctions[]) (void) = {test_periodic_DD,test_periodic_DDrppi,test_wp};
  int total_tests=0;
  
  if(argc==1) {
	//nothing was passed at the command-line run all tests
	for(int i=0;i<ntests;i++) {
	  status = (*allfunctions[i])();
	  total_tests++;
	  if(status != EXIT_SUCCESS)  {
		failed++;
	  }
	}
  } else {
	//run specific tests
	for(int i=1;i<argc;i++){
	  int this_test_num = atoi(argv[i]);
	  if(this_test_num >= 0 && this_test_num < ntests) {
		total_tests++;
		status = (*allfunctions[this_test_num])();
		if(status != EXIT_SUCCESS)  {
		  failed++;
		}
	  }
	}
  }

  
  gettimeofday(&t1,NULL);
  double total_time = ADD_DIFF_TIME(t0,t1);
  if(failed > 0) {
	fprintf(stderr,ANSI_COLOR_RED "FAILED %d out of %d tests. Total time = %8.2lf seconds " ANSI_COLOR_RESET "\n", failed, total_tests, total_time);
  } else {
	fprintf(stderr,ANSI_COLOR_GREEN "PASSED: ALL %d tests. Total time = %8.2lf seconds " ANSI_COLOR_RESET "\n", total_tests, total_time);
	char execstring[MAXLEN];
	my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
	system(execstring);
  }
  free(x1);free(y1);free(z1);
  return EXIT_SUCCESS;
}
