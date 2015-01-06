/*
	Example code to show how to use the correlation function libraries
	Author: Manodeep Sinha <manodeep@gmail.com>
	Date: At some point in early 2015. 



*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "function_precision.h"
#include "io.h"
#include "defs.h"
#include "utils.h"

#include "countpairs.h"
#include "countpairs_rp_pi.h"
#include "countpairs_wp.h"

#ifndef MAXLEN
#define MAXLEN 500
#endif

//Just to output some colors

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_BLUE    "\x1b[34m"

/* #define ANSI_COLOR_YELLOW  "\x1b[33m" */
/* #define ANSI_COLOR_MAGENTA "\x1b[35m" */
/* #define ANSI_COLOR_CYAN    "\x1b[36m" */



void Printhelp(void);

void Printhelp(void)
{
  fprintf(stderr,ANSI_COLOR_RED "=========================================================================\n") ;
  fprintf(stderr,"   --- run_correlations file format binfile boxsize [Nthreads] \n") ;
  fprintf(stderr,"   --- Measure the auto-correlation functions DD(r), DD(rp, pi) and wp(rp) for a single file\n");
  fprintf(stderr,"     * fille        = name of data file\n") ;
  fprintf(stderr,"     * format       = format of data file  (a=ascii, f=fast-food)\n") ;
  fprintf(stderr,"     * binfile      = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
  fprintf(stderr,"     * boxsize      = BoxSize (in same units as X/Y/Z of the data)\n");
#ifdef USE_OMP
  fprintf(stderr,"     * numthreads   = number of threads to use\n");
#endif
  fprintf(stderr,"=========================================================================" ANSI_COLOR_RESET "\n") ;	
}

int main(int argc, char **argv)
{
	char file[MAXLEN];
	char fileformat[MAXLEN];
	char binfile[MAXLEN];
	DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;
	double boxsize;
	struct timeval t0,t1;
		
#ifndef USE_OMP
	const char argnames[][30]={"file","format","binfile","boxsize"};
#else
	int nthreads=4;//default to 4 threads
	const char argnames[][30]={"file","format","binfile","boxsize","Nthreads"};
#endif
	int nargs=sizeof(argnames)/(sizeof(char)*30);
	
	if(argc > 1) {
		//Command-line options were supplied - check that they are correct
		if(argc < (nargs + 1) ) {
			//Not enough options were supplied
			Printhelp();
			exit(EXIT_FAILURE);
		} else {
			//Correct number of options - let's parse them.
			my_snprintf(file,MAXLEN, "%s",argv[1]);
			my_snprintf(fileformat,MAXLEN, "%s",argv[2]);
			my_snprintf(binfile,MAXLEN,"%s",argv[3]);
			boxsize=atof(argv[4]);
#ifdef USE_OMP
			nthreads = atoi(argv[5]);
#endif			
		}
	} else {
		my_snprintf(file, MAXLEN, "%s", "../tests/data/gals_Mr19.ff");
		my_snprintf(fileformat, MAXLEN, "%s","f");
		my_snprintf(binfile, MAXLEN,"%s","../tests/bins");
		boxsize=420.0;
	}

	fprintf(stderr,ANSI_COLOR_BLUE  "Running `%s' with the parameters \n",argv[0]);
	fprintf(stderr,"\n\t\t -------------------------------------\n");
	fprintf(stderr,"\t\t %-10s = %s \n",argnames[0],file);
	fprintf(stderr,"\t\t %-10s = %s \n",argnames[1],fileformat);
	fprintf(stderr,"\t\t %-10s = %s \n",argnames[2],binfile);
	fprintf(stderr,"\t\t %-10s = %lf\n",argnames[3],boxsize);
#ifdef USE_OMP
	fprintf(stderr,"\t\t %-10s = %d\n",argnames[4],nthreads);
#endif	
	fprintf(stderr,"\t\t -------------------------------------" ANSI_COLOR_RESET "\n");

	//Read-in the data
	const int64_t ND1 = read_positions(file,fileformat,(void **) &x1,(void **) &y1,(void **) &z1,sizeof(DOUBLE));

	int autocorr=1;
	DOUBLE *x2 = x1;
	DOUBLE *y2 = y1;
	DOUBLE *z2 = z1;
	int64_t ND2 = ND1;
	
	//Do the straight-up DD counts
	{
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
		fprintf(stderr,ANSI_COLOR_GREEN "Done 3-d auto-correlation. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", ND1, pair_time);
		//The results structure contains the pair-counts
		
		
		//free the result structure
		free_results(&results);
	}
	
	//Do the DD(rp, pi) counts
	{
		gettimeofday(&t0,NULL);
		DOUBLE pimax=40.0;
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
		fprintf(stderr,ANSI_COLOR_GREEN "Done DD(rp,pi) auto-correlation. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", ND1, pair_time);		
		

		//free the result structure
		free_results_rp_pi(&results);
	}


			
	//Do the wp counts
	{
		gettimeofday(&t0,NULL);
		DOUBLE pimax=40.0;
		results_countpairs_wp *results = countpairs_wp(ND1,x1,y1,z1,
																									 boxsize,
#ifdef USE_OMP
																											nthreads,
#endif
																									 binfile,
																									 pimax);
		gettimeofday(&t1,NULL);
		double pair_time = ADD_DIFF_TIME(t0,t1);
		fprintf(stderr,ANSI_COLOR_GREEN "Done wp. Ngalaxies = %12"PRId64" Time taken = %8.2lf seconds" ANSI_COLOR_RESET "\n", ND1, pair_time);

		//free the result structure
		free_results_wp(&results);
	}

			
	
	free(x1);free(y1);free(z1);
	return EXIT_SUCCESS;
}
