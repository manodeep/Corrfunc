#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities

void Printhelp(void);

int main(int argc, char *argv[])
{

  /*---Data-variables--------------------*/
  int ND1,ND2 ;

  DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;
  DOUBLE *x2=NULL,*y2=NULL,*z2=NULL;

  char *file1=NULL,*file2=NULL;
  char *fileformat1=NULL,*fileformat2=NULL;
	char *binfile=NULL;


  /*---Corrfunc-variables----------------*/
#ifndef USE_OMP
  const char argnames[][30]={"file1","format1","file2","format2","binfile"};
#else
  int nthreads=2;
  const char argnames[][30]={"file1","format1","file2","format2","binfile","Nthreads"};
#endif
  int nargs=sizeof(argnames)/(sizeof(char)*30);
  
  struct timeval t_end,t_start,t0,t1;
  double read_time=0.0;

#ifdef USE_ISPC
#error "ISPC features has not been implemented correctly yet. Please disable this Makefile option"
#endif


  gettimeofday(&t_start,NULL);
  
  /*---Read-arguments-----------------------------------*/
  if(argc< (nargs+1)) {
    Printhelp() ;
    fprintf(stderr,"\nFound: %d parameters\n ",argc-1);
		int i;
    for(i=1;i<argc;i++) {
      if(i <= nargs)
				fprintf(stderr,"\t\t %s = `%s' \n",argnames[i-1],argv[i]);
      else
				fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
    }
    if(i <= nargs) {
      fprintf(stderr,"\nMissing required parameters \n");
      for(i=argc;i<=nargs;i++)
				fprintf(stderr,"\t\t %s = `?'\n",argnames[i-1]);
    }
    return EXIT_FAILURE;
  }

  file1=argv[1];
  fileformat1=argv[2];
  file2=argv[3];
  fileformat2=argv[4];
	binfile=argv[5];
	
#ifdef USE_OMP
  nthreads=atoi(argv[6]);
  assert(nthreads >= 1 && "Number of threads must be at least 1");
#endif
  
  
  fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
  fprintf(stderr,"\n\t\t -------------------------------------\n");
  for(int i=1;i<argc;i++) {
    if(i <= nargs) {
      fprintf(stderr,"\t\t %-10s = %s \n",argnames[i-1],argv[i]);
    }  else {
      fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
    }
  }
  fprintf(stderr,"\t\t -------------------------------------\n");
  
  

  /*---Read-data1-file----------------------------------*/
  gettimeofday(&t0,NULL);
  ND1=read_positions(file1,fileformat1,(void **) &x1,(void **) &y1,(void **) &z1,sizeof(DOUBLE));
  gettimeofday(&t1,NULL);
  read_time += ADD_DIFF_TIME(t0,t1);

  int autocorr=0;
  if( strcmp(file1,file2)==0) {
    autocorr=1;
  }
  
  gettimeofday(&t0,NULL);  
  if (autocorr==0) {
    /*---Read-data2-file----------------------------------*/
		ND2=read_positions(file2,fileformat2,(void **) &x2,(void **) &y2,(void **) &z2,sizeof(DOUBLE));
    gettimeofday(&t1,NULL);
    read_time += ADD_DIFF_TIME(t0,t1);
  } else {
    ND2 = ND1;
    x2 = x1;
    y2 = y1;
    z2 = z1;
  }
  
  /*---Count-pairs--------------------------------------*/
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
  free(x1);free(y1);free(z1);
  if(autocorr == 0) {
	  free(x2);free(y2);free(z2);
  }

	double rlow=results->rupp[0];
  for(int i=1;i<results->nbin;i++) {
#ifdef OUTPUT_RPAVG	  
    fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->rpavg[i],rlow,results->rupp[i]);
#else
		fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],0.0,    rlow,results->rupp[i]);
#endif	
    rlow=results->rupp[i];
  }

	//free the memory in the results structx
	free_results(&results);
  gettimeofday(&t_end,NULL);
  fprintf(stderr,"xi> Done -  ND1=%d ND2=%d. Time taken = %6.2lf seconds. read-in time = %6.2lf seconds sec pair-counting time = %6.2lf sec\n",
		  ND1,ND2,ADD_DIFF_TIME(t_start,t_end),read_time,pair_time);
  return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{
  fprintf(stderr,"=========================================================================\n") ;
  fprintf(stderr,"   --- DDrppi data1 format1 data2 format2 binfile > DDfile\n") ;
  fprintf(stderr,"   --- Measure the cross-correlation function xi(rp,pi) for two different\n") ;
  fprintf(stderr,"       data files (or autocorrelation if data1=data2).\n") ;
  fprintf(stderr,"     * data1         = name of first data file\n") ;
  fprintf(stderr,"     * format1       = format of first data file\n") ;
  fprintf(stderr,"     * data2         = name of second data file\n") ;
  fprintf(stderr,"     * format2       = format of second data file\n") ;
  fprintf(stderr,"     * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
#ifdef USE_OMP
  fprintf(stderr,"     * numthreads    = number of threads to use\n");
#endif
  fprintf(stderr,"     > DDfile        = name of output file <npairs logrp pi pairs>\n") ;
	fprintf(stderr,"\n\tCompile options: \n");
#ifdef PERIODIC
	fprintf(stderr,"Periodic = True\n");
#else
	fprintf(stderr,"Periodic = False\n");
#endif

#ifdef OUTPUT_RPAVG
	fprintf(stderr,"Output RPAVG = True\n");
#else
	fprintf(stderr,"Output RPAVG = False\n");
#endif	

#ifdef DOUBLE_PREC
	fprintf(stderr,"Precision = double\n");
#else
	fprintf(stderr,"Precision = float\n");
#endif
	
#ifdef USE_AVX
	fprintf(stderr,"Use AVX = True\n");
#else	
	fprintf(stderr,"Use AVX = False\n");
#endif

#ifdef USE_OMP
	fprintf(stderr,"Use OMP = True\n");
#else	
	fprintf(stderr,"Use OMP = False\n");
#endif

	fprintf(stderr,"=========================================================================\n") ;
	
}




