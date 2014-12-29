/* PROGRAM DDrppi

   --- DDrppi rpmin rpmax nrpbin data1 data2 [pimax] > wpfile
   --- Measure the cross-correlation function xi(rp,pi) for two different
       data files (or autocorrelation if data1=data2).

      * rpmin   = inner radius of smallest bin (in Mpc/h)
      * rpmax   = outer radius of largest bin
      * nrpbin  = number of bins (logarithmically spaced in r)
      * data1   = name of first data file
      * data2   = name of second data file
      *[pimax]  = maximum value of line-of-sight separation (default=40 Mpc/h)
      > DDfile  = name of output file <logrp log(<rp>) pi pairs>
      ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "defs.h" //for ADD_DIFF_TIME
#include "function_precision.h" //definition of DOUBLE
#include "countpairs_wp.h" //function proto-type for countpairs
#include "io.h" //function proto-type for file input
#include "utils.h" //general utilities

#include "sglib.h"

void Printhelp(void);

int main(int argc, char *argv[])
{

  /*---Arguments-------------------------*/
	double rpmin,rpmax;
  DOUBLE pimax ;
	char *file=NULL,*fileformat=NULL;
	double boxsize;
	
  /*---Data-variables--------------------*/
  int ND1=0,ND2=0;

  DOUBLE *x1=NULL,*y1=NULL,*z1=NULL;


  /*---Corrfunc-variables----------------*/
#ifndef USE_OMP
	const char argnames[][30]={"boxsize","file","format","binfile","pimax"};
#else
	int nthreads=2;
	const char argnames[][30]={"boxsize","file","format","binfile","pimax","Nthreads"};
#endif
  int nargs=sizeof(argnames)/(sizeof(char)*30);
  
  struct timeval t_end,t_start,t0,t1;
  double read_time=0.0;
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
    if(i < nargs) {
      fprintf(stderr,"\nMissing required parameters \n");
      for(i=argc;i<=nargs;i++)
				fprintf(stderr,"\t\t %s = `?'\n",argnames[i-1]);
    }
    return EXIT_FAILURE;
  }
	boxsize=atof(argv[1]);
  file=argv[2];
  fileformat=argv[3];

  /***********************
   *initializing the  bins
   ************************/
	double *rupp;
	int nbin;
	setup_bins(argv[4],&rpmin,&rpmax,&nbin,&rupp);
	assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
	assert(nbin > 0 && "Number of rp bins is valid");
	pimax=40.0;

#ifdef DOUBLE_PREC
	sscanf(argv[5],"%lf",&pimax) ;
#else    
	sscanf(argv[5],"%f",&pimax) ;
#endif    
		
		
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
  
  
  gettimeofday(&t0,NULL);
  /*---Read-data1-file----------------------------------*/
  ND1=read_positions(file,fileformat,(void **) &x1,(void **) &y1,(void **) &z1,sizeof(DOUBLE));
  gettimeofday(&t1,NULL);
  read_time += ADD_DIFF_TIME(t0,t1);

	//check that theee positions are within limits
  for(int i=0;i<ND1;i++) {
    assert(x1[i] >= 0.0 && x1[i] <= boxsize && "xpos is within limits [0, boxsize]");
		assert(y1[i] >= 0.0 && y1[i] <= boxsize && "ypos is within limits [0, boxsize]");
		assert(z1[i] >= 0.0 && z1[i] <= boxsize && "zpos is within limits [0, boxsize]");
  }

	//Sort the arrays on z
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,x1,i,j); \
	SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,y1,i,j); \
	SGLIB_ARRAY_ELEMENTS_EXCHANGER(DOUBLE,z1,i,j) }

	SGLIB_ARRAY_QUICK_SORT(DOUBLE, z1, ND1, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
	
  /*---Count-pairs--------------------------------------*/
	uint64_t *npairs = my_calloc(sizeof(*npairs),nbin);
#ifdef OUTPUT_RPAVG
	DOUBLE *rpavg = my_calloc(sizeof(*rpavg),nbin);
#endif	
	
  gettimeofday(&t0,NULL);
	countpairs_wp(ND1, x1, y1, z1,
								boxsize, pimax,
#ifdef USE_OMP
								nthreads,
#endif
								npairs,
#ifdef OUTPUT_RPAVG
								rpavg,
#endif
								rupp, nbin);
	
	gettimeofday(&t1,NULL);
  double pair_time = ADD_DIFF_TIME(t0,t1);

	//Output the results
	const DOUBLE avgweight2 = 1.0, avgweight1 = 1.0;
  const DOUBLE density=0.5*avgweight2*ND1/(boxsize*boxsize*boxsize);//pairs are not double-counted

	DOUBLE rlow=0.0 ;
	DOUBLE prefac_density_DD=avgweight1*ND1*density;
	DOUBLE twice_pimax = 2.0*pimax;
	DOUBLE xi_full[nbin];
	
	for (int kbin=0;kbin<nbin;kbin++)  {      /* loop over radial bins */
		const DOUBLE weight0 = (DOUBLE) npairs[kbin];

		/* compute xi, dividing summed weight by that expected for a random set */
		const DOUBLE vol=M_PI*(rupp[kbin]*rupp[kbin]-rlow*rlow)*twice_pimax;
		const DOUBLE weightrandom = prefac_density_DD*vol;
		xi_full[kbin] = (weight0/weightrandom-1)*twice_pimax;
		rlow=rupp[kbin] ;
	}                                     /* next radial bin */

	/* Note: we discard the first bin, to mimic the fact that close pairs
	 * are disregarded in SDSS data.
	 */
	rlow=rupp[0];
	for(int i=1;i<nbin;++i) {
#ifdef OUTPUT_RPAVG
		fprintf(stdout,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",xi_full[i],rpavg[i],rlow,rupp[i],npairs[i]);
#else		
		fprintf(stdout,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",xi_full[i],0.0,rlow,rupp[i],npairs[i]);
#endif		
		rlow=rupp[i];
	}
	


	free(x1);free(y1);free(z1);
	free(rupp);

	free(npairs);
#ifdef OUTPUT_RPAVG
	free(rpavg);
#endif	
	
  gettimeofday(&t_end,NULL);
  fprintf(stderr,"xi_rp_pi> Done -  ND1=%d ND2=%d. Time taken = %6.2lf seconds. read-in time = %6.2lf seconds pair-counting time = %6.2lf sec\n",
	  ND1,ND2,ADD_DIFF_TIME(t_start,t_end),read_time,pair_time);
  return EXIT_SUCCESS;
}

/*---Print-help-information---------------------------*/
void Printhelp(void)
{
  fprintf(stderr,"=========================================================================\n") ;
  fprintf(stderr,"   --- wp boxsize file format binfile pimax [Nthreads] > wpfile\n") ;
  fprintf(stderr,"   --- Measure the projected auto-correlation function wp(rp) for a single periodic box\n") ;
  fprintf(stderr,"       data files (or autocorrelation if data1=data2).\n") ;
	fprintf(stderr,"     * boxsize      = BoxSize (in same units as X/Y/Z of the data)\n") ;
  fprintf(stderr,"     * file         = name of first data file\n") ;
  fprintf(stderr,"     * format       = format of first data file  (a=ascii, f=fast-food)\n") ;
  fprintf(stderr,"     * binfile       = name of ascii file containing the r-bins (rmin rmax for each bin)\n") ;
  fprintf(stderr,"     * pimax         = maximum line-of-sight-separation\n") ;
#ifdef USE_OMP
	fprintf(stderr,"     * numthreads    = number of threads to use\n");
#endif

#ifdef OUTPUT_RPAVG	
  fprintf(stderr,"     > wpfile        = name of output file. Contains <wp  rpavg  rmin rmax npairs>\n") ;
#else
	fprintf(stderr,"     > wpfile        = name of output file. Contains <wp [rpavg=0.0] rmin rmax npairs>\n") ;
#endif	
  fprintf(stderr,"=========================================================================\n") ;
}



