/* PROGRAM logbins

   --- logbins xmin xmax nbin
   --- Outputs the bins for a logarithmic bin spacing

      * xmin   = inner radius of smallest bin
      * xmax   = outer radius of largest bin
      * nbin  = number of bins (logarithmically spaced in theta)
      ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

int main(int argc, char *argv[])
{
  /*---Arguments-------------------------*/
  int nbin ;
  double xmin,xmax,logxmin,logxmax,dlogx;

  if(argc<4) {
    fprintf(stderr,"Usage `%s' rmin rmax nbins \n",argv[0]);
    exit(EXIT_FAILURE);
  } else {
    /*---Read-arguments-----------------------------------*/
    sscanf(argv[1],"%lf",&xmin) ;
    sscanf(argv[2],"%lf",&xmax) ;
    sscanf(argv[3],"%d",&nbin) ;
  }
  /*---Setup-Pairs-arrays-------------------------------*/
  logxmin = log10(xmin) ;
  logxmax = log10(xmax) ;
  dlogx = (logxmax-logxmin)/(double)nbin ;

  double left_edge=logxmin;
  for(int i=0;i<nbin;i++){
    double right_edge=left_edge+dlogx; 
    fprintf(stdout,"%20.12lf %20.12lf\n",pow(10.0,left_edge),pow(10.0,right_edge)) ;
    left_edge=right_edge;
  }
  exit(EXIT_SUCCESS);
}




