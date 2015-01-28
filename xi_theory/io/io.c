/* File: io.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#include "io.h"
#include "ftread.h"
#include "utils.h"

#ifndef MEMORY_INCREASE_FAC
#define MEMORY_INCREASE_FAC 1.1
#endif

int64_t read_positions(const char *filename, const char *format, void **xpos, void **ypos, void **zpos, const size_t size)
{
  int64_t np;
  void *x,*y,*z;
  if(strncmp(format,"f",1)==0) { /*Read-in fast-food file*/
	//read fast-food file
	size_t bytes=sizeof(int) + sizeof(float)*9 + sizeof(int);//skip fdat
	bytes += sizeof(int) + sizeof(float) + sizeof(int); //skip znow
	int idat[5];
	FILE *fp = my_fopen(filename,"r");
	my_ftread(idat,sizeof(int),5,fp);
	np = (int64_t) idat[1]; //idat[1] is int
	
	assert((size == 4 || size == 8) && "Size of each position element can be either 4 (float) or 8 (double)");
	
	x = my_malloc(size,np);
	y = my_malloc(size,np);
	z = my_malloc(size,np);
	
	my_fseek(fp,bytes,SEEK_CUR);
	//Check that the file was written with the requested precision
	unsigned int dummy;
	my_fread(&dummy,sizeof(dummy), 1, fp);
	//so rewind by 4 bytes  prepare for calls to ftread
	my_fseek(fp, -sizeof(dummy), SEEK_CUR);
	dummy /= np;
	assert((dummy == 4 || dummy == 8) && "File must contain either 4 byte (float) or 8 byte(double) precision");
	
	if(dummy == size) {
	  my_ftread(x,size, np, fp);
	  my_ftread(y,size, np, fp);
	  my_ftread(z,size, np, fp);
	} else {
	  fprintf(stderr,"WARNING: File was written in a different precision than requested (file precision = %u requested precision = %zu)\n",dummy,size);
	  //Okay so the file was written in a different precision.
	  //First, print a warning message and then read-in correctly with the
	  //requested precision
	  if(dummy == 4) {
		assert(size == 8 && "Expected to be storing to doubles");
		float *tmp = my_malloc(dummy,np);
		double *tmp_x = (double *) x;
		double *tmp_y = (double *) y;
		double *tmp_z = (double *) z;
		//read-in x
		my_ftread(tmp, dummy, np, fp);
		for(int64_t i=0;i<np;i++) tmp_x[i] = tmp[i];
			  
		//read-in y
		my_ftread(tmp, dummy, np, fp);
		for(int64_t i=0;i<np;i++) tmp_y[i] = tmp[i];
				
		//read-in z
		my_ftread(tmp, dummy, np, fp);
		for(int64_t i=0;i<np;i++) tmp_z[i] = tmp[i];
				
		//free memory
		free(tmp);
	  } else {
		assert(size == 4 && "Expected to be storing to doubles");
		double *tmp = my_malloc(dummy,np);
		float *tmp_x = (float *) x;
		float *tmp_y = (float *) y;
		float *tmp_z = (float *) z;
				
				
		//read-in x
		my_ftread(tmp, dummy, np, fp);
		for(int64_t i=0;i<np;i++) tmp_x[i] = tmp[i];

		//read-in y
		my_ftread(tmp, dummy, np, fp);
		for(int64_t i=0;i<np;i++) tmp_y[i] = tmp[i];

		//read-in z
		my_ftread(tmp, dummy, np, fp);
		for(int64_t i=0;i<np;i++) tmp_z[i] = tmp[i];
				
		//free memory
		free(tmp);
	  }
	}

	fclose(fp);
  } else if(strncmp(format,"a",1)==0) { /* Read in ascii file*/
	int64_t i;
	int64_t nmax=300000;
	const int nitems=3;
	int nread;
	const int MAXBUFSIZE=10000;
	char buffer[MAXBUFSIZE];
	x = my_malloc(size,nmax);
	y = my_malloc(size,nmax);
	z = my_malloc(size,nmax);
	FILE *fp = my_fopen(filename,"r");

	i = 0 ;
	while(fgets(buffer,MAXBUFSIZE,fp) != NULL) {
	  if(size==8) {
		double dbl_x,dbl_y,dbl_z;
		nread = sscanf(buffer,"%lf %lf %lf",&dbl_x,&dbl_y,&dbl_z);
		if(nread == nitems) {
		  if(isfinite(dbl_x) && isfinite(dbl_y) && isfinite(dbl_z)){
			((double *) x)[i] = dbl_x;
			((double *) y)[i] = dbl_y;
			((double *) z)[i] = dbl_z;
			i++ ;
		  } else {
			fprintf(stderr,"WARNING: Nans found in data for i=%"PRId64" file = `%s'\n",i,filename);
		  }
		}
	  } else {
		float flt_x,flt_y,flt_z;
		nread = sscanf(buffer,"%f %f %f"   ,&flt_x,&flt_y,&flt_z);
		if(isfinite(flt_x) && isfinite(flt_y) && isfinite(flt_z)){
		  ((float *) x)[i] = flt_x;
		  ((float *) y)[i] = flt_y;
		  ((float *) z)[i] = flt_z;
		  i++ ;
		} else {
		  fprintf(stderr,"WARNING: Nans found in data for i=%"PRId64" file = `%s'\n",i,filename);
		}
	  }
	  if(i==nmax) {
		nmax *= MEMORY_INCREASE_FAC;
		while(nmax==i)
		  nmax += 5;
			
		x = my_realloc(x,size,nmax,"x");
		y = my_realloc(y,size,nmax,"y");
		z = my_realloc(z,size,nmax,"z");
	  }
	}
	np=i ;
	nmax=np;
		
	//release the extra memory
	x = my_realloc(x,size,nmax,"x");
	y = my_realloc(y,size,nmax,"y");
	z = my_realloc(z,size,nmax,"z");
	fclose(fp);
  } else if(strncmp(format,"c",1)==0) { /* Read in ascii csv file*/
	int64_t i;
	int64_t nmax=300000;
	const int nitems=3;
	int nread;
	const int MAXBUFSIZE=10000;
	char buffer[MAXBUFSIZE];
	x = my_malloc(size,nmax);
	y = my_malloc(size,nmax);
	z = my_malloc(size,nmax);
	FILE *fp = my_fopen(filename,"r");

	i = 0 ;
	while(fgets(buffer,MAXBUFSIZE,fp) != NULL) {
	  if(size==8) {
		double dbl_x,dbl_y,dbl_z;
		nread = sscanf(buffer,"%lf,%lf,%lf",&dbl_x,&dbl_y,&dbl_z);
		if(nread == nitems) {
		  if(isfinite(dbl_x) && isfinite(dbl_y) && isfinite(dbl_z)){
			((double *) x)[i] = dbl_x;
			((double *) y)[i] = dbl_y;
			((double *) z)[i] = dbl_z;
			i++ ;
		  } else {
			fprintf(stderr,"WARNING: Nans found in data for i=%"PRId64" file = `%s'\n",i,filename);
		  }
		}
	  } else {
		float flt_x,flt_y,flt_z;
		nread = sscanf(buffer,"%f,%f,%f"   ,&flt_x,&flt_y,&flt_z);
		if(isfinite(flt_x) && isfinite(flt_y) && isfinite(flt_z)){
		  ((float *) x)[i] = flt_x;
		  ((float *) y)[i] = flt_y;
		  ((float *) z)[i] = flt_z;
		  i++ ;
		} else {
		  fprintf(stderr,"WARNING: Nans found in data for i=%"PRId64" file = `%s'\n",i,filename);
		}
	  }
	  if(i==nmax) {
		nmax *= MEMORY_INCREASE_FAC;
		while(nmax==i)
		  nmax += 5;
			
		x = my_realloc(x,size,nmax,"x");
		y = my_realloc(y,size,nmax,"y");
		z = my_realloc(z,size,nmax,"z");
	  }
	}
	np=i ;
	nmax=np;
		
	//release the extra memory
	x = my_realloc(x,size,nmax,"x");
	y = my_realloc(y,size,nmax,"y");
	z = my_realloc(z,size,nmax,"z");
	fclose(fp);
  } else {
	fprintf(stderr,"ERROR: Unknown format `%s'\n",format);
	exit(EXIT_FAILURE);
  }

  *xpos = x;
  *ypos = y;
  *zpos = z;
  return np;
}

