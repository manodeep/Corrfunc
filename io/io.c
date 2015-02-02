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
#include "function_precision.h"

#ifndef MEMORY_INCREASE_FAC
#define MEMORY_INCREASE_FAC 1.1
#endif

int64_t read_positions(const char *filename, const char *format, const size_t size, const int num_fields, ...)
{
  int64_t np;
  assert(num_fields >= 1 && "You have to request at least one field to read-in");
  assert(size == sizeof(DOUBLE) && "Requested size of an item does not match sizeof(DOUBLE)");

  DOUBLE *data[num_fields];
  if(strncmp(format,"f",1)==0) { /*Read-in fast-food file*/
	//read fast-food file
	size_t bytes=sizeof(int) + sizeof(float)*9 + sizeof(int);//skip fdat
	bytes += sizeof(int) + sizeof(float) + sizeof(int); //skip znow
	int idat[5];
	FILE *fp = my_fopen(filename,"r");
	my_ftread(idat,sizeof(int),5,fp);
	np = (int64_t) idat[1]; //idat[1] is int. 
	
	assert((size == 4 || size == 8) && "Size of each position element can be either 4 (float) or 8 (double)");
	
	for(int i=0;i<num_fields;i++) {
	  data[i] = my_malloc(size,np);
	}
	
	my_fseek(fp,bytes,SEEK_CUR);
	//Check that the file was written with the requested precision
	unsigned int dummy;
	my_fread(&dummy,sizeof(dummy), 1, fp);
	//so rewind by 4 bytes  prepare for calls to ftread
	my_fseek(fp, -sizeof(dummy), SEEK_CUR);
	dummy /= np;
	assert((dummy == 4 || dummy == 8) && "File must contain either 4 byte (float) or 8 byte(double) precision");
	
	if(dummy == size) {
	  for(int i=0;i<num_fields;i++) {
		my_ftread(data[i],size, np, fp);
	  }
	} else {
#ifndef SILENT
	  fprintf(stderr,"WARNING: File was written in a different precision than requested (file precision = %u requested precision = %zu)\n",dummy,size);
#endif
	  //Okay so the file was written in a different precision.
	  //First, print a warning message and then read-in correctly with the
	  //requested precision
	  if(dummy == 4) {
		assert(size == 8 && "Expected to be storing to doubles");
		float *tmp = my_malloc(dummy,np);
		//read-in the fields
		for(int i=0;i<num_fields;i++) {
		  my_ftread(tmp, dummy, np, fp);
		  double *tmp_pos = (double *) data[i];
		  for(int64_t j=0;j<np;j++) tmp_pos[j] = tmp[j];
		}  
				
		//free memory
		free(tmp);
	  } else {
		assert(size == 4 && "Expected to be storing to doubles");
		double *tmp = my_malloc(dummy,np);
				
		//read-in the fields
		for(int i=0;i<num_fields;i++) {
		  my_ftread(tmp, dummy, np, fp);
		  float *tmp_pos = (float *) data[i];
		  for(int64_t j=0;j<np;j++) tmp_pos[j] = tmp[j];
		}
		//free memory
		free(tmp);
	  }
	}

	fclose(fp);
  } else if(strncmp(format,"a",1)==0 || strncmp(format,"c",1)==0) { /* Read in ascii (white-space/comma) separated file*/
	int64_t i;
	int64_t nmax=300000;
	int nread;
	const int MAXBUFSIZE=10000;
	char buffer[MAXBUFSIZE];
	char delimiters[]=" ,\t";//delimiters are white-space, comma and tab

	for(i=0;i<num_fields;i++) {
	  data[i] = my_malloc(size,nmax);
	}

	FILE *fp = my_fopen(filename,"r");
	i = 0 ;
	while(fgets(buffer,MAXBUFSIZE,fp) != NULL) {
	  DOUBLE tmp;
	  char *token,*saveptr;
	  int flag = 1;
		char *copy=buffer;
	  for(int j=0;j<num_fields;j++,copy=NULL) {
		token = strtok_r(copy,delimiters,&saveptr);
		nread = sscanf(token,"%"DOUBLE_FORMAT,&tmp);
		if(nread == 1) {
		  (data[j])[i] = tmp;
		} else {
		  flag = 0;
		}
	  }
	  if(flag == 1) {
		i++;
	  }

	  if(i==nmax) {
		nmax *= MEMORY_INCREASE_FAC;
		while(nmax==i)
		  nmax += 5;
		
		for(int j=0;j<num_fields;j++) {
		  char varname[20];
		  snprintf(varname,20,"data[%d]",j);
		  data[j] = my_realloc(data[j],size,nmax,varname);
		}
	  }
	}
	np=i ;
	nmax=np;

	//release the extra memory.
	for(int j=0;j<num_fields;j++) {
	  char varname[20];
	  snprintf(varname,20,"data[%d]",j);
	  data[j] = my_realloc(data[j],size,nmax,varname);
	}
	fclose(fp);
  } else {
	fprintf(stderr,"ERROR: Unknown format `%s'\n",format);
	exit(EXIT_FAILURE);
  }

  va_list ap;
  va_start(ap,num_fields);

  for(int i=0;i<num_fields;i++) {
	DOUBLE **source = va_arg(ap, DOUBLE **);
	*source = data[i];
  }
  va_end(ap);

  return np;
}

