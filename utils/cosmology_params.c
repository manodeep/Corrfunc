/* File: cosmology_params.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#include "cosmology_params.h"

//Declare the variables
double OMEGA_M;
double OMEGA_B;
double OMEGA_L;
double HUBBLE;
double LITTLE_H;
double SIGMA_8;
double NS;
int cosmology=-1;

void init_cosmology(const int which_cosmology)
{

	switch(which_cosmology)
	{
	case 1: 
		//LasDamas Cosmology
		OMEGA_M=0.25;
		OMEGA_B=0.04;
		OMEGA_L=0.75;
		LITTLE_H=0.7;
		SIGMA_8=0.8;
		NS=1.0;
		break;
	case 2: 
    //Planck cosmology
		OMEGA_M=0.302;
		OMEGA_B=0.048;
		OMEGA_L=1.0-OMEGA_M;
		LITTLE_H=0.681;
		SIGMA_8=0.828;
		NS=0.96;
		break;

	default:
		fprintf(stderr,"ERROR: Cosmology=%d not implemented..exiting\n",which_cosmology);
		exit(EXIT_FAILURE);
	}
	
	HUBBLE=100.0*LITTLE_H;
	cosmology=which_cosmology;
}
