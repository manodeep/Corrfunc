#include "cosmology_params.h"

//Declare the variables
double OMEGA_M;
double OMEGA_B;
double OMEGA_L;
double HUBBLE;
double LITTLE_H;
double SIGMA_8;
double NS;
int cosmology_initialized=0;

void init_cosmology(const int lasdamas_cosmology)
{	
	if (lasdamas_cosmology == 1) {
		//LasDamas Cosmology
		OMEGA_M=0.25;
		OMEGA_B=0.04;
		OMEGA_L=0.75;
		LITTLE_H=0.7;
		SIGMA_8=0.8;
		NS=1.0;
	}	else {
    //Planck cosmology
		OMEGA_M=0.302;
		OMEGA_B=0.048;
		OMEGA_L=1.0-OMEGA_M;
		LITTLE_H=0.681;
		SIGMA_8=0.828;
		NS=0.96;
	}
	
	HUBBLE=100.0*LITTLE_H;
	cosmology_initialized=1;
}
