#ifndef __COSMOLOGY_PARAMS_H
#define __COSMOLOGY_PARAMS_H
extern int cosmology_initialized;
extern double OMEGA_M;
extern double OMEGA_B;
extern double OMEGA_L;
extern double HUBBLE;
extern double LITTLE_H;
extern double SIGMA_8;
extern double NS;

void init_cosmology(const int lasdamas_cosmology);
#endif
