//
// Generates random Gaussian field
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <gsl/gsl_rng.h>

#include "msg.h"
#include "mem.h"
#include "config.h"
#include "cosmology.h"
#include "power.h"
#include "fft.h"
#include "lpt.h"

using namespace std;

static unsigned int*
generate_seedtable(const int nc, gsl_rng* random_generator);

void gaussian_generate(const unsigned long seed,
		       PowerSpectrum* const ps,
		       const Float boxsize,
		       const bool fix_amplitude,
		       FFT* const fft_delta)

{
  const size_t nc= fft_delta->nc; assert(nc > 0);
  const size_t local_nx= fft_delta->local_nx;
  const size_t local_ix0= fft_delta->local_ix0;


  // Generates random Gaussian Field delta_k
  // from N-GenIC by Volker Springel
  msg_printf(msg_verbose, "Generating delta_k...\n");
  msg_printf(msg_info, "Random Seed = %lu\n", seed);
  if(fix_amplitude)
    msg_printf(msg_info, "Amplitude fixed");

  complex_t* delta_k= fft_delta->fk;
  
  const size_t nckz= nc/2 + 1;
  const double dk= 2.0*M_PI/boxsize;
  const double knq= nc*M_PI/boxsize; // Nyquist frequency
  const double fac= pow(2*M_PI/boxsize, 1.5);
  const double fac_2pi3= 1.0/(8.0*M_PI*M_PI*M_PI);
  
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);
  unsigned int* const seedtable=
    generate_seedtable(nc, random_generator);

  // clean the delta_k grid
  for(size_t ix=0; ix<local_nx; ix++)
   for(size_t iy=0; iy<nc; iy++)
    for(size_t iz=0; iz<nckz; iz++)
      for(int i=0; i<3; i++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	delta_k[index][0] = 0;
	delta_k[index][1] = 0;
      }

  double kvec[3];
  for(size_t ix=0; ix<nc; ix++) {
    size_t iix = nc - ix;
    if(iix == nc)
      iix = 0;

    if(!((local_ix0 <= ix  && ix  < (local_ix0 + local_nx)) ||
	 (local_ix0 <= iix && iix < (local_ix0 + local_nx))))
      continue;
    
    for(size_t iy=0; iy<nc; iy++) {
      gsl_rng_set(random_generator, seedtable[ix*nc + iy]);
      
      for(size_t iz=0; iz<nc/2; iz++) {
	double phase= gsl_rng_uniform(random_generator)*2*M_PI;
	double ampl;

	do
	  ampl = gsl_rng_uniform(random_generator);
	while(ampl == 0.0);

	if(ix == nc/2 || iy == nc/2 || iz == nc/2)
	  continue;
	if(ix == 0 && iy == 0 && iz == 0)
	  continue;
	
	if(ix < nc/2)
	  kvec[0]= dk*ix;
	else
	  kvec[0]= -dk*(nc - ix);
	
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);
	
	double kmag2 = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
	double kmag = sqrt(kmag2);
	
#ifdef SPHEREMODE
	// select a sphere in k-space
	if(kmag > knq)
	  continue;
#else
	if(fabs(kvec[0]) > knq)
	  continue;
	if(fabs(kvec[1]) > knq)
	  continue;
	if(fabs(kvec[2]) > knq)
	  continue;
#endif

	double delta2;
	if(fix_amplitude)
	  delta2= fac_2pi3*ps->P(kmag);
	else
	  delta2= -log(ampl)*fac_2pi3*ps->P(kmag);
	
	double delta_k_mag= fac*sqrt(delta2);
	// delta_k_mag -- |delta_k| extrapolated to a=1
	// Displacement is extrapolated to a=1

	if(iz > 0) {
	  if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
	    size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
	    delta_k[index][0]= delta_k_mag*cos(phase);
	    delta_k[index][1]= delta_k_mag*sin(phase);
	    /*
	    for(int i=0; i<3; i++) {
	      psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
	      psi_k[i][index][1]=  kvec[i]/kmag2*delta_k_mag*cos(phase);
	    }
	    */
	  }
	}
	else { // k=0 plane needs special treatment
	  if(ix == 0) {
	    if(iy >= nc/2)
	      continue;
	    else {
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t iiy= nc - iy; // note: j!=0 surely holds at this point
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;				size_t iindex= ((ix - local_ix0)*nc + iiy)*nckz + iz;

		delta_k[index][0]= delta_k_mag*cos(phase);
		delta_k[index][1]= delta_k_mag*sin(phase);

		// complex conjugate delta_k(-k) = delta_k(k)*
		delta_k[iindex][0]= delta_k_mag*cos(phase);
		delta_k[iindex][1]= -delta_k_mag*sin(phase);
	      }
	    }
	  }
	  else { // here comes i!=0 : conjugate can be on other processor!
	    if(ix >= nc/2)
	      continue;
	    else {
	      iix = nc - ix;
	      if(iix == nc)
		iix = 0;
	      int iiy = nc - iy;
	      if(iiy == nc)
		iiy = 0;
	      
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
		delta_k[index][0]= delta_k_mag*cos(phase);
		delta_k[index][1]= delta_k_mag*sin(phase);
	      }
	      
	      if(local_ix0 <= iix && iix < (local_ix0 + local_nx)) {
		size_t index= ((iix - local_ix0)*nc + iiy)*nckz + iz;
		delta_k[index][0]= delta_k_mag*cos(phase);
		delta_k[index][1]= -delta_k_mag*sin(phase);
	      }
	    }
	  }
	}
      }
    }
  }

  fft_delta->mode= fft_mode_k;
  
  gsl_rng_free(random_generator);
  
  free(seedtable);
}

//void gaussian_init(const int nc_, const double boxsize_, Mem* mem)
//{
  // nc_: number of particles per dimension
  // boxsize_: box length on a side
  // Mem: Memory object for LPT, can be 0.
  //      If mem = 0, memory is allocated exclusively for LPT
  /*
  if(mem)
    mem->use_from_zero(0);
  
  for(int i=0; i<3; i++)
    fft_psi[i]= new FFT("Psi_i", nc, mem, 0);

  for(int i=0; i<6; i++)
    fft_psi_ij[i]= new FFT("Psi_ij", nc, mem, 0);

  for(int i=0; i<3; i++)
    fft_psi2[i]= fft_psi_ij[i];

  fft_div_psi2= fft_psi_ij[3];
  
  */
//}

unsigned int* generate_seedtable(const int nc, gsl_rng* random_generator)
{
  unsigned int* const stable=
    (unsigned int *) malloc(nc*nc*sizeof(unsigned int));
  assert(stable);

  // from N-GenIC

  for(int i=0; i<nc/2; i++) {
    for(int j=0; j<i; j++)
      stable[i*nc + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[j*nc + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[(nc - 1 - i)*nc + j] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[(nc - 1 - j)*nc + i] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[i*nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[j*nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[(nc - 1 - i)*nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[(nc - 1 - j)*nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);
  }

  return stable;
}
