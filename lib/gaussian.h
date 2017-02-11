#ifndef GAUSSIAN_H
#define GAUSSIAN_H 1

#include "power.h"
#include "fft.h"

void gaussian_generate(const unsigned long seed,
		       PowerSpectrum* const ps,
		       const Float boxsize,
		       const bool fix_amplitude,
		       FFT* const fft_delta);

#endif
