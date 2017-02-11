//
// Using Boost program options
//   style:   ./options [options] <required arg>
//   example: ./options --x=3 filename
//

#include <iostream>
#include <cstdio>
#include <string>

#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "fs.h"

using namespace std;
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  comm_mpi_init(&argc, &argv);
    
  //
  // command-line options (Boost program_options)
  //
  options_description opt("mock_gaussian [options] <P(k) filename>");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "P(k) filename")
    ("nc", value<int>()->default_value(64), "number of grids per dim")
    ("boxsize", value<double>()->default_value(1000.0, "1000"),
     "length of the box on a side [1/h Mpc]")
    ("nbar", value<double>()->default_value(1.0e-4, "1.0e-4"),
     "output particle number density")
    ("scale-factor,a", value<double>()->default_value(1.0, "1"),
     "rescale power spectrum amplitude from a=1 to this scale factor")
    ("random-seed", value<size_t>()->default_value(1))
    ("fixed-amplitude", value<bool>()->default_value(false),
     "Fix the amplitude, only phase is random")
    ("mock-file,o", value<string>()->default_value("mock_gaussian.txt"),
     "output mock file name")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt; 
    return 0;
  }

  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  const double boxsize= vm["boxsize"].as<double>();
  const string filename= vm["filename"].as<string>();
  const size_t random_seed= vm["random-seed"].as<size_t>();
  const bool fix_amplitude= vm["fixed-amplitude"].as<bool>();
  
  PowerSpectrum* ps= new PowerSpectrum(filename.c_str());
  FFT* const fft_delta = new FFT("delta", nc, 0, false);
  
  gaussian_generate(random_seed, ps, boxsize, fix_amplitude, fft_delta);

  fft_delta->execute_inverse();

  Float* const delta= fft_delta->fx;

  //
  // Generate mock
  //
  string ofilename= vm["mock-file"].as<string>();
  FILE* const fp= fopen(ofilename.c_str(), "w");
  if(fp == 0) {
    msg_printf(msg_error, "Unable to write to file: %s\n", ofilename.c_str());
    comm_abort();
  }
  double sigma2= 0.0;
  const size_t ngrid= nc*nc*nc;
  
  

  const double nbar= vm["nbar"].as<double>();
  const double dx= boxsize/nc;
  const double nmean = nbar*dx*dx*dx;

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, random_seed + 1);
  //

  const double fac= 0.1; /// debug!!!

  for(int ix=0; ix<nc; ++ix) {
   double x0= ix*dx;
   for(int iy=0; iy<nc; ++iy) {
    double y0= iy*dx;
    for(int iz=0; iz<nc; ++iz) {
     double z0= iz*dx;
     size_t index= (ix*nc + iy)*nc + iz;
     double d= fac*delta[index];
     sigma2 += d*d;
     const double n= nmean*(1.0 + d);

     if(n > 0.0) {
       const int nout= gsl_ran_poisson(rng, n);
       for(int i=0; i<nout; ++i) {
	 double x = x0 + dx*gsl_rng_uniform(rng);
	 double y = y0 + dx*gsl_rng_uniform(rng);
	 double z = z0 + dx*gsl_rng_uniform(rng);

	 fprintf(fp, "%e %e %e\n", x, y, z);
       }
     }
     else
       msg_printf(msg_warn, "Warning: density is negative\n");
    }
   }
  }

  fclose(fp);
  msg_printf(msg_info, "mock written: %s\n", ofilename.c_str());

  msg_printf(msg_info, "sigma = %e\n", sqrt(sigma2/ngrid));
    

  comm_mpi_finalise();
  return 0;
}
