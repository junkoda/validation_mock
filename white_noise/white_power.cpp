//
// Compute power spectrum of white noise
// including the mass assignment smoothing and aliasing
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  //
  // Command-line options (Boost program_options)
  //
  options_description opt("gaussian_rsd [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("boxsize", value<double>()->default_value(1000.0),   "boxsize")
    ("nc", value<size_t>()->default_value(64), "number of grids per dimension")
    ("nbar", value<double>()->default_value(2.5e-4), "number density")
    ;
  
  positional_options_description p;
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help")) {
    cout << opt; 
    return 0;
  }

  const double nbar= vm["nbar"].as<double>();
  const double boxsize= vm["boxsize"].as<double>();
  const size_t nc= vm["nc"].as<size_t>();

  const double knq= (M_PI*nc)/boxsize;
  const double knq_inv= boxsize/(M_PI*nc);


  //
  const double dk_fundamental= (M_PI)/boxsize;
  const double k_min= 0.5*dk_fundamental;
  const double k_max= knq;
  const double dk= 2.0*dk_fundamental;
  const int nbin= (int) round((k_max - k_min)/dk);
  
  //
  // Allocate memory
  //
  int* const nmodes= (int*) calloc(sizeof(int), nbin);
  double* const kmean= (double*) calloc(sizeof(double), nbin);
  double* const P= (double*) calloc(sizeof(double), nbin);
  
  const double pk_fac= (boxsize*boxsize*boxsize)/pow((double)nc, 6);
  const double sin_fac= 0.5*boxsize/nc;

  const double fac= 2.0*M_PI/boxsize;
  const size_t ncz= nc/2+1;



  assert(nbar > 0.0);
  const double nbar_inv= 1.0/nbar;

  for(int ix=0; ix<nc; ++ix) {
   double kx= ix <= nc/2 ? fac*ix : fac*(ix-nc);
   double sintx= sin(sin_fac*kx);
   double c1x= 1.0 - 2.0/3.0*sintx*sintx;

   for(int iy=0; iy<nc; ++iy) {
    double ky= iy <= nc/2 ? fac*iy : fac*(iy-nc);
    double sinty= sin(sin_fac*ky);
    double c1y= 1.0 - 2.0/3.0*sinty*sinty;

    int iz0 = !(kx > 0.0f || (kx == 0.0 && ky > 0.0));

    // Avoid double counting on kz=plain
    // k=(0,0,0) dropped because this is 0
    // iz0= 0 if kx>0 or (kx == 0 and ky > 0)
    //      1 otherwize
    //
      
    for(int iz=iz0; iz<nc/2+1; ++iz) {
      double kz= fac*iz;
      double sintz= sin(sin_fac*kz);
      double c1z= 1.0 - 2.0/3.0*sintz*sintz;
      
      double k= sqrt(kx*kx + ky*ky + kz*kz);
      double shot_noise= c1x*c1y*c1z*nbar_inv; // C1 function in Jing 2005

      int i= (int) floor((k - k_min)/dk);
      
      if(0 <= i && i < nbin) {
	nmodes[i]++;
	kmean[i] += k;
	P[i] += shot_noise;
      }
    }  
   }
  }

  for(int i=0; i<nbin; ++i) {
    if(nmodes[i] > 0) {
      printf("%e %e %d %e\n",
	     k_min + (i + 0.5)*dk,
	     kmean[i]/nmodes[i],
	     nmodes[i],
	     P[i]/nmodes[i]);
    }
  }
  
  return 0;
}

