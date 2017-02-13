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
  // command-line options (Boost program_options)
  //
  options_description opt("gaussian_rsd [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("boxsize", value<double>()->default_value(1000.0),   "boxsize")
    ("nbar", value<double>()->default_value(0.01), "number density of points")
    ("random_seed", value<size_t>()->default_value(1))
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

  const size_t n= nbar*boxsize*boxsize*boxsize;

  const size_t random_seed= vm["random_seed"].as<size_t>();
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, random_seed);

  for(size_t i=0; i<n; ++i) {
    double x = boxsize*gsl_rng_uniform(rng);
    double y = boxsize*gsl_rng_uniform(rng);
    double z = boxsize*gsl_rng_uniform(rng);

    printf("%e %e %e\n", x, y, z);
  }

  gsl_rng_free(rng);
  
  return 0;
}

