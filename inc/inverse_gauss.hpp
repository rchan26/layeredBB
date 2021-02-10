#ifndef INVERSE_GAUSS
#define INVERSE_GAUSS

#include <Rcpp.h>
#include <random>
#include <cmath>

// forward declaration for simulation of Inverse Gaussian variable
double inv_gauss_sampler(const double &mu, const double &lambda);

#endif