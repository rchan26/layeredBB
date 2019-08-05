#ifndef BESSEL_LAYER_SIM
#define BESSEL_LAYER_SIM

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace Rcpp;

List bessel_layer_simulation(const double &x, const double &y,
                             const double &s, const double &t,
                             Rcpp::NumericVector &a);

#endif