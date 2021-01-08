#ifndef LAYERED_BROWNIAN_BRIDGE
#define LAYERED_BROWNIAN_BRIDGE

#include <RcppArmadillo.h>

using namespace Rcpp;

// forward declaration for function that finds the maximum value in a vector
double find_max(const Rcpp::NumericVector vect);

// forward declaration for function that finds the minimum value in a vector
double find_min(const Rcpp::NumericVector vect);

// forward declaration layared brownian bridge sampler
Rcpp::NumericMatrix layered_brownian_bridge(const double &x, 
                                            const double &y,
                                            const double &s, 
                                            const double &t,
                                            const Rcpp::List &bessel_layer,
                                            const Rcpp::NumericVector &times,
                                            const bool &remove_m);

// forward declaration multi-dimensional layared brownian bridge sampler
Rcpp::NumericMatrix multi_layered_brownian_bridge(const int &dim,
                                                  const arma::vec &x,
                                                  const arma::vec &y,
                                                  const double &s,
                                                  const double &t,
                                                  const Rcpp::List &bessel_layers,
                                                  Rcpp::NumericVector times);

#endif