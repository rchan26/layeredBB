#ifndef LAYERED_BROWNIAN_BRIDGE
#define LAYERED_BROWNIAN_BRIDGE

#include <RcppArmadillo.h>

using namespace Rcpp;

// forward declaration layared brownian bridge sampler
Rcpp::List layered_brownian_bridge(const double &x,
                                   const double &y,
                                   const double &s,
                                   const double &t,
                                   const Rcpp::List &bessel_layer,
                                   const Rcpp::NumericVector &times);

// forward declaration multi-dimensional layared brownian bridge sampler
Rcpp::List multi_layered_brownian_bridge(const int &dim,
                                         const arma::vec &x,
                                         const arma::vec &y,
                                         const double &s,
                                         const double &t,
                                         const Rcpp::List &bessel_layers,
                                         const Rcpp::NumericVector &times);

#endif