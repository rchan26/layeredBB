#ifndef BROWNIAN_BRIDGE_MIN
#define BROWNIAN_BRIDGE_MIN

#include <Rcpp.h>
#include <random>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// forward declaration for Brownian bridge path samplers
Rcpp::NumericMatrix Brownian_bridge_path_sampler(const double &x, 
                                                 const double &y,
                                                 const double &s, 
                                                 const double &t,
                                                 Rcpp::NumericVector times);

// forward declaration for multi-dimensional Brownian bridge path samplers
Rcpp::NumericMatrix multi_brownian_bridge(const double &dim,
                                          const Rcpp::NumericVector &x,
                                          const Rcpp::NumericVector &y,
                                          const double &s,
                                          const double &t,
                                          const Rcpp::List &layers,
                                          Rcpp::NumericVector times);

// forward declaration for M_function that's used in min_sampler
double M_function(const double &a,
                  const double &x, 
                  const double &y, 
                  const double &s, 
                  const double &t);

// forward declaration for Brownian Bridge minimum point sampler
Rcpp::NumericVector min_sampler(const double &x,
                                const double &y,
                                const double &s,
                                const double &t,
                                const double &low_bound, 
                                const double &up_bound);

// forward declaration for minimum Bessel Bridge simulation at time q
double min_Bessel_bridge_sampler(const double &x, 
                                 const double &y,
                                 const double &s, 
                                 const double &t,
                                 const double &min, 
                                 const double &tau,
                                 const double &q);

// forward declaration for minimum Bessel Bridge path simulation
Rcpp::NumericMatrix min_Bessel_bridge_path_sampler(const double &x,
                                                   const double &y,
                                                   const double &s, 
                                                   const double &t,
                                                   const double &min,
                                                   const double &tau,
                                                   Rcpp::NumericVector times);

// forward declaration for maximum Bessel Bridge simulation at time q
Rcpp::NumericVector max_sampler(const double &x, 
                                const double &y,
                                const double &s, 
                                const double &t,
                                const double &low_bound, 
                                const double &up_bound);

// forward declaration for maximum Bessel Bridge simulation at time q
double max_Bessel_bridge_sampler(const double &x, 
                                 const double &y,
                                 const double &s, 
                                 const double &t,
                                 const double &max, 
                                 const double &tau,
                                 const double &q);

// forward declaration for maximum Bessel Bridge path simulation
Rcpp::NumericMatrix max_Bessel_bridge_path_sampler(const double &x,
                                                   const double &y,
                                                   const double &s, 
                                                   const double &t,
                                                   const double &max,
                                                   const double &tau,
                                                   Rcpp::NumericVector times);

#endif
