#ifndef BROWNIAN_BRIDGE
#define BROWNIAN_BRIDGE

#include <Rcpp.h>
#include <random>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// forward declaration for Brownian motion path samplers
Rcpp::NumericMatrix Brownian_motion_path_sampler(const double &x,
                                                 const Rcpp::NumericVector &times);

// forward declaration for Brownian bridge path samplers
Rcpp::List Brownian_bridge_path_sampler(const double &x,
                                        const double &y,
                                        const double &s,
                                        const double &t,
                                        const Rcpp::NumericVector &times);

// forward declaration for multi-dimensional Brownian bridge path samplers
Rcpp::List multi_brownian_bridge(const int &dim,
                                 const Rcpp::NumericVector &x,
                                 const Rcpp::NumericVector &y,
                                 const double &s,
                                 const double &t,
                                 const Rcpp::NumericVector &times);

// forward declaration for Brownian Bridge minimum point sampler
Rcpp::NumericVector min_sampler(const double &x,
                                const double &y,
                                const double &s,
                                const double &t,
                                const double &low_bound,
                                const double &up_bound,
                                const bool &checks);

// forward declaration for minimum Bessel Bridge simulation at time q
double min_Bessel_bridge_sampler(const double &x,
                                 const double &y,
                                 const double &s,
                                 const double &t,
                                 const double &m,
                                 const double &tau,
                                 const double &q,
                                 const bool &checks);

// forward declaration for minimum Bessel Bridge path simulation
Rcpp::List min_Bessel_bridge_path_sampler(const double &x,
                                          const double &y,
                                          const double &s,
                                          const double &t,
                                          const double &m,
                                          const double &tau,
                                          const Rcpp::NumericVector &times,
                                          const bool &checks);

// forward declaration for maximum Bessel Bridge simulation at time q
Rcpp::NumericVector max_sampler(const double &x,
                                const double &y,
                                const double &s,
                                const double &t,
                                const double &low_bound,
                                const double &up_bound,
                                const bool &checks);

// forward declaration for maximum Bessel Bridge simulation at time q
double max_Bessel_bridge_sampler(const double &x,
                                 const double &y,
                                 const double &s,
                                 const double &t,
                                 const double &m,
                                 const double &tau,
                                 const double &q,
                                 const bool &checks);

// forward declaration for maximum Bessel Bridge path simulation
Rcpp::List max_Bessel_bridge_path_sampler(const double &x,
                                          const double &y,
                                          const double &s,
                                          const double &t,
                                          const double &m,
                                          const double &tau,
                                          const Rcpp::NumericVector &times,
                                          const bool &checks);

#endif
