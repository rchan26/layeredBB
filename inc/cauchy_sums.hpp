#ifndef CAUCHY_SUMS
#define CAUCHY_SUMS

#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// forward declaration for functions needed to calculate the Cauchy sums
double sigma_bar(const double &j, const double &x, const double &y,
                 const double &s, const double &t,
                 const double &l, const double &v);

double sigma(const double &j, const double &x, const double &y, 
             const double &s, const double &t,
             const double &l, const double &v);

double phi_bar(const double &j, const double &x, const double &y,
               const double &s, const double &t,
               const double &l, const double &v);

double phi(const double &j, const double &x, const double &y, 
           const double &s, const double &t,
           const double &l, const double &v);

double psi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v);

double chi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v);


// forward declarations for calculating the Cauchy sums (S_{2k+1}, S_{2k})
Rcpp::NumericVector calc_SgammaK_intervals(const int &k, const double &x, const double &y,
                                           const double &s, const double &t,
                                           const double &l, const double &v);

Rcpp::NumericVector calc_SdeltaK_1_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v);

Rcpp::NumericVector calc_SdeltaK_2_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v);

Rcpp::NumericVector calc_SdeltaK_intervals(const int &k, const double &x, const double &y,
                                           const double &s, const double &t,
                                           const double &min, const double &v);

#endif