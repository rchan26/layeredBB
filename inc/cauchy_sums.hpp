#ifndef CAUCHY_SUMS
#define CAUCHY_SUMS

#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// forward declaration for functions needed to calculate the Cauchy sums
double easigma_bar(const double &j,
                   const double &x,
                   const double &y,
                   const double &s,
                   const double &t,
                   const double &l,
                   const double &v);

double easigma(const double &j,
               const double &x,
               const double &y,
               const double &s,
               const double &t,
               const double &l,
               const double &v);

double eaphi_bar(const double &j,
                 const double &x,
                 const double &y,
                 const double &s,
                 const double &t,
                 const double &l,
                 const double &v);

double eaphi(const double &j,
             const double &x,
             const double &y,
             const double &s,
             const double &t,
             const double &l,
             const double &v);

double eapsi(const double &j,
             const double &xoy,
             const double &s,
             const double &t,
             const double &min,
             const double &v);

double eachi(const double &j,
             const double &xoy,
             const double &s,
             const double &t,
             const double &min,
             const double &v);

double eagamma(const int &n,
               const double &x,
               const double &y,
               const double &s,
               const double &t,
               const double &l,
               const double &v);

double eadelta1(const int &n,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &min,
                const double &v);

double eadelta2(const int &n,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &min,
                const double &v);

double eadelta(const int &n,
               const double &x,
               const double &y,
               const double &s,
               const double &t,
               const double &min,
               const double &v);

Rcpp::NumericVector eagamma_intervals(const int &k,
                                      const double &x,
                                      const double &y,
                                      const double &s,
                                      const double &t,
                                      const double &l,
                                      const double &v);

Rcpp::NumericVector eadelta1_intervals(const int &k,
                                       const double &x,
                                       const double &y,
                                       const double &s,
                                       const double &t,
                                       const double &min,
                                       const double &v);

Rcpp::NumericVector eadelta2_intervals(const int &k,
                                       const double &x,
                                       const double &y,
                                       const double &s,
                                       const double &t,
                                       const double &min,
                                       const double &v);

Rcpp::NumericVector eadelta_intervals(const int &k,
                                      const double &x,
                                      const double &y,
                                      const double &s,
                                      const double &t,
                                      const double &min,
                                      const double &v);

#endif