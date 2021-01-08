#ifndef COINS_CAUCHY
#define COINS_CAUCHY

#include <Rcpp.h>
#include <random>
#include <iomanip>

using namespace Rcpp;

double product_vector_elements(const Rcpp::NumericVector &vect);

bool gamma_coin(int k,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &l,
                const double &v);

bool gamma_coin_intervals(int k,
                          const Rcpp::NumericVector &X,
                          const Rcpp::NumericVector &times,
                          const double &l,
                          const double &v);

bool delta_coin(int k,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &min,
                const double &v);

bool delta_coin_intervals(int k,
                          const Rcpp::NumericVector &X,
                          const Rcpp::NumericVector &times,
                          const double &min,
                          const double &v);

#endif