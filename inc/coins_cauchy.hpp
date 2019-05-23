#ifndef COINS_CAUCHY
#define COINS_CAUCHY

#include <random>

using namespace Rcpp;

// forward declaration for taking products of vector elements
double product_vector_elements(const Rcpp::NumericVector &vect);

// forward declaration for flipping gamma coin
bool gamma_coin(const double &x, const double &y, 
                const double &s, const double &t,
                const double &l, const double &v,
                int k);

// forward declaration for flipping gamma coin between each interval
bool gamma_coin_intervals(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, const Rcpp::NumericVector &t,
                          const double &l, const double &v, int k);

// forward declaration for flipping delta coin
bool delta_coin(const double &x, const double &y, 
                const double &s, const double &t,
                const double &min, const double &v,
                int k);

// forward decalartion for flipping delta coin between each interval
bool delta_coin_intervals(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, const Rcpp::NumericVector &t,
                          const double &l, const double &v, int k);

#endif