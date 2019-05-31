#include "../inc/coins_cauchy.hpp"
#include "../inc/cauchy_sums.hpp"
#include <random>
#include <Rcpp.h>
#include <iomanip>

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
double product_vector_elements(const Rcpp::NumericVector &vect) {
  double prod = 1;
  for (const auto &element: vect) {
    prod *= element;
  }
  return prod;
}

// [[Rcpp::export]]
bool gamma_coin(const double &x, const double &y, 
                const double &s, const double &t,
                const double &l, const double &v,
                int k)
{
  // simulating from a uniform distribution
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  
  // calculate the current interval (S_{2k+1}^{gamma}, S_{2k}^{gamma})
  Rcpp::NumericVector current = calc_SgammaK_intervals(k,x,y,s,t,l,v);
  double left = current.at(0);
  double right = current.at(1);
  
  // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
  while (left < u && u < right) {
    k = k+1;
    right = right - (sigma(k,x,y,s,t,l,v) - phi(k,x,y,s,t,l,v));
    left = right - sigma(k+1,x,y,s,t,l,v);
  }
  
  if (u <= left) {
    return true;
  } else {
    return false;
  }
}

// [[Rcpp::export]]
bool gamma_coin_intervals(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, const Rcpp::NumericVector &t,
                          const double &l, const double &v, int k)
{
  // check if vector lengths are all the same
  if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
    stop("error in gamma_coin_intervals: vector lengths are not equal");
  }
  
  // simulating from a uniform distribution
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  
  Rcpp::NumericVector left, right;
  for (int i=0; i < x.size(); ++i) {
    Rcpp::NumericVector current = calc_SgammaK_intervals(k,x.at(i),y.at(i),s.at(i),t.at(i),l,v);
    left.push_back(current.at(0));
    right.push_back(current.at(1));
  }
  
  double left_product = product_vector_elements(left);
  double right_product = product_vector_elements(right);
  while (left_product < u && u < right_product) {
    k = k+1;
    for (int i=0; i < x.size(); ++i) {
      right.at(i) = right.at(i) - (sigma(k,x.at(i),y.at(i),s.at(i),t.at(i),l,v) - 
        phi(k,x.at(i),y.at(i),s.at(i),t.at(i),l,v));
      left.at(i) = right.at(i) - sigma(k+1,x.at(i),y.at(i),s.at(i),t.at(i),l,v);
      // recalculate products
      left_product = product_vector_elements(left);
      right_product = product_vector_elements(right);
    }
  }
  
  if (u <= left_product) {
    return true;
  } else {
    return false;
  }
}

// [[Rcpp::export]]
bool delta_coin(const double &x, const double &y, 
                const double &s, const double &t,
                const double &min, const double &v,
                int k)
{
  // simulating from a uniform distribution
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  
  // calculate the current interval (S_{2k+1}^{delta}, S_{2k}^{delta})
  Rcpp::NumericVector current = calc_SdeltaK_intervals(k,x,y,s,t,min,v);
  double left = current.at(0);
  double right = current.at(1);
  
  // initialise denom variable incase we need to carry on calculating Cauchy sums
  double denom = 1 - exp(-2.0*(x-min)*(y-min)/(t-s));
  
  // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
  while (left < u && u < right) {
    k = k+1;
    if (std::min(x,y) > min) {
      right = right - (sigma(k,x,y,s,t,min,v) - phi(k,x,y,s,t,min,v))/denom;
      left = right - sigma(k+1,x,y,s,t,min,v)/denom;
    } else if (std::min(x,y) == min) {
      right = right - (psi(k,x,y,s,t,min,v) - chi(k,x,y,s,t,min,v))/fabs(x-y);
      left = right - psi(k+1,x,y,s,t,min,v)/fabs(x-y);
    }   
  }
  
  if (u <= left) {
    return true;
  } else {
    return false;
  }
}

// [[Rcpp::export]]
bool delta_coin_intervals(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, const Rcpp::NumericVector &t,
                          const double &min, const double &v, int k)
{
  // check if vector lengths are all the same
  if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
    stop("error in gamma_coin_intervals: vector lengths are not equal");
  }
  
  // simulating from a uniform distribution
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  
  Rcpp::NumericVector left, right;
  for (int i=0; i < x.size(); ++i) {
    Rcpp::NumericVector current = calc_SdeltaK_intervals(k,x.at(i),y.at(i),s.at(i),t.at(i),min,v);
    left.push_back(current.at(0));
    right.push_back(current.at(1));
  }
  
  // initialise denom variable incase we need to carry on calculating Cauchy sums
  Rcpp::NumericVector denom;
  for (int i=0; i < x.size(); ++i) {
    if (std::min(x.at(i),y.at(i)) > min) {
      denom.push_back(1 - exp(-2.0*(x.at(i)-min)*(y.at(i)-min)/(t.at(i)-s.at(i))));
    } else {
      denom.push_back(1);
    }
  }
  
  // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
  double left_product = product_vector_elements(left);
  double right_product = product_vector_elements(right);
  while (left_product < u && u < right_product) {
    k = k+1;
    for (int i=0; i < x.size(); ++i) {
      if (std::min(x.at(i), y.at(i)) > min) {
        right.at(i) = right.at(i) - (sigma(k,x.at(i),y.at(i),s.at(i),t.at(i),min,v) -
          phi(k,x.at(i),y.at(i),s.at(i),t.at(i),min,v))/denom.at(i);
        left.at(i) = right.at(i) - sigma(k+1,x.at(i),y.at(i),s.at(i),t.at(i),min,v)/denom.at(i);
      } else if (std::min(x.at(i), y.at(i)) == min) {
        right.at(i) = right.at(i) - (psi(k,x.at(i),y.at(i),s.at(i),t.at(i),min,v) - 
          chi(k,x.at(i),y.at(i),s.at(i),t.at(i),min,v))/fabs(x.at(i)-y.at(i));
        left.at(i) = right.at(i) - psi(k+1,x.at(i),y.at(i),s.at(i),t.at(i),min,v)/fabs(x.at(i)-y.at(i));
      }
      // recalculate products
      left_product = product_vector_elements(left);
      right_product = product_vector_elements(right);
    }
  }
  
  if (u <= left_product) {
    return true;
  } else {
    return false;
  }
}
