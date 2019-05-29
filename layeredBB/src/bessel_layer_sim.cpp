#include "../inc/bessel_layer_sim.hpp"
#include "../inc/cauchy_sums.hpp"
#include "../inc/coins_cauchy.hpp"
#include <vector>
#include <algorithm>
#include <random>

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
List bessel_layer_simulation(const double &x, const double &y,
                             const double &s, const double &t,
                             Rcpp::NumericVector &a)
{
  // initialise Bessel layer
  int l = 0;
  
  while (true) {
    // flip gamma coin to determine if BB stays within (min(x,y)-a[l], max(x,y)+a[l])
    if (!(x==y && a.at(l)==0)) {
      if (gamma_coin(x,y,s,t,std::min(x,y)-a.at(l),std::max(x,y)+a.at(l),0)) {
        // if true, then return current Bessel layer (l) and current vector (a) as a list
        // we return (l+1), since indicies start from 0 in C++ but start from 1 in R
        return List::create(_["a"] = a, _["l"] = l+1);
      }
    }
    
    // if false, l = l+1
    l += 1;
    // need to check that we are not at the end of the sequence (a)
    // in this case, we extend the sequence of (a), and try carry on to find a layer
    if (l >= a.size()) {
      // find the size of the vector (a) and set it to (size)
      auto size = a.size();
      // find the last element of the vector (a) and set it to (last)
      // since indices start from 0 in C++, this is the (size-1)-th element
      double last = a[(size-1)];
      for (int i = 0; i < size; ++i) {
        a.push_back(a[i] + last);
      }
    }
  }
}
