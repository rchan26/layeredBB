#include "../inc/bessel_layer_sim.hpp"
#include "../inc/cauchy_sums.hpp"
#include "../inc/coins_cauchy.hpp"

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

int bessel_layer_simulation(const double &x, const double &y,
                            const double &s, const double &t,
                            Rcpp::NumericVector &a)
{
    // initialise Bessel layer
    int l = 0;

    while (true) {
        // flip gamma coin to determine if BB stays within (min(x,y)-a[l], max(x,y)+a[l])
        if (!(x==y && a.at(l)==0)) {
            if (gamma_coin(x,y,s,t,std::min(x,y)-a.at(l),std::max(x,y)+a.at(l),0)) {
                // if true, then return current Bessel layer (l)
                return l;
            }
        }
        
        // if false, l = l+1
        l += 1;
        // need to check that we are not at the end of the sequence a
        // in this case, we extend the sequence of a, and try carry on to find a layer
        if (l >= a.size()) {
            double last = a.back();
            auto size = a.size();
            for (int i = 0; i < size; ++i) {
                a.push_back(a[i] + last);
            }
        }
    }
}
                                
