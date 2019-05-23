#include <Rcpp.h>
#include "../inc/cauchy_sums.hpp"

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
double sigma_bar(const double &j, const double &x, const double &y,
                 const double &s, const double &t,
                 const double &l, const double &v)
{
    return exp((-2.0/(t-s))*((fabs(v-l)*j)+std::min(l,v)-x)*((fabs(v-l)*j)+std::min(l,v)-y));
}

// [[Rcpp::export]]
double sigma(const double &j, const double &x, const double &y, 
                         const double &s, const double &t,
                         const double &l, const double &v)
{
    return sigma_bar(j,x,y,s,t,l,v) + sigma_bar(j,-x,-y,s,t,-l,-v);
}

// [[Rcpp::export]]
double phi_bar(const double &j, const double &x, const double &y,
                           const double &s, const double &t,
                           const double &l, const double &v)
{
    return exp(-(2.0*j/(t-s)) * (fabs(v-l)*fabs(v-l)*j + fabs(v-l)*(x-y))); 
}

// [[Rcpp::export]]
double phi(const double &j, const double &x, const double &y, 
                       const double &s, const double &t,
                       const double &l, const double &v)
{
    return phi_bar(j,x,y,s,t,l,v) + phi_bar(j,-x,-y,s,t,-l,-v);
}

// [[Rcpp::export]]
double psi(const double &j, const double &x, const double &y,
                       const double &s, const double &t,
                       const double &min, const double &v)
{
    return ((2*fabs(v-min)*j)-std::max(x,y)-min)*exp(-(2.0*fabs(v-min)*j/(t-s))*((fabs(v-min)*j)-(std::max(x,y)-min)));
}

// [[Rcpp::export]]
double chi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v)
{
    return ((2*fabs(v-min)*j)+std::max(x,y)-min)*exp(-(2.0*fabs(v-min)*j/(t-s))*((fabs(v-min)*j)+(std::max(x,y)-min)));
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SgammaK_intervals(const int &k, const double &x, const double &y,
                                           const double &s, const double &t,
                                           const double &l, const double &v)
{
    // function calculates (S_{2k+1}^{gamma}, S_{2k}^{gamma}) for a given k

    // if k==0, then we need to calculate (S_{1}^{gamma}, S_{0}^{gamma}) = (1-sigma(1), 1)
    if (k == 0) {
        return Rcpp::NumericVector::create(1-sigma(1.0,x,y,s,t,l,v), 1);
    } 

    // calculating S_{2k}^{gamma} = 1 - sum_{j=1}^{k} sigma(j)-phi(j)
    double S_2k = 1;
    for (int j=1; j <= k; ++j) {
        S_2k -= (sigma(j,x,y,s,t,l,v) - phi(j,x,y,s,t,l,v));
    }
    // calculating S_{2k+1}^{gamma} = S_{2k} - sigma(k+1)
    double S_2k_plus_1 = S_2k - sigma(k+1,x,y,s,t,l,v);

    return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_1_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v)
{
    // function calculates (S_{2k+1}^{delta,1}, S_{2k}^{delta,1}) for a given k
    // S_{k}^{delta,1} = S_{k}^{gamma}/denom
    double denom = 1 - exp(-2.0*(x-min)*(y-min)/(t-s));
    // use calc_SgammaK_intervals and divide the result by denom
    Rcpp::NumericVector SgammaK = calc_SgammaK_intervals(k,x,y,s,t,min,v);

    return Rcpp::NumericVector::create(SgammaK.at(0)/denom, SgammaK.at(1)/denom);
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_2_intervals(const int &k, const double &x, const double &y,
                                               const double &s, const double &t,
                                               const double &min, const double &v)
{
    // function calculates (S_{2k+1}^{delta,2}, S_{2k}^{delta,2}) for a given k
    // checking k is large enough for this to be valid
    double K = sqrt((t-s)+(fabs(v-min)*fabs(v-min)))/(2*fabs(v-min));
    if (k < K) {
        Rcout << "given k is too small. \n";
    }

    // calculating S_{2k+1}^{delta,2} = 1 - (psi(j)-chi(j)) / abs(x-y)
    double S_2k = 1;
    for (int j=1; j <= k; ++j) {
        S_2k -= (psi(j,x,y,s,t,min,v) - chi(j,x,y,s,t,min,v));
    }
    double S_2k_plus_1 = S_2k - psi(k+1,x,y,s,t,min,v)/fabs(x-y);

    return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v)
{
    if (std::min(x,y) > min) {
        return calc_SdeltaK_1_intervals(k,x,y,s,t,min,v);
    } else if (std::min(x,y) == min) {
        return calc_SdeltaK_2_intervals(k,x,y,s,t,min,v);
    } else {
        Rcout << "min(x,y) < min: given minimum point is not the minimum of the BB.";
    }
}

