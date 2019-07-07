#include "../inc/cauchy_sums.hpp"
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//' Sigma_Bar
//'
//' This function evaluates the sigma_bar function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: sigma_bar evaluated at point j
//'
//' @examples
//' sigma_bar(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double sigma_bar(const double &j, const double &x, const double &y,
                 const double &s, const double &t,
                 const double &l, const double &v)
{
  return exp((-2.0/(t-s))*((fabs(v-l)*j)+std::min(l,v)-x)*((fabs(v-l)*j)+std::min(l,v)-y));
}

//' Sigma
//'
//' This function evaluates the sigma function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: sigma evaluated at point j
//'
//' @examples
//' sigma(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double sigma(const double &j, const double &x, const double &y, 
             const double &s, const double &t,
             const double &l, const double &v)
{
  return sigma_bar(j,x,y,s,t,l,v) + sigma_bar(j,-x,-y,s,t,-l,-v);
}

//' Phi_Bar
//'
//' This function evaluates the phi_bar function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: phi_bar evaluated at point j
//'
//' @examples
//' phi_bar(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double phi_bar(const double &j, const double &x, const double &y,
               const double &s, const double &t,
               const double &l, const double &v)
{
  return exp(-(2.0*j/(t-s)) * (fabs(v-l)*fabs(v-l)*j + fabs(v-l)*(x-y))); 
}

//' Phi
//'
//' This function evaluates the phi function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: phi evaluated at point j
//'
//' @examples
//' phi(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//' 
//' @export
// [[Rcpp::export]] 
double phi(const double &j, const double &x, const double &y, 
           const double &s, const double &t,
           const double &l, const double &v)
{
  return phi_bar(j,x,y,s,t,l,v) + phi_bar(j,-x,-y,s,t,-l,-v);
}

//' Psi
//'
//' This function evaluates the psi function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: psi evaluated at point j
//'
//' @examples
//' psi(j = 1, x = 0, y = 0, s = 0, t = 1, min = -2, v = 1)
//' 
//' @export
// [[Rcpp::export]] 
double psi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v)
{
  return (((2*fabs(v-min)*j)-(std::max(x,y)-min))*(exp(-(2.0*fabs(v-min)*j/(t-s))*((fabs(v-min)*j)-(std::max(x,y)-min)))));
}

//' Chi
//'
//' This function evaluates the psi function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: chi evaluated at point j
//'
//' @examples
//' chi(j = 1, x = 0, y = 0, s = 0, t = 1, min = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double chi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v)
{
  return (((2*fabs(v-min)*j)+(std::max(x,y)-min))*(exp(-(2.0*fabs(v-min)*j/(t-s))*((fabs(v-min)*j)+(std::max(x,y)-min)))));
}

//' Calculate interval: [S^{gamma}_{2k+1}, S^{gamma}_{2k}]
//'
//' This function calculates the interval [S^{gamma}_{2k+1}, S^{gamma}_{2k}] for given k
//'
//' @param k integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @return vector of two values, S^{gamma}_{2k+1} and S^{gamma}_{2k}
//'
//' @examples
//' calc_SgammaK_intervals(k = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
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

//' Calculate interval: [S^{delta,1}_{2k+1}, S^{delta,1}_{2k}]
//'
//' This function calculates the interval [S^{delta,1}_{2k+1}, S^{delta, 1}_{2k}] (case where min(x,y) > min)
//'
//' @param k integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @return vector of two values, S^{delta,1}_{2k+1} and S^{delta,1}_{2k}
//'
//' @examples
//' calc_SdeltaK_1_intervals(k = 1, x = 0, y = 0, s = 0, t = 1, min = -2, v = 1)
//'
//' @export
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

//' Calculate interval: [S^{delta,2}_{2k+1}, S^{delta,2}_{2k}]
//'
//' This function calculates the interval [S^{delta,2}_{2k+1}, S^{delta, 2}_{2k}] (case where min(x,y) == min)
//'
//' @param k integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @return vector of two values, S^{delta,2}_{2k+1} and S^{delta,2}_{2k}
//'
//' @examples
//' K = ceiling(sqrt((1)+(abs(1-(-2))*abs(1-(-2))))/(2*abs(1-(-2))))
//' calc_SdeltaK_2_intervals(k = K, x = -2, y = 0, s = 0, t = 0, min = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_2_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v)
{
  // function calculates (S_{2k+1}^{delta,2}, S_{2k}^{delta,2}) for a given k
  // checking k is large enough for this to be valid
  double K = sqrt((t-s)+(fabs(v-min)*fabs(v-min)))/(2*fabs(v-min));
  if (k < K) {
    stop("error in calc_SdeltaK_2_intervals: given k is too small");
  }
  
  // calculating S_{2k+1}^{delta,2} = 1 - (psi(j)-chi(j)) / abs(x-y)
  double S_2k = 1;
  for (int j=1; j <= k; ++j) {
    S_2k -= ((psi(j,x,y,s,t,min,v)-chi(j,x,y,s,t,min,v))/fabs(x-y));
  }
  double S_2k_plus_1 = S_2k - (psi(k+1,x,y,s,t,min,v)/fabs(x-y));
   
  return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
}

//' Calculate interval: [S^{delta}_{2k+1}, S^{delta}_{2k}]
//'
//' This function calculates the interval [S^{delta}_{2k+1}, S^{delta}_{2k}] (case where min(x,y) > min or where min(x,y) == min)
//'
//' @param k integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @return vector of two values, S^{delta}_{2k+1} and S^{delta}_{2k}
//'
//' @examples
//' # case where min(x,y ) > min
//' calc_SdeltaK_1_intervals(k = 1, x = 0, y = 0, s = 0, t = 1, min = -2, v = 1)
//'
//' # case where min(x,y) == min
//' K = ceiling(sqrt((1)+(abs(1-(-2))*abs(1-(-2))))/(2*abs(1-(-2))))
//' calc_SdeltaK_intervals(k = K, x = -2, y = 0, s = 0, t = 0, min = -2, v = 1)
//'
//' @export
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
    stop("error in calc_SdeltaK_intervals: min(x,y) < min - given minimum point is not the minimum of the BB");
  }
}

