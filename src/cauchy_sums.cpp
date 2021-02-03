#include "../inc/cauchy_sums.hpp"

using namespace Rcpp;

//' Sigma_Bar (Equation 141 in ST329)
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
//' @return real value: easigma_bar evaluated at point j
//'
//' @examples
//' easigma_bar(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double easigma_bar(const double &j,
                   const double &x,
                   const double &y,
                   const double &s,
                   const double &t,
                   const double &l,
                   const double &v)
{
  double P = -2.0/(t-s);
  double D = std::abs(v-l);
  double landv = std::min(l,v);
  return exp(P * (D*j + landv - x) * (D*j + landv - y));
}

//' Sigma (Equation 140 in ST329)
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
//' easigma(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double easigma(const double &j,
               const double &x,
               const double &y,
               const double &s,
               const double &t,
               const double &l,
               const double &v)
{
  return easigma_bar(j,x,y,s,t,l,v) + easigma_bar(j,-x,-y,s,t,-l,-v);
}

//' Phi_Bar (Equation 142 in ST329)
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
//' eaphi_bar(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double eaphi_bar(const double &j,
                 const double &x,
                 const double &y,
                 const double &s,
                 const double &t,
                 const double &l,
                 const double &v)
{
  double P = -2.0*j/(t-s);
  double D = std::abs(v-l);
  return exp(P * (D*D*j + D*(x-y)));
}

//' Phi (Equation 140 in ST329)
//'
//' This function evaluates the eaphi function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: eaphi evaluated at point j
//'
//' @examples
//' eaphi(j = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//' 
//' @export
// [[Rcpp::export]] 
double eaphi(const double &j,
             const double &x,
             const double &y,
             const double &s,
             const double &t,
             const double &l,
             const double &v)
{
  return eaphi_bar(j,x,y,s,t,l,v) + eaphi_bar(j,-x,-y,s,t,-l,-v);
}

//' Psi (Equation 155 and 162 in ST329)
//'
//' This function evaluates the eapsi function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param xoy maximum of x and y where x and y are the start and end values of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: eapsi evaluated at point j
//'
//' @examples
//' eapsi(j = 1, xoy = 0, s = 0, t = 1, min = -2, v = 1)
//' 
//' @export
// [[Rcpp::export]] 
double eapsi(const double &j,
             const double &xoy,
             const double &s,
             const double &t,
             const double &min,
             const double &v)
{
  double P = 2 * std::abs(v-min) * j;
  return (P - (xoy - min)) * exp(-(P/(t-s)) * (std::abs(v-min)*j - (xoy-min)));
}

//' Chi (Equation 156 and 163 in ST329)
//'
//' This function evaluates the chi function used to simulate Bessel Layers in the infinite sums
//'
//' @param j real value
//' @param xoy maximum of x and y where x and y are the start and end values of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: chi evaluated at point j
//'
//' @examples
//' eachi(j = 1, xoy = 0, s = 0, t = 1, min = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
double eachi(const double &j,
             const double &xoy,
             const double &s,
             const double &t,
             const double &min,
             const double &v)
{
  double P = 2 * std::abs(v-min) * j;
  return (P + (xoy - min)) * exp(-(P/(t-s)) * (std::abs(v-min)*j + (xoy-min)));
}

// ---------- Calculate S_{n}^{gamma} and S_{n}^{delta}

//' Gamma (Corollary 2 and Algorithm 26 in ST329) 
//'
//' This function evaluates the gamma function, S_{n}^{gamma},
//' used to simulate Bessel Layers in the infinite sums
//'
//' @param n integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: S_{n}^{gamma} evaluated at n
//'
//' @examples
//' eagamma(n = 1, x = 0, y = 0, s = 0, t = 1, l = -0.4, v = 0.8)
//'
//' @export
// [[Rcpp::export]]
double eagamma(const int &n,
               const double &x,
               const double &y,
               const double &s,
               const double &t,
               const double &l,
               const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < l) {
    return 0;
  } else if (std::max(x,y) > v) {
    return 0;
  }
  if (n % 2 == 0) {
    int k = n / 2;
    double zeta = 0;
    for (int j=1; j <= k; ++j) {
      zeta += (easigma(j,x,y,s,t,l,v) - eaphi(j,x,y,s,t,l,v));
    }
    return 1 - zeta;
  } else {
    if (n > 1) {
      int k = (n-1)/2;
      double zeta = 0;
      for (int j=1; j <= k; ++j) {
        zeta += (easigma(j,x,y,s,t,l,v) - eaphi(j,x,y,s,t,l,v));
      }
      return 1 - zeta - easigma(k+1,x,y,s,t,l,v);
    } else {
      return 1 - easigma(1,x,y,s,t,l,v);
    }
  }
}

//' Delta_1 (Corollary 3 in ST329)
//'
//' This function evaluates the delta_1 function, S_{n}^{delta,1}
//' used to simulate Bessel Layers in the infinite sums
//'
//' @param n integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: S_{n}^{delta,1} evaluated at n
//'
//' @examples
//' eadelta1(n = 1, x = 0, y = 0.8, s = 0, t = 1, min = -2, v = 2)
//'
//' @export
// [[Rcpp::export]]
double eadelta1(const int &n,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &min,
                const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < min) {
    return 0;
  } else if (std::max(x,y) > v) {
    return 0;
  }
  // check if min(x,y) == min. If so, should be using eadelta2, so return NaN
  if (std::min(x,y) == min) {
    return R_NaN;
  }
  return eagamma(n,x,y,s,t,min,v) / (1 - exp(-2.0*(x-min)*(y-min)/(t-s)));
}

//' Delta_2 (Corollary 4 in ST329)
//'
//' This function evaluates the delta_2function, S_{n}^{delta,2}
//' used to simulate Bessel Layers in the infinite sums
//'
//' @param n integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: S_{n}^{delta,2} evaluated at n
//'
//' @examples
//' eadelta2(n = 1, x = -2, y = 0.8, s = 0, t = 1, min = -2, v = 2)
//'
//' @export
// [[Rcpp::export]]
double eadelta2(const int &n,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &min,
                const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < min) {
    return 0;
  } else if (std::max(x,y) > v) {
    return 0;
  }
  // check if min(x,y) > min. If so, should be using eadelta1, so return NaN
  if (std::min(x,y) > min) {
    return R_NaN;
  }
  double xoy = std::max(x,y);
  double denom = xoy - min;
  if (n % 2 == 0) {
    int k = n / 2;
    double sum = 0;
    for (int j=1; j <= k; ++j) {
      sum += (eapsi(j,xoy,s,t,min,v) - eachi(j,xoy,s,t,min,v)) / denom;
    }
    return 1 - sum;
  } else {
    if (n > 1) {
      int k = (n-1)/2;
      double sum = 0;
      for (int j=1; j <= k; ++j) {
        sum += (eapsi(j,xoy,s,t,min,v) - eachi(j,xoy,s,t,min,v)) / denom;
      }
      return 1 - sum - (eapsi(k+1,xoy,s,t,min,v) / denom);
    } else {
      return 1 - (eapsi(1,xoy,s,t,min,v) / denom);
    }
  } 
}

//' Delta (Corollary 3 and 4 in ST329)
//'
//' This function evaluates the delta function, S_{n}^{delta}
//' used to simulate Bessel Layers in the infinite sums
//'
//' @param n integer value
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' 
//' @return real value: S_{n}^{delta} evaluated at n
//'
//' @examples
//' # example where min(x,y) > min (delta_1 case)
//' eadelta(n = 1, x = 0,  y = 0.8, s = 0, t = 1, min = -2, v = 2)
//' eadelta1(n = 1, x = 0, y = 0.8, s = 0, t = 1, min = -2, v = 2)
//' # example where min(x,y) == min
//' eadelta(n = 1, x = -2, y = 0.8, s = 0, t = 1, min = -2, v = 2)
//' eadelta2(n = 1, x = -2, y = 0.8, s = 0, t = 1, min = -2, v = 2)
//'
//' @export
// [[Rcpp::export]]
double eadelta(const int &n,
               const double &x,
               const double &y,
               const double &s,
               const double &t,
               const double &min,
               const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < min) {
    return 0;
  } else if (std::max(x,y) > v) {
    return 0;
  }
  // decide whether to use delta_1 (if min(x,y) > min) or delta_2 (if min(x,y) == min)
  if (std::min(x,y) > min) {
    return eadelta1(n,x,y,s,t,min,v);
  } else {
    return eadelta2(n,x,y,s,t,min,v);
  }
}

// ---------- Calculate intervals needed for EA. i.e. calculating [S_{2k+1}, S_{2k}]

//' Calculate interval: [S^{gamma}_{2k+1}, S^{gamma}_{2k}] (Corollary 2 and Algorithm 26 in ST329)
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
//' eagamma_intervals(k = 1, x = 0, y = 0, s = 0, t = 1, l = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eagamma_intervals(const int &k,
                                      const double &x,
                                      const double &y,
                                      const double &s,
                                      const double &t,
                                      const double &l,
                                      const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < l) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  } else if (std::max(x,y) > v) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  }
  if (k == 0) {
    return Rcpp::NumericVector::create(1-easigma(1.0,x,y,s,t,l,v), 1);
  } else {
    double zeta = 0;
    for (int j=1; j <= k; ++j) {
      zeta += (easigma(j,x,y,s,t,l,v) - eaphi(j,x,y,s,t,l,v));
    }
    double S_2k = 1 - zeta;
    double S_2k_plus_1 = S_2k - easigma(k+1,x,y,s,t,l,v);
    return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
  }
}

//' Calculate interval: [S^{delta,1}_{2k+1}, S^{delta,1}_{2k}] (Corollary 3 in ST329)
//'
//' This function calculates the interval [S^{delta,1}_{2k+1}, S^{delta, 1}_{2k}] 
//' (case where min(x,y) > min)
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
//' eadelta1_intervals(k = 1, x = 0, y = 0, s = 0, t = 1, min = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eadelta1_intervals(const int &k,
                                       const double &x,
                                       const double &y,
                                       const double &s,
                                       const double &t,
                                       const double &min,
                                       const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < min) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  } else if (std::max(x,y) > v) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  }
  // check if min(x,y) == min. If so, should be using eadelta2_intervals, so return c(NaN, NaN)
  if (std::min(x,y) == min) {
    return Rcpp::NumericVector::create(R_NaN, R_NaN);
  }
  double denom = 1 - exp(-2.0*(x-min)*(y-min)/(t-s));
  Rcpp::NumericVector eagamma = eagamma_intervals(k,x,y,s,t,min,v);
  return Rcpp::NumericVector::create(eagamma.at(0)/denom, eagamma.at(1)/denom);
}

//' Calculate interval: [S^{delta,2}_{2k+1}, S^{delta,2}_{2k}] (Corollary 4 in ST329)
//'
//' This function calculates the interval [S^{delta,2}_{2k+1}, S^{delta, 2}_{2k}] 
//' (case where min(x,y) == min)
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
//' eadelta2_intervals(k = K, x = -2, y = 0, s = 0, t = 0, min = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eadelta2_intervals(const int &k,
                                       const double &x,
                                       const double &y,
                                       const double &s,
                                       const double &t,
                                       const double &min,
                                       const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < min) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  } else if (std::max(x,y) > v) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  }
  // check if min(x,y) > min. If so, should be using eadelta1_intervals, so return c(NaN, NaN)
  if (std::min(x,y) > min) {
    return Rcpp::NumericVector::create(R_NaN, R_NaN);
  }
  // check that given k is large enough
  double D = std::abs(v-min);
  double K = sqrt(t-s + D*D) / (2*D);
  if (k < K) {
    stop("layeredBB::eadelta2_intervals: given k is too small");
  }
  double sum = 0;
  double xoy = std::max(x,y);
  double denom = xoy - min;
  for (int j=1; j <= k; ++j) {
    sum += (eapsi(j,xoy,s,t,min,v) - eachi(j,xoy,s,t,min,v)) / denom;
  }
  double S_2k = 1- sum;
  double S_2k_plus_1 = S_2k - (eapsi(k+1,xoy,s,t,min,v) / denom);
  return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
}

//' Calculate interval: [S^{delta}_{2k+1}, S^{delta}_{2k}] (Algorithm 28 in ST329)
//'
//' This function calculates the interval [S^{delta}_{2k+1}, S^{delta}_{2k}] 
//' (case where min(x,y) > min or where min(x,y) == min)
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
//' eadelta1_intervals(k = 1, x = 0, y = 0, s = 0, t = 1, min = -2, v = 1)
//'
//' # case where min(x,y) == min
//' K = ceiling(sqrt((1)+(abs(1-(-2))*abs(1-(-2))))/(2*abs(1-(-2))))
//' eadelta_intervals(k = K, x = -2, y = 0, s = 0, t = 0, min = -2, v = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eadelta_intervals(const int &k,
                                      const double &x,
                                      const double &y,
                                      const double &s,
                                      const double &t,
                                      const double &min,
                                      const double &v)
{
  // check if x or y already outside bounds
  if (std::min(x,y) < min) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  } else if (std::max(x,y) > v) {
    return Rcpp::NumericVector::create(0.0, 0.0);
  }
  // decide whether to use delta_1 (if min(x,y) > min) or delta_2 (if min(x,y) == min)
  if (std::min(x,y) > min) {
    return eadelta1_intervals(k,x,y,s,t,min,v);
  } else {
    return eadelta2_intervals(k,x,y,s,t,min,v);
  }
}
