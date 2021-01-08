#include "../inc/coins_cauchy.hpp"
#include "../inc/cauchy_sums.hpp"

using namespace Rcpp;

//' Find product of a vector
//'
//' This function product of the elements of a numerical vector
//'
//' @param vect numerical vector
//'
//' @return product of element in vector given
//'
//' @examples
//' # returns 120
//' find_max(c(1,2,3,4,5)) 
//'
//' @export
// [[Rcpp::export]]
double product_vector_elements(const Rcpp::NumericVector &vect) {
  double prod = 1;
  for (const auto &element: vect) {
    prod *= element;
  }
  return prod;
}

//' Gamma coin flipper (Algorithm 26 in ST329)
//'
//' Flips 'Gamma coin'; uses the Cauchy sequence S^{gamma}_{k} to 
//' determine whether or not the Brownian bridge starting at x, ending at y, between [s,t]
//' remains in interval [l,v]
//'
//' @param k integer value starting index for calculating the intervals
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @examples
//' gamma_coin(k = 0, x = 0, y = 0, s = 0, t = 1, l = -0.5, v = 0.5)
//'
//' @return boolean value: if T, accept probability that Brownian bridge remains 
//'         in [l,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool gamma_coin(int k,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &l,
                const double &v)
{
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  // calculate the current interval (S_{2k+1}^{gamma}, S_{2k}^{gamma}) at k
  Rcpp::NumericVector current = eagamma_intervals(k,x,y,s,t,l,v);
  // set left = S_{2k+1}^{gamma} and right = S_{2k}^{gamma}
  double left = current.at(0);
  double right = current.at(1);
  // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
  // keep increasing k until u no longer falls in the interval between left and right
  // i.e. between (S_{2k+1}^{gamma}, S_{2k}^{gamma})
  while (left < u && u < right) {
    k = k+1;
    right = right - (easigma(k,x,y,s,t,l,v) - eaphi(k,x,y,s,t,l,v));
    left = right - easigma(k+1,x,y,s,t,l,v);
  }
  if (u <= left) {
    return true;
  } else {
    return false;
  }
}

// //' Gamma coin flipper for intervals
// //'
// //' Flips 'Gamma coin' for intervals; takes the product of the Cauchy sequence S^{gamma}_{k} to
// //' determine whether or not the Brownian bridge remains in the interval [l,v]
// //' Vectors x, y, s, t should all be the same length, L, where for i = 1, ..., L, the Brownian Bridge skeleton
// //' we have is broken up so that x[i] goes to y[i] between s[i] and t[i] - see example
// //'
// //' @param k integer value starting index for calculating the intervals
// //' @param x vector of values
// //' @param y vector of values
// //' @param s vector of values
// //' @param t vector of values
// //' @param l lower bound of Brownian bridge
// //' @param v upper bound of Brownian bridge
// //'
// //' @examples
// //' # setting up vectors
// //' x_vect <- c(); y_vect <- c(); s_vect <- c(); t_vect <- c()
// //' brownian_bridge <- matrix(c(0, 0, -0.2, 0.4, 1, 1), ncol = 3, nrow = 2)
// //' for (i in 1:(ncol(brownian_bridge)-1)) {
// //'   x_vect[i] <- brownian_bridge[1,i]
// //'   y_vect[i] <- brownian_bridge[1,(i+1)]
// //'   s_vect[i] <- brownian_bridge[2,i]
// //'   t_vect[i] <- brownian_bridge[2,(i+1)]
// //' }
// //'
// //' # flip gamma coin whether or not Brownian bridge remains in [-0.5, 1.5]
// //' gamma_coin_intervals(x = x_vect, y = y_vect, s = s_vect, t = t_vect, l = -0.5, v = 1.5, k = 1)
// //'
// //' @return boolean value: if T, accept probability that Brownian bridge remains in [l,v], otherwise reject
// //'
// //' @export
// // [[Rcpp::export]]
// bool gamma_coin_intervals(int k,
//                           const Rcpp::NumericVector &x,
//                           const Rcpp::NumericVector &y,
//                           const Rcpp::NumericVector &s,
//                           const Rcpp::NumericVector &t,
//                           const double &l,
//                           const double &v)
// {
//   // check if vector lengths are all the same
//   if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
//     stop("layeredBB::gamma_coin_intervals: vector lengths are not equal");
//   }
//   double u = Rcpp::runif(1, 0.0, 1.0)[0];
//   Rcpp::NumericVector left(x.size()), right(x.size());
//   for (int i=0; i < x.size(); ++i) {
//     Rcpp::NumericVector current = eagamma_intervals(k,x.at(i),y.at(i),s.at(i),t.at(i),l,v);
//     left.at(i) = current.at(0);
//     right.at(i) = current.at(1);
//   }
//   double left_product = product_vector_elements(left);
//   double right_product = product_vector_elements(right);
//   while (left_product < u && u < right_product) {
//     k = k+1;
//     for (int i=0; i < x.size(); ++i) {
//       double sigma_at_k = easigma(k,x.at(i),y.at(i),s.at(i),t.at(i),l,v);
//       double phi_at_k = eaphi(k,x.at(i),y.at(i),s.at(i),t.at(i),l,v);
//       double sigma_at_k_plus_1 = easigma(k+1,x.at(i),y.at(i),s.at(i),t.at(i),l,v);
//       right.at(i) = right.at(i) - (sigma_at_k - phi_at_k);
//       left.at(i) = right.at(i) - sigma_at_k_plus_1;
//     }
//     // recalculate products
//     left_product = product_vector_elements(left);
//     right_product = product_vector_elements(right);
//   }
//   if (u <= left_product) {
//     return true;
//   } else {
//     return false;
//   }
// }

//' Delta coin flipper (Algorithm 28 in ST329)
//'
//' Flips 'Delta coin'; uses the Cauchy sequence S^{delta}_{k} to 
//' determine whether or not the Brownian bridge with minimum, min, 
//' starting at x, ending at y, between [s,t] remains in interval [l,v]
//'
//' @param k integer value starting index for calculating the intervals
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @examples 
//' delta_coin_intervals(k = 0,
//'                      x = 0.1,
//'                      y = 0.4,
//'                      s = 0,
//'                      t = 1,
//'                      min = -0.2,
//'                      v = 1.5)
//'                      
//' @return boolean value: if T, accept probability that Brownian bridge with 
//'         minimum, min, remains in [l,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool delta_coin(int k,
                const double &x,
                const double &y,
                const double &s,
                const double &t,
                const double &min,
                const double &v)
{
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  // calculate the current interval (S_{2k+1}^{delta}, S_{2k}^{delta})
  Rcpp::NumericVector current = eadelta_intervals(k,x,y,s,t,min,v);
  // set left = S_{2k+1}^{delta} and right = S_{2k}^{delta}
  double left = current.at(0);
  double right = current.at(1);
  // initialise variables incase we need to carry on calculating Cauchy sums
  double denom;
  double xoy;
  if (std::min(x,y) > min) {
    denom = 1 - exp(-2.0*(x-min)*(y-min)/(t-s));
  } else if (std::min(x,y) == min) {
    xoy = std::max(x,y);
    denom = std::abs(x-y);
  }
  // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
  // keep increasing k until u no longer falls in the interval between left and right
  // i.e. between (S_{2k+1}^{delta}, S_{2k}^{delta})
  while (left < u && u < right) {
    k = k+1;
    if (std::min(x,y) > min) {
      right = right - ((easigma(k,x,y,s,t,min,v) - eaphi(k,x,y,s,t,min,v)) / denom);
      left = right - (easigma(k+1,x,y,s,t,min,v) / denom);
    } else if (std::min(x,y) == min) {
      right = right - ((eapsi(k,xoy,s,t,min,v) - eachi(k,xoy,s,t,min,v)) / denom);
      left = right - (eapsi(k+1,xoy,s,t,min,v) / denom);
    }   
  }
  if (u <= left) {
    return true;
  } else {
    return false;
  }
}

//' Delta coin flipper for intervals (used for Algorithm 33 in ST329)
//'
//' Flips 'Delta coin' for intervals; takes the product of the Cauchy sequence S^{delta}_{k} to 
//' determine whether or not the Brownian bridge with minimum, min, remains in the interval [l,v]
//' Vectors x, y, s, t should all be the same length, L, where for i = 1, ..., L, the Brownian Bridge skeleton 
//' we have is broken up so that x[i] goes to y[i] between s[i] and t[i] - see example
//'
//' @param k integer value starting index for calculating the intervals
//' @param X vector of values of Brownian bridge
//' @param times vector of times
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//'
//' @examples
//' # setting up Brownian bridge variable
//' brownian_bridge <- matrix(c(0, 0, -0.2, 0.4, 0.3, 0.5, 1, 1),
//'                           ncol = 4, nrow = 2)
//' 
//' # flip delta coin whether or not Brownian bridge remains in [-0.2, 1.5]
//' d <- abs(1.5 - -0.2)
//' k <- ceiling(sqrt(1 + d^2)/(2*d))
//' delta_coin_intervals(k = k,
//'                      X = brownian_bridge[1,],
//'                      times = brownian_bridge[2,],
//'                      min = -0.2,
//'                      v = 1.5)
//'
//' @return boolean value: if T, accept probability that Brownian bridge with 
//'         minimum, min, remains in [min,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool delta_coin_intervals(int k,
                          const Rcpp::NumericVector &X,
                          const Rcpp::NumericVector &times,
                          const double &min,
                          const double &v)
{
  // check if vector lengths are all the same
  if (X.size() != times.size()) {
    stop("layeredBB::gamma_coin_intervals: vector lengths are not equal");
  }
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  int n = X.size()-1;
  Rcpp::NumericVector left(n), right(n);
  for (int i=0; i < n; ++i) {
    Rcpp::NumericVector current = eadelta_intervals(k,X.at(i),X.at(i+1),times.at(i),times.at(i+1),min,v);
    left.at(i) = current.at(0);
    right.at(i) = current.at(1);
  }
  // initialise variables incase we need to carry on calculating Cauchy sums
  Rcpp::NumericVector denom(n);
  Rcpp::NumericVector xoy(n);
  for (int i=0; i < n; ++i) {
    xoy.at(i) = std::max(X.at(i), X.at(i+1));
    if (std::min(X.at(i),X.at(i+1)) > min) {
      denom.at(i) = 1 - exp(-2.0*(X.at(i)-min)*(X.at(i+1)-min)/(times.at(i+1)-times.at(i)));
    } else if (std::min(X.at(i), X.at(i+1)) == min) {
      denom.at(i) = std::abs(X.at(i)-X.at(i+1));
    }
  }
  double left_product = product_vector_elements(left);
  double right_product = product_vector_elements(right);
  // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
  while (left_product < u && u < right_product) {
    k += 1;
    for (int i=0; i < n; ++i) {
      if (std::min(X.at(i), X.at(i+1)) > min) {
        double sigma_at_k = easigma(k, X.at(i), X.at(i+1), times.at(i), times.at(i+1), min, v);
        double phi_at_k = eaphi(k, X.at(i), X.at(i+1), times.at(i), times.at(i+1), min, v);
        double sigma_at_k_plus_1 = easigma(k+1, X.at(i), X.at(i+1), times.at(i), times.at(i+1), min, v);
        right.at(i) = right.at(i) - ((sigma_at_k - phi_at_k) / denom.at(i));
        left.at(i) = right.at(i) - (sigma_at_k_plus_1 / denom.at(i));
      } else if (std::min(X.at(i), X.at(i+1)) == min) {
        double eapsi_at_k = eapsi(k, xoy.at(i), times.at(i), times.at(i+1), min, v);
        double chi_at_k = eachi(k, xoy.at(i), times.at(i), times.at(i+1), min, v);
        double eapsi_at_k_plus_1 = eapsi(k+1, xoy.at(i), times.at(i), times.at(i+1), min, v);
        right.at(i) = right.at(i) - ((eapsi_at_k - chi_at_k) / denom.at(i));
        left.at(i) = right.at(i) - (eapsi_at_k_plus_1 / denom.at(i));
      }
    }
    // recalculate products
    left_product = product_vector_elements(left);
    right_product = product_vector_elements(right);
  }
  if (u <= left_product) {
    return true;
  } else {
    return false;
  }
}

// // [[Rcpp::export]]
// bool delta_coin_intervals_basic(int k,
//                                 const Rcpp::NumericVector &X,
//                                 const Rcpp::NumericVector &times,
//                                 const double &min,
//                                 const double &v)
// {
//   // check if vector lengths are all the same
//   if (X.size() != times.size()) {
//     stop("layeredBB::gamma_coin_intervals: vector lengths are not equal");
//   }
//   double u = Rcpp::runif(1, 0.0, 1.0)[0];
//   int n = X.size()-1;
//   Rcpp::NumericVector left(n), right(n);
//   for (int i=0; i < n; ++i) {
//     Rcpp::NumericVector current = eadelta_intervals(k,X.at(i),X.at(i+1),times.at(i),times.at(i+1),min,v);
//     left.at(i) = current.at(0);
//     right.at(i) = current.at(1);
//   }
//   // while u \in (S_{2k+1}, S_{2k}) = (left,right), we keep calculating Cauchy sums
//   double left_product = product_vector_elements(left);
//   double right_product = product_vector_elements(right);
//   while (left_product < u && u < right_product) {
//     // increase k and update the intervals
//     k += 1;
//     for (int i=0; i < n; ++i) {
//       Rcpp::NumericVector current = eadelta_intervals(k,X.at(i),X.at(i+1),times.at(i),times.at(i+1),min,v);
//       left.at(i) = current.at(0);
//       right.at(i) = current.at(1);
//     }
//     // recalculate products
//     left_product = product_vector_elements(left);
//     right_product = product_vector_elements(right);
//   }
//   if (u <= left_product) {
//     return true;
//   } else {
//     return false;
//   }
// }
