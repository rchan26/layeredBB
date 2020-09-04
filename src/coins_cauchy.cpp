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

//' Gamma coin flipper
//'
//' Flips 'Gamma coin'; uses the Cauchy sequence S^{gamma}_{k} to 
//' determine whether or not the Brownian bridge starting at x, ending at y, between [s,t]
//' remains in interval [l,v]
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' @param k integer value
//'
//' @examples
//' gamma_coin(x = 0, y = 0, s = 0, t = 1, l = -0.5, v = 0.5, k = 1)
//'
//' @return boolean value: if T, accept probability that Brownian bridge remains in [l,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool gamma_coin(const double &x, 
                const double &y, 
                const double &s, 
                const double &t,
                const double &l, 
                const double &v,
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

//' Gamma coin flipper for intervals
//'
//' Flips 'Gamma coin' for intervals; takes the product of the Cauchy sequence S^{gamma}_{k} to 
//' determine whether or not the Brownian bridge remains in the interval [l,v]
//' Vectors x, y, s, t should all be the same length, L, where for i = 1, ..., L, the Brownian Bridge skeleton 
//' we have is broken up so that x[i] goes to y[i] between s[i] and t[i] - see example
//'
//' @param x vector of values
//' @param y vector of values
//' @param s vector of values
//' @param t vector of values
//' @param l lower bound of Brownian bridge
//' @param v upper bound of Brownian bridge
//' @param k integer value
//'
//' @examples
//' # setting up vectors
//' x_vect <- c(); y_vect <- c(); s_vect <- c(); t_vect <- c()
//' brownian_bridge <- matrix(c(0, 0, -0.2, 0.4, 1, 1), ncol = 3, nrow = 2)
//' for (i in 1:(ncol(brownian_bridge)-1)) {
//'   x_vect[i] <- brownian_bridge[1,i]
//'   y_vect[i] <- brownian_bridge[1,(i+1)]
//'   s_vect[i] <- brownian_bridge[2,i]
//'   t_vect[i] <- brownian_bridge[2,(i+1)]
//' }
//' 
//' # flip gamma coin whether or not Brownian bridge remains in [-0.5, 1.5]
//' gamma_coin_intervals(x = x_vect, y = y_vect, s = s_vect, t = t_vect, l = -0.5, v = 1.5, k = 1)
//'
//' @return boolean value: if T, accept probability that Brownian bridge remains in [l,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool gamma_coin_intervals(const Rcpp::NumericVector &x, 
                          const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, 
                          const Rcpp::NumericVector &t,
                          const double &l, 
                          const double &v, 
                          int k)
{
  // check if vector lengths are all the same
  if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
    stop("layeredBB::gamma_coin_intervals: vector lengths are not equal");
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

//' Delta coin flipper
//'
//' Flips 'Delta coin'; uses the Cauchy sequence S^{delta}_{k} to 
//' determine whether or not the Brownian bridge with minimum, min, 
//' starting at x, ending at y, between [s,t] remains in interval [l,v]
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' @param k integer value
//'
//' @examples
//'
//' @return boolean value: if T, accept probability that Brownian bridge with minimum, min, remains in [l,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool delta_coin(const double &x, 
                const double &y, 
                const double &s, 
                const double &t,
                const double &min, 
                const double &v,
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

//' Delta coin flipper for intervals
//'
//' Flips 'Delta coin' for intervals; takes the product of the Cauchy sequence S^{delta}_{k} to 
//' determine whether or not the Brownian bridge with minimum, min, remains in the interval [l,v]
//' Vectors x, y, s, t should all be the same length, L, where for i = 1, ..., L, the Brownian Bridge skeleton 
//' we have is broken up so that x[i] goes to y[i] between s[i] and t[i] - see example
//'
//' @param x vector of values
//' @param y vector of values
//' @param s vector of values
//' @param t vector of values
//' @param min minimum of Brownian bridge
//' @param v upper bound of Brownian bridge
//' @param k integer value
//'
//' @examples
//' # setting up vectors
//' x_vect <- c(); y_vect <- c(); s_vect <- c(); t_vect <- c()
//' brownian_bridge <- matrix(c(0, 0, -0.2, 0.4, 1, 1), ncol = 3, nrow = 2)
//' for (i in 1:(ncol(brownian_bridge)-1)) {
//'   x_vect[i] <- brownian_bridge[1,i]
//'   y_vect[i] <- brownian_bridge[1,(i+1)]
//'   s_vect[i] <- brownian_bridge[2,i]
//'   t_vect[i] <- brownian_bridge[2,(i+1)]
//' }
//' 
//' # flip delta coin whether or not Brownian bridge remains in [-0.2, 1.5]
//' delta_coin_intervals(x = x_vect, y = y_vect, s = s_vect, t = t_vect, min = -0.2, v = 1.5, k = 1)
//'
//' @return boolean value: if T, accept probability that Brownian bridge with minimum, min, remains in [l,v], otherwise reject
//'
//' @export
// [[Rcpp::export]]
bool delta_coin_intervals(const Rcpp::NumericVector &x, 
                          const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, 
                          const Rcpp::NumericVector &t,
                          const double &min, 
                          const double &v, 
                          int k)
{
  // check if vector lengths are all the same
  if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
    stop("layeredBB::gamma_coin_intervals: vector lengths are not equal");
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
