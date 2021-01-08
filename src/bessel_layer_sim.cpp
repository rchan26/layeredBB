#include "../inc/bessel_layer_sim.hpp"
#include "../inc/cauchy_sums.hpp"
#include "../inc/coins_cauchy.hpp"

using namespace Rcpp;

//' Bessel Layer simulation
//'
//' Simulates a Bessel layer l for a given sequence a
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start time of Brownian bridge
//' @param t end time of Brownian bridge
//' @param a vector/sequence of numbers
//'
//' @examples
//' bessel_layer_simulation(x = 0, y = 0, s = 0, t = 1, mult = 0.5)
//' 
//' @return 
//' A list with the following items:
//' \describe{
//'   \item{L}{Hard lower bound}
//'   \item{l}{Soft lower bound}
//'   \item{u}{Soft upper bound}
//'   \item{U}{Hard upper bound}
//' }
//' where the Bessel layer is [L, U] and either the minimum occurs in [L, l] or
//' the maximum occurs in [u, U] 
//'
//' @export
// [[Rcpp::export]]
Rcpp::List bessel_layer_simulation(const double &x,
                                   const double &y,
                                   const double &s,
                                   const double &t,
                                   const double &mult = 1)
{
  int l = 1;
  double xandy = std::min(x, y);
  double xoy = std::max(x, y);
  double layer_size = sqrt(t-s)*mult;
  while (true) {
    // flip gamma coin to determine if BB stays within the interval:
    // [min(x,y) - layer_size*l, max(x,y) + layer_size*l]
    if (gamma_coin(1, x, y, s, t, xandy - layer_size*l, xoy + layer_size*l)) {
      return List::create(Named("L", xandy - layer_size*l),
                          Named("l", xandy - layer_size*(l-1)),
                          Named("u", xoy + layer_size*(l-1)),
                          Named("U", xoy + layer_size*l));
    }
    l += 1;
  }
}

//' Multiple Bessel Layer simulation
//'
//' Simulates a Bessel layer l for a given sequence a for each component of the Brownian bridge
//'
//' @param dim dimension of Brownian bridge
//' @param x vector of start values of Brownian bridge (length of vector must be equal to dim)
//' @param y vector of end values of Brownian bridge (length of vector must be equal to dim)
//' @param s start time of Brownian bridge
//' @param t end time of Brownian bridge
//' @param a vector/sequence of numbers
//'
//' @examples
//' # simulate layer information for two-dimensional Brownian bridge starting 
//' # and ending at (0,0) in time [0,1]
//' multi_bessel_layer_simulation(dim = 2,
//'                               x = c(0, 0),
//'                               y = c(0, 0),
//'                               s = 0,
//'                               t = 1,
//'                               mult = 0.5)
//' 
//' @return 
//' A list of length dim where list[i] is the Bessel layer for component i,
//' which is represented in a list with the following items:
//' \describe{
//'   \item{L}{Hard lower bound}
//'   \item{l}{Soft lower bound}
//'   \item{u}{Soft upper bound}
//'   \item{U}{Hard upper bound}
//' }
//' where the Bessel layer for compnent i is [L, U] and either the minimum 
//' occurs in [L, l] or the maximum occurs in [u, U] 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List multi_bessel_layer_simulation(const int &dim,
                                         const arma::vec &x,
                                         const arma::vec &y,
                                         const double &s,
                                         const double &t,
                                         const double &mult = 1)
{
  // check that x and y match the dimensions of dim
  if (x.size() != dim) {
    stop("layeredBB::multi_bessel_layer_simulation: size of x is not equal to dim");
  } else if (y.size() != dim) {
    stop("layeredBB::multi_bessel_layer_simulation: size of y is not equal to dim");
  }
  // for each component, we simulate a Bessel layer and store as value in list
  Rcpp::List layers(dim);
  for (int i=0; i < dim; ++i) {
    layers.at(i) = bessel_layer_simulation(x.at(i), y.at(i), s, t, mult);
  }
  return(layers);
}
