#include "../inc/inverse_gauss.hpp"
#include <Rcpp.h>
#include <random>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

//' Inverse Gaussian Sampler
//'
//' This function returns 1 sample from an Inverse Gaussian distribution with mean mu and shape lambda
//'
//' @param mu mean
//' @param lambda shape
//' 
//' @return real value: simulated point from Inverse Gaussian distribution with mean mu and shape lambda
//'
//' @examples
//' curve(statmod::dinvgauss(x, mean = 1, shape = 1), 0, 4)
//' samples <- sapply(1:10000, function(i) inv_gauss_sampler(mu = 1, lambda =1))
//' lines(density(x = samples, adjust = 0.5), col = 'blue')
//'
//' @export
// [[Rcpp::export]]
double inv_gauss_sampler(const double &mu, const double &lambda)
{
  // function simulates from an inverse Gaussian distribution with mean mu and scale lambda
  
  // simulate from normal and uniform distributions
  double gauss_sim = Rcpp::rnorm(1, 0.0, 1.0)[0];
  double unif_sim = Rcpp::runif(1, 0.0, 1.0)[0];
  
  double y = gauss_sim*gauss_sim;
  double x = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*y + mu*mu*y*y);
  if (unif_sim <= (mu/(mu+x))) {
    return x;
  } else {
    return (mu*mu/x);
  }
}
