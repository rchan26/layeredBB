#include <Rcpp.h>
#include "../inc/brownian_bridge_min.hpp"
#include "../inc/inverse_gauss.hpp"

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
double M_function(const double &a, const double &x, const double &y, const double &s, const double &t) {
  // M function that is used to simulate a minimum of a Brownian bridge
  return exp(-2.0 * (a-x) * (a-y) / (t-s));
}

// [[Rcpp::export]]
Rcpp::NumericVector min_sampler(const double &x, const double &y,
                          const double &s, const double &t,
                          const double &low_bound, const double &up_bound)
{
  // function simulates a minimum of a Brownian bridge between (low_bound) and (up_bound)
  // first element returned is the simulated minimum
  // second element returned is the simulated time which the minimum occurs

  // calculate bounds for u1
  double low_M = M_function(low_bound, x, y, s, t);
  double up_M = M_function(up_bound, x, y, s, t);

  // simulate uniform random variables
  double u1 = Rcpp::runif(1, low_M, up_M)[0];
  double u2 = Rcpp::runif(1, 0.0, 1.0)[0];

  // set simulated minimum value
  double min = x - (0.5*(sqrt((y-x)*(y-x) - 2.0*(t-s)*log(u1)) - y + x));

  // condition for setting V
  double condition = (x-min)/(x+y-(2.0*min));
  // simulating from Inverse Gaussian
  double mu, lambda, V;
  if (u2 < condition) {
    mu = (y-min)/(x-min);
    lambda = (y-min)*(y-min)/(t-s);
    V = inv_gauss_sampler(mu, lambda);
  } else {
    mu = (x-min)/(y-min);
    lambda = (x-min)*(x-min)/(t-s);
    V = 1.0 / inv_gauss_sampler(mu, lambda);
  }

  // set tau (time of simualted minimum)
  double tau = ((s*V)+t)/(1.0+V);
  // setting simulated minimum and tau in array
  Rcpp::NumericVector simulated_min = Rcpp::NumericVector::create(min, tau);

  return simulated_min;
}

// [[Rcpp::export]]
double min_Bessel_bridge_sampler(const double &x, const double &y,
                                 const double &s, const double &t,
                                 const double &min, const double &tau,
                                 const double &q)
{
  // function simulates a Bessel bridge at a given time (q) with minimum (min) at time (tau)

  // cases where we already know the location at time q
  if (q == s) {
    return x;
  } else if (q == t) {
    return y;
  } else if (q == tau) {
    return min;
  }

  double r;
  double Wr;
  if (q < tau) {
    r = s;
    Wr = x;
  } else {
    r = t;
    Wr = y;
  }

  // simulate normal random variables
  double std_dev = sqrt(fabs(tau-q)*fabs(q-r))/fabs(tau-r);
  Rcpp::NumericVector b(3);
  for (int i=0; i<=2; ++i) {
    b[i] = rnorm(1, 0.0, std_dev)[0];
  }

  // set simulated value
  double term1 = ((Wr-min)*fabs(tau-q)/(pow(fabs(tau-r), 1.5))) + b[0];
  double W = min + sqrt(fabs(tau-r)*(term1*term1 + b[1]*b[1] + b[2]*b[2]));

  return W;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix min_Bessel_bridge_path_sampler(const double &x, const double &y,
                                                   const double &s, const double &t,
                                                   const double &min, const double &tau,
                                                   Rcpp::NumericVector &times)
{
  // function simulates a Bessel bridge sample path with minimum (min) at time (tau)
  // (times) get altered in this function to match the indices of the simulated path
  // i.e. after using the function, simulated_bb[i] is the path at times[i]
  // this is because (times) may not include the times (s), (t), (tau)

  // collect all times into one vector
  times.insert(times.end(), s);
  times.insert(times.end(), t);
  times.insert(times.end(), tau);
  // sort the vector 'time' forward in time
  times.sort();
  // delete any duplicates
  times.erase(std::unique(times.begin(), times.end()), times.end());
  
  // create vector to store the simulated Bessel bridge path
  NumericVector simulated_bb(times.size());
  for (int i = 0; i < times.size(); ++i) {
    // simulate the point at each time
    simulated_bb[i] = min_Bessel_bridge_sampler(x, y, s, t, min, tau, times[i]);
  }
  
  Rcpp::NumericMatrix bb(2, simulated_bb.size());
  bb(0, _) = simulated_bb;
  bb(1, _) = times;
  
  return bb;
}

// [[Rcpp::export]]
Rcpp::NumericVector max_sampler(const double &x, const double &y,
                                const double &s, const double &t,
                                const double &low_bound, const double &up_bound)
{
  // function simulates a maximum of a Brownian bridge between (low_bound) and (up_bound)
  // first element returned is the simulated maximum
  // second element returned is the simulated time which the maximum occurs

  // reflect the problem to simulate a minimum
  Rcpp::NumericVector sim_min = min_sampler(-x, -y, s, t, -up_bound, -low_bound);

  // reflect on x-axis
  return Rcpp::NumericVector::create(-sim_min[0], sim_min[1]);
}

// [[Rcpp::export]]
double max_Bessel_bridge_sampler(const double &x, const double &y,
                                 const double &s, const double &t,
                                 const double &max, const double &tau,
                                 const double &q)
{
  // function simulates a Bessel bridge at a given time (q) with minimum (min) at time (tau)

  // reflect the problem to simulate a Bessel bridge with a given minimum point
  double W = min_Bessel_bridge_sampler(-x, -y, s, t, -max, tau, q);

  // reflect on x-axis
  return -W;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix max_Bessel_bridge_path_sampler(const double &x, const double &y,
                                                   const double &s, const double &t,
                                                   const double &max, const double &tau,
                                                   Rcpp::NumericVector &times)
{
  // function simulates a Bessel bridge sample path with maximum (max) at time (tau)
  // (times) get altered in this function to match the indices of the simulated path
  // i.e. after using the function, simulated_bb[i] is the path at times[i]
  // this is because (times) may not include the times (s), (t), (tau)

  // reflect the problem to simulate a Bessel bright path with a given minimum point
  Rcpp::NumericMatrix sim_path = min_Bessel_bridge_path_sampler(-x, -y, s, t, -max, tau, times);

  // reflect on x-axis
  for (int i=0; i < sim_path.ncol(); ++i) {
    sim_path(0, i) = -sim_path(0, i);
  }
  return sim_path;
}

