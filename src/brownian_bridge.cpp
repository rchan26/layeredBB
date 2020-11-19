#include "../inc/brownian_bridge_min.hpp"
#include "../inc/inverse_gauss.hpp"

using namespace Rcpp;

//' Brownian Bridge path sampler
//'
//' This function simulates a path of a Brownian bridge at given times
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param times vector of real numbers to simulate Brownian bridge
//'
//' @return matrix of the simulated Brownian bridge path, first row is points X, second row are corresponding times
//'
//' @example 
//' # simulate a Brownian bridge path starting at 0 and ending at 0 in time [0,1]
//' Brownian_bridge_path_sampler(0, 0, 0, 1, seq(0, 1, 0.01))
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix Brownian_bridge_path_sampler(const double &x,
                                                 const double &y,
                                                 const double &s, 
                                                 const double &t,
                                                 Rcpp::NumericVector times)
{
  // collect all times into one vector
  times.insert(times.end(), s);
  times.insert(times.end(), t);
  // sort the vector 'times' forward in time
  times.sort();
  // delete any duplicates
  times.erase(std::unique(times.begin(), times.end()), times.end());
  
  // create vector to store the simulated Bessel bridge path
  Rcpp::NumericVector simulated_bb(times.size());
  
  // when simulating the path, want to work from left to right
  for (int i = 0; i < times.size(); ++i) {
    // simulate the point at each time
    if (times.at(i) == s) {
      simulated_bb.at(i) = x;
    } else if (times.at(i) == t) {
      simulated_bb.at(i) = y;
    } else {
      double l = times.at(i-1), q = times.at(i), r = t;
      double W_l = simulated_bb.at(i-1), W_r = y;
      double M = W_l + ((q-l)*(W_r-W_l)/(r-l));
      double S = sqrt((r-q)*(q-l)/(r-l));
      simulated_bb.at(i) = Rcpp::rnorm(1, M, S)[0];
    }
  }
  
  // creating matrix to store the path and the corresponding times
  Rcpp::NumericMatrix bb(2, simulated_bb.size());
  bb(0, _) = simulated_bb;
  bb(1, _) = times;
  
  // setting rownames
  rownames(bb) = CharacterVector::create("X", "time");
  
  return bb;
}

//' Multi-dimensional Brownian Bridge path sampler
//'
//' This function simulates a multi-dimensional Brownian Bridge, at given times
//'
//' @param dim dimension of Brownian bridge
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param times vector of real numbers to simulate Brownian bridge
//' 
//' @return matrix of the simulated layered Brownian bridge path, first dim rows are points for X in each component, 
//'         last row are corresponding times
//'
//' @examples
//' # simulate two-dimensional Brownian bridge starting and ending at (0,0) in time [0,1]
//' multi_brownian_bridge(dim = 2, x = c(0,0), y = c(0,0), s = 0, t = 1, times = seq(0.2, 0.8, 0.2))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix multi_brownian_bridge(const int &dim,
                                          const Rcpp::NumericVector &x,
                                          const Rcpp::NumericVector &y,
                                          const double &s,
                                          const double &t,
                                          Rcpp::NumericVector times) {
  // check that x and y match the dimensions of dim
  if (x.size() != dim) {
    stop("multi_brownian_bridge: size of x is not equal to dim");
  } else if (y.size() != dim) {
    stop("multi_brownian_bridge: size of y is not equal to dim");
  } 
  
  // collect all times into one vector
  times.insert(times.end(), s);
  times.insert(times.end(), t);
  // sort the vector 'times' forward in time
  times.sort();
  // delete any duplicates
  times.erase(std::unique(times.begin(), times.end()), times.end());
  
  // for each component, we simulate a Brownian bridge
  // multi_BB is a matrix with dimensions (dim+1) x times.size()
  Rcpp::NumericMatrix multi_BB(dim+1, times.size());
  multi_BB(dim, _) = times;
  
  // loop through the components and simulate a Brownian bridge
  for (int i=0; i < dim; ++i) {
    Rcpp::NumericMatrix component_BB = Brownian_bridge_path_sampler(x.at(i), y.at(i), s, t, times);
    multi_BB(i, _) = component_BB(0, _);
  }
  
  return(multi_BB);
}

double M_func(const double &a,
              const double &x, 
              const double &y, 
              const double &s, 
              const double &t) {
  // M function that is used to simulate a minimum of a Brownian bridge
  return exp(-2.0*(a-x)*(a-y)/(t-s));
}

//' Brownian Bridge minimum point sampler
//'
//' This function simulates a minimum point of a Brownian bridge
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param low_bound Lower bound of minimum point (low_bound < up_bound <= min(x,y))
//' @param up_bound Upper bound of minimum point (low_bound < up_bound <= min(x,y))
//'
//' @return vector: the simulated minimum, 'min', and time where minimum occurs, 'tau'
//'
//' @examples
//' # simulate a minimum between 0 and 1 of a Brownian bridge starting at 0 and ending at 0 in time [0,1]
//' min_sampler(x=0, y=0, s=0, t=1, low_bound = -1, up_bound = 0)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector min_sampler(const double &x, 
                                const double &y,
                                const double &s, 
                                const double &t,
                                const double &low_bound, 
                                const double &up_bound)
{
  // function simulates a minimum of a Brownian bridge between (low_bound) and (up_bound)
  // first element returned is the simulated minimum
  // second element returned is the simulated time which the minimum occurs
  if (low_bound > up_bound) {
    stop("layeredBB::min_sampler: low_bound > up_bound");
  } else if (up_bound > std::min(x, y)) {
    stop("layeredBB::min_sampler: up_bound > min(x, y)");
  }
  
  // set simulated minimum value
  double min = x - (0.5*(sqrt((y-x)*(y-x)-2.0*(t-s)*log(Rcpp::runif(1, M_func(low_bound,x,y,s,t), M_func(up_bound,x,y,s,t))[0])) - y + x));
  
  // simulating from Inverse Gaussian to set V and tau
  double mu, lambda, V;
  if (Rcpp::runif(1, 0.0, 1.0)[0] < (x-min)/(x+y-(2.0*min))) {
    mu = (y-min)/(x-min);
    lambda = (y-min)*(y-min)/(t-s);
    V = inv_gauss_sampler(mu, lambda);
  } else {
    mu = (x-min)/(y-min);
    lambda = (x-min)*(x-min)/(t-s);
    V = (1.0 / inv_gauss_sampler(mu, lambda));
  }
  
  // // setting simulated minimum and tau in array
  return Rcpp::NumericVector::create(Named("min", min), 
                                     Named("tau", ((s*V)+t)/(1.0+V)));
}

//' Bessel Bridge point sampler given minimum
//'
//' This function simulates a point of a Bessel bridge at time q, given minimum occurs at time tuu
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start value of Bessel bridge
//' @param t end value of Bessel bridge
//' @param min minumum point
//' @param tau time of minimum point
//' @param q time of simulation
//'
//' @return simulated point of the Bessel bridge at time q
//'
//' @examples
//' # simulating a point at q=0.2 for a Bessel bridge starting at 0 and ending at 0 in time [0,1] given minimum is at -0.4 at time 0.6
//' min_Bessel_bridge_sampler(x = 0, y = 0, s = 0, t = 1, min = -0.4, tau = 0.6, q = 0.2)
//'
//' @export
// [[Rcpp::export]]
double min_Bessel_bridge_sampler(const double &x, 
                                 const double &y,
                                 const double &s, 
                                 const double &t,
                                 const double &min,
                                 const double &tau,
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
  
  // set variable r
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
  
  // set simulated value and return
  double term1 = ((Wr-min)*fabs(tau-q)/(pow(fabs(tau-r), 1.5))) + b.at(0);
  return (min + sqrt(fabs(tau-r)*(term1*term1 + b.at(1)*b.at(1) + b.at(2)*b.at(2))));
}

//' Bessel Bridge path sampler given minimum
//'
//' This function simulates a path of a Bessel bridge at given times, given minimum occurs at time tau
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start value of Bessel bridge
//' @param t end value of Bessel bridge
//' @param min minumum point
//' @param tau time of minimum point
//' @param times vector of real numbers to simulate Bessel bridge
//' @param keep_min if TRUE (default), the minimum point is returned in the sample path
//'
//' @return matrix of the simulated Bessel bridge path, first row is points X, second row are corresponding times
//'
//' @examples
//' # simulating a path at times=c(0.2, 0.4, 0.8) for a Bessel bridge starting at 0 and ending at 0 in time [0,1] given minimum is at -0.4 at time 0.6
//' min_Bessel_bridge_path_sampler(x = 0, y = 0, s = 0, t = 1, min = -0.4, tau = 0.6, times = c(0.2, 0.4, 0.8))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix min_Bessel_bridge_path_sampler(const double &x, 
                                                   const double &y,
                                                   const double &s, 
                                                   const double &t,
                                                   const double &min,
                                                   const double &tau,
                                                   Rcpp::NumericVector times,
                                                   const bool &keep_min = true)
{
  // function simulates a Bessel bridge sample path with minimum (min) at time (tau)
  // (times) get altered in this function to match the indices of the simulated path
  // i.e. after using the function, simulated_bb[i] is the path at times[i]
  // this is because (times) may not include the times (s), (t), (tau)
  
  // collect all times into one vector
  times.insert(times.end(), s);
  times.insert(times.end(), t);
  // might want to not keep the simulated minimum
  if (keep_min) {
    times.insert(times.end(), tau);
  }
  
  // sort the vector 'times' forward in time
  times.sort();
  // delete any duplicates
  times.erase(std::unique(times.begin(), times.end()), times.end());
  
  // create vector to store the simulated Bessel bridge path
  Rcpp::NumericVector simulated_bb(times.size());

  // when simulating the path, want to work from left to right when left of the min
  // and work from right to left when right of the min
  // this is so we can use the Markov property

  for (int i = 0; times.at(i) <= tau; ++i) {
    // simulate the point at each time
    if (times.at(i) == s) {
      simulated_bb.at(i) = x;
    } else if (times.at(i) == tau) {
      simulated_bb.at(i) = min;
    } else {
      // if left of tau, then we simulate the next point normally, starting from the previous point to the end
      simulated_bb.at(i) = min_Bessel_bridge_sampler(simulated_bb.at(i-1), y, times.at(i-1), t, min, tau, times.at(i));
    }
  }
  
  for (int i = times.size()-1; times.at(i) > tau; --i) {
    // simulate the point at each time
    if (times.at(i) == t) {
      simulated_bb.at(i) = y;
    } else {
      // reflect it by taking the absolute value of the (time wanted - t)
      simulated_bb.at(i) = min_Bessel_bridge_sampler(simulated_bb.at(i+1), x, fabs(times.at(i+1)-t), fabs(s-t), 
                                                     min, fabs(tau-t), fabs(times.at(i)-t));
    }
  }

  // creating matrix to store the path and the corresponding times
  Rcpp::NumericMatrix bb(2, simulated_bb.size());
  bb(0, _) = simulated_bb;
  bb(1, _) = times;
  
  // setting rownames
  rownames(bb) = CharacterVector::create("X", "time");
  
  return bb;
}

//' Brownian Bridge maximum point sampler
//'
//' This function simulates a maximum point of a Brownian bridge
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param low_bound Lower bound of maximum point (max(x,y) <= low_bound < up_bound)
//' @param up_bound Upper bound of maximum point (max(x,y) <= low_bound < up_bound)
//'
//' @return vector: the simulated maximum, 'max', and time where maximum occurs, 'tau'
//'
//' @examples
//' # simulate a maximum between 0 and 1 of a Brownian bridge starting at 0 and ending at 0 in time [0,1]
//' max_sampler(x=0, y=0, s=0, t=1, low_bound = 0, up_bound = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector max_sampler(const double &x, 
                                const double &y,
                                const double &s, 
                                const double &t,
                                const double &low_bound, 
                                const double &up_bound)
{
  // function simulates a maximum of a Brownian bridge between (low_bound) and (up_bound)
  // first element returned is the simulated maximum
  // second element returned is the simulated time which the maximum occurs
  
  // reflect the problem to simulate a minimum
  Rcpp::NumericVector sim_min = min_sampler(-x, -y, s, t, -up_bound, -low_bound);
  
  // reflect on x-axis
  return Rcpp::NumericVector::create(Named("max", -sim_min["min"]), 
                                     Named("tau", sim_min["tau"]));
}

//' Bessel Bridge point sampler given maximum
//'
//' This function simulates a point of a Bessel bridge at time q, given maximum occurs at time tuu
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start value of Bessel bridge
//' @param t end value of Bessel bridge
//' @param max maxumum point 
//' @param tau time of maximum point
//' @param q time of simulation
//' 
//' @return simulated point of the Bessel bridge at time q
//'
//' @examples
//' # simulating a point at q=0.2 for a Bessel bridge starting at 0 and ending at 0 in time [0,1] given maximum is at 0.4 at time 0.6
//' max_Bessel_bridge_sampler(x = 0, y = 0, s = 0, t = 1, max = 0.4, tau = 0.6, q = 0.2)
//'
//' @export
// [[Rcpp::export]]
double max_Bessel_bridge_sampler(const double &x, 
                                 const double &y,
                                 const double &s, 
                                 const double &t,
                                 const double &max,
                                 const double &tau,
                                 const double &q)
{
  // function simulates a Bessel bridge at a given time (q) with minimum (min) at time (tau)
  
  // reflect the problem to simulate a Bessel bridge with a given minimum point
  // reflect on x-axis
  return -min_Bessel_bridge_sampler(-x, -y, s, t, -max, tau, q);
}

//' Bessel Bridge path sampler given maximum
//'
//' This function simulates a path of a Bessel bridge at given times, given maximum occurs at time tuu
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start value of Bessel bridge
//' @param t end value of Bessel bridge
//' @param max maxumum point 
//' @param tau time of maximum point
//' @param times vector of real numbers to simulate Bessel bridge
//' @param keep_max if TRUE (default), the maximum point is returned in the sample path
//' 
//' @return matrix of the simulated Bessel bridge path, first row is points X, second row are corresponding times
//'
//' @examples
//' # simulating a path at times=c(0.2, 0.4, 0.8) for a Bessel bridge starting at 0 and ending at 0 in time [0,1] given maximum is at 0.4 at time 0.6
//' max_Bessel_bridge_path_sampler(x = 0, y = 0, s = 0, t = 1, max = 0.4, tau = 0.6, times = c(0.2, 0.4, 0.8))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix max_Bessel_bridge_path_sampler(const double &x,
                                                   const double &y,
                                                   const double &s, 
                                                   const double &t,
                                                   const double &max, 
                                                   const double &tau,
                                                   Rcpp::NumericVector times,
                                                   const bool &keep_max = true)
{
  // function simulates a Bessel bridge sample path with maximum (max) at time (tau)
  // (times) get altered in this function to match the indices of the simulated path
  // i.e. after using the function, simulated_bb[i] is the path at times[i]
  // this is because (times) may not include the times (s), (t), (tau)
  
  // reflect the problem to simulate a Bessel bright path with a given minimum point
  Rcpp::NumericMatrix sim_path = min_Bessel_bridge_path_sampler(-x, -y, s, t, -max, tau, times, keep_max);
  
  // reflect on x-axis
  for (int i=0; i < sim_path.ncol(); ++i) {
    sim_path(0, i) = -sim_path(0, i);
  }
  
  return sim_path;
}

