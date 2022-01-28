#include "../inc/brownian_bridge.hpp"
#include "../inc/inverse_gauss.hpp"
#include "../inc/double_comparison.hpp"

using namespace Rcpp;

//' Brownian Motion path sampler 
//'
//' Simulation of a path of a Brownian motion at given times
//'
//' @param x start value of Brownian motion
//' @param times vector of real numbers to simulate Brownian motion
//'
//' @return Matrix of the simulated Brownian motion path at all
//'         included time points. The times are sorted. 
//'         The first row are the points of the Brownian motion (named 'X')
//'         second row are corresponding times (named 'times')
//'
//' @examples
//' # simulating path for Brownian motion starting at 0 between [0,1]
//' path <- Brownian_motion_path_sampler(x = 0, times = seq(0, 1, 0.01))
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X')
//' lines(x = path['time',], y = path['X',])
//' 
//' # comparing the simulated distribution of simulated points to the
//' # theoretical distribution of simulated points
//' # set variables
//' x <- 0
//' start_time <- 1.8
//' end_time <- 5
//' replicates <- 10000
//' paths <- list()
//' # repeatedly simulate Brownian bridge 
//' for (i in 1:replicates) {
//'   paths[[i]] <- Brownian_motion_path_sampler(x = x, times = seq(start_time, end_time, 0.01))
//' }
//' # select the points at the specified time q
//' index <- which(seq(start_time, end_time, 0.01)==end_time)
//' simulated_points <- sapply(1:replicates, function(i) paths[[i]]['X', index])
//' # calculate the theoretical mean and standard deviation of the simulated points at time q
//' theoretical_mean <- x
//' theoretical_sd <- sqrt(end_time-start_time)
//' # plot distribution of the simulated points and the theoretical distribution
//' plot(density(simulated_points))
//' curve(dnorm(x, theoretical_mean, theoretical_sd), add = T, col = 'red')
//' print(paste('Theoretical variance is', end_time-start_time, 'and sample variance is', var(simulated_points)))
// [[Rcpp::export]]
Rcpp::NumericMatrix Brownian_motion_path_sampler(const double &x,
                                                 const Rcpp::NumericVector &times)
{
  if (times.size() < 2) {
    stop("layeredBB::Brownian_motion_path_sampler: times must be a vector of length at least 2");
  }
  // ----- sort times
  Rcpp::NumericVector sorted_times = Rcpp::sort_unique(times);
  if (sorted_times.size() < 2) {
    stop("layeredBB::Brownian_motion_path_sampler: times must be a vector with at least 2 unique values in");
  }
  // ----- create matrix variable to store Brownian motion
  Rcpp::NumericMatrix bm(2, sorted_times.size());
  bm(0,0) = x;
  bm.row(1) = sorted_times;
  rownames(bm) = CharacterVector::create("X", "time");
  // ----- when simulating the path, want to work from left to right
  for (int i = 1; i < bm.ncol(); ++i) {
    bm(0, i) = R::rnorm(bm(0,i-1), sqrt(sorted_times[i]-sorted_times[i-1]));
  }
  return(bm);
}

//' Multi-dimensional Brownian Motion path sampler
//'
//' Simulation of a multi-dimensional Brownian Motion, at given times
//'
//' @param dim dimension of Brownian motion
//' @param x start value of Brownian motion
//' @param times vector of real numbers to simulate Brownian motion
//' 
//' @return Matrix of the simulated Brownian motion path at all
//'         included time points. The times are sorted. 
//'         The first dim rows are the points of the Brownian motion
//'         dim+1 row are corresponding times
//'
//' @examples
//' # simulate two-dimensional Brownian bridge starting and ending 
//' # at (0,0) in time [0,1]
//' multi_brownian_motion(dim = 2,
//'                       x = c(0,0),
//'                       times = c(0.1, 0.2, 0.4, 0.6, 0.8))
//'                       
//' # note that the times are sorted and duplicates are removed
//' multi_brownian_motion(dim = 2,
//'                       x = c(0.5,1),
//'                       times = c(0.1, 0.2, 0.4, 0.6, 0.6, 0.8, 0.1))
// [[Rcpp::export]]
Rcpp::NumericMatrix multi_brownian_motion(const int &dim,
                                          const Rcpp::NumericVector &x,
                                          const Rcpp::NumericVector &times)
{
  if (times.size() < 2) {
    stop("layeredBB::Brownian_motion_path_sampler: times must be a vector of length at least 2");
  }
  // ----- sort times
  Rcpp::NumericVector sorted_times = Rcpp::sort_unique(times);
  if (sorted_times.size() < 2) {
    stop("layeredBB::Brownian_motion_path_sampler: times must be a vector with at least 2 unique values in");
  }
  // ----- create matrix variable to store Brownian motion
  Rcpp::NumericMatrix bm(dim+1, sorted_times.size());
  bm.row(dim) = sorted_times;
  // loop through the components and simulate a Brownian bridge
  for (int d=0; d < dim; ++d) {
    Rcpp::NumericMatrix component_BM = Brownian_motion_path_sampler(x.at(d), times);
    bm.row(d) = component_BM.row(0);
  }
  return(bm);
}


//' Brownian Bridge path sampler (Algorithm 13 in ST329)
//'
//' Simulation of a path of a Brownian bridge at given times
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start time of Brownian bridge
//' @param t end time of Brownian bridge
//' @param times vector of real numbers to simulate Brownian bridge
//'
//' @return A list with the following components
//' \describe{
//'   \item{full_path}{Matrix of the simulated Brownian bridge path at all 
//'                    included time points, i.e. s, t and times. The times
//'                    are sorted and duplicates are removed. The first row
//'                    are the points of the Brownian bridge (named 'X') 
//'                    second row are corresponding times (named 'time')}
//'   \item{simulated_path}{Matrix of the simulated Brownian bridge path only at 
//'                         the specified times passed into the function, i.e. 
//'                         the times vector. The times are not sorted and
//'                         duplicates are not removed. The first row
//'                         are the points of the Brownian bridge (named 'X') 
//'                         second row are corresponding times (named 'time')}
//' }
//'
//' @examples 
//' # simulating paths for time [0,1] and plotting them
//' start <- runif(1, -1, 1)
//' end <- runif(1, -1, 1)
//' path <- Brownian_bridge_path_sampler(x = start,
//'                                      y = end,
//'                                      s = 0,
//'                                      t = 1,
//'                                      times = seq(0, 1, 0.01))$full_path
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X')
//' lines(x = path['time',], y = path['X',])
//' 
//' # notice that simulated_path only includes points that are included in times vector
//' # note that simulated_path does not remove duplicates passed into times
//' Brownian_bridge_path_sampler(x = 0,
//'                              y = 1,
//'                              s = 0,
//'                              t = 1,
//'                              times = c(0.1, 0.2, 0.4, 0.6, 0.6, 0.8, 0.1))
//'
//' # comparing the simulated distribution of simulated points to the
//' # theoretical distribution of simulated points
//' # set variables
//' x <- 0.53
//' y <- 4.32
//' s <- 0.53
//' t <- 2.91
//' q <- 1.72
//' replicates <- 10000
//' paths <- list()
//' # repeatedly simulate Brownian bridge 
//' for (i in 1:replicates) {
//'   paths[[i]] <- Brownian_bridge_path_sampler(x = x,
//'                                              y = y,
//'                                              s = s,
//'                                              t = t,
//'                                              times = seq(s, t, 0.01))
//' }
//' # select the points at the specified time q
//' index <- which(seq(s, t, 0.01)==q)
//' simulated_points <- sapply(1:replicates, function(i) paths[[i]]$full_path['X', index])
//' # calculate the theoretical mean and standard deviation of the simulated points at time q
//' theoretical_mean <- x + (q-s)*(y-x)/(t-s)
//' theoretical_sd <- sqrt((t-q)*(q-s)/(t-s))
//' # plot distribution of the simulated points and the theoretical distribution
//' plot(density(simulated_points))
//' curve(dnorm(x, theoretical_mean, theoretical_sd), add = T, col = 'red')
// [[Rcpp::export]]
Rcpp::List Brownian_bridge_path_sampler(const double &x,
                                        const double &y,
                                        const double &s,
                                        const double &t,
                                        const Rcpp::NumericVector &times)
{
  if (t <= s) {
    stop("layeredBB::Brownian_bridge_path_sampler: t <= s. Must have s < t");
  } if (Rcpp::min(times) < s) {
    stop("layeredBB::Brownian_bridge_path_sampler: minimum of specified times is less than s");
  } else if (Rcpp::max(times) > t) {
    stop("layeredBB::Brownian_bridge_path_sampler: maximum of specified times is greater than t");
  }
  // ----- collect all times into one vector
  // and remove duplicates and sort full_times vector
  Rcpp::NumericVector full_times = times;
  full_times.push_front(s);
  full_times.push_back(t);
  full_times = Rcpp::sort_unique(full_times);
  // ----- create two paths:
  // BB at all times included
  Rcpp::NumericMatrix full_bb(2, full_times.size());
  full_bb.row(1) = full_times;
  rownames(full_bb) = CharacterVector::create("X", "time");
  // BB at specified times
  Rcpp::NumericMatrix simulated_bb(2, times.size());
  simulated_bb.row(1) = times;
  rownames(simulated_bb) = CharacterVector::create("X", "time");
  // ----- when simulating the path, want to work from left to right
  // simulate points using the full path and fill in simulated_bb 
  for (int i = 0; i < full_times.size(); ++i) {
    // simulate path
    if (full_times.at(i) == s) {
      full_bb(0, i) = x;
    } else if (full_times.at(i) == t) {
      full_bb(0, i) = y;
    } else {
      double l = full_times.at(i-1), q = full_times.at(i), r = t;
      double W_l = full_bb(0, i-1), W_r = y;
      double M = W_l + ((q-l)*(W_r-W_l)/(r-l));
      double S = sqrt((r-q)*(q-l)/(r-l));
      full_bb(0, i) = R::rnorm(M, S);
    }
    // loop through times vector, if times.at(j) == full_times.at(i), 
    // then we know that we just simulated a point that was specified in the
    // times vector and add it to the simulated_bb path
    for (int j = 0; j < times.size(); ++j) {
      if (times.at(j) == full_times.at(i)) {
        simulated_bb(0, j) = full_bb(0, i);
      }
    }
  }
  return Rcpp::List::create(Named("full_path", full_bb),
                            Named("simulated_path", simulated_bb));
}

//' Multi-dimensional Brownian Bridge path sampler
//'
//' Simulation of a multi-dimensional Brownian Bridge, at given times
//'
//' @param dim dimension of Brownian bridge
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start time of Brownian bridge
//' @param t end time of Brownian bridge
//' @param times vector of real numbers to simulate Brownian bridge
//' 
//' @return A list with the following components
//' \describe{
//'   \item{full_path}{Matrix of the simulated Brownian bridge path at all 
//'                    included time points, i.e. s, t and times. The times
//'                    are sorted and duplicates are removed. The first dim rows
//'                    are the points of the Brownian bridge in each component,
//'                    last row gives the corresponding times}
//'   \item{simulated_path}{Matrix of the simulated Brownian bridge path only at 
//'                         the specified times passed into the function, i.e. 
//'                         the times vector. The times are not sorted and
//'                         duplicates are not removed. The first dim rows
//'                         are the points of the Brownian bridge in each component,
//'                         last row gives the corresponding times}
//' }
//'
//' @examples
//' # simulate two-dimensional Brownian bridge starting and ending 
//' # at (0,0) in time [0,1]
//' multi_brownian_bridge(dim = 2,
//'                       x = c(0,0),
//'                       y = c(0,0),
//'                       s = 0,
//'                       t = 1,
//'                       times = c(0.1, 0.2, 0.4, 0.6, 0.8))
//'                       
//' # note that simulated_path does not remove duplicates passed into times
//' multi_brownian_bridge(dim = 2,
//'                       x = c(0,0),
//'                       y = c(0,0),
//'                       s = 0,
//'                       t = 1,
//'                       times = c(0.1, 0.2, 0.4, 0.6, 0.6, 0.8, 0.1))
// [[Rcpp::export]]
Rcpp::List multi_brownian_bridge(const int &dim,
                                 const Rcpp::NumericVector &x,
                                 const Rcpp::NumericVector &y,
                                 const double &s,
                                 const double &t,
                                 const Rcpp::NumericVector &times)
{
  // check that x and y match the dimensions of dim
  if (t <= s) {
    stop("layeredBB::multi_brownian_bridge: t <= s. Must have s < t");
  } else if (x.size() != dim) {
    stop("layeredBB::multi_brownian_bridge: size of x is not equal to dim");
  } else if (y.size() != dim) {
    stop("layeredBB::multi_brownian_bridge: size of y is not equal to dim");
  } else if (Rcpp::min(times) < s) {
    stop("layeredBB::multi_brownian_bridge: minimum of specified times is less than s");
  } else if (Rcpp::max(times) > t) {
    stop("layeredBB::multi_brownian_bridge: maximum of specified times is greater than t");
  }
  // ----- collect all times into one vector
  // and remove duplicates and sort full_times vector
  Rcpp::NumericVector full_times = times;
  full_times.push_front(s);
  full_times.push_back(t);
  full_times = Rcpp::sort_unique(full_times);
  // ----- for each component, we simulate a Brownian bridge
  // BB at all times included
  Rcpp::NumericMatrix multi_full_bb(dim+1, full_times.size());
  multi_full_bb.row(dim) = full_times;
  // BB at specified times
  Rcpp::NumericMatrix multi_simulated_bb(dim+1, times.size());
  multi_simulated_bb.row(dim) = times;
  // loop through the components and simulate a Brownian bridge
  for (int d=0; d < dim; ++d) {
    Rcpp::List component_BB = Brownian_bridge_path_sampler(x.at(d), y.at(d), s, t, times);
    // fill in full path BB
    Rcpp::NumericMatrix full_path = component_BB["full_path"];
    multi_full_bb.row(d) = full_path.row(0);
    // fill in simulated BB
    Rcpp::NumericMatrix simulated_path = component_BB["simulated_path"];
    multi_simulated_bb.row(d) = simulated_path.row(0);
  }
  return Rcpp::List::create(Named("full_path", multi_full_bb),
                            Named("simulated_path", multi_simulated_bb));
}

double M_func(const double &a,
              const double &x,
              const double &y,
              const double &s,
              const double &t)
{
  // M function that is used to simulate a minimum of a Brownian bridge
  return exp(-2.0*(a-x)*(a-y)/(t-s));
}

//' Brownian Bridge minimum point sampler (Algorithm 14 in ST329)
//'
//' Simulation of a minimum point of a Brownian bridge
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start time of Brownian bridge
//' @param t end time of Brownian bridge
//' @param low_bound Lower bound of minimum point (low_bound < up_bound <= min(x,y))
//' @param up_bound Upper bound of minimum point (low_bound < up_bound <= min(x,y))
//' @param checks logical value to determine if arguments that are passed into 
//'        the function are checked. Things that are checked include that 
//'        low_bound < up_bound <= min(x,y) and that s < t
//'
//' @return vector: the simulated minimum, 'min', and time where minimum occurs, 'tau'
//'
//' @examples
//' # simulate a minimum between 0 and 1 of a Brownian bridge starting 
//' # at 0 and ending at 0 in time [0,1]
//' min_sampler(x = 0, y = 0, s = 0, t = 1, low_bound = -1, up_bound = 0)
//'
//' # plotting multiple simulated minimums and their times
//' minimums <- sapply(1:5000, function(i) {
//'   min_sampler(x = 0, y = 0, s = 0, t = 1, low_bound = -10, up_bound = 0)
//' })
//' plot(x = minimums[2,], y = minimums[1,], pch = 20, lwd = 0.1,
//'      xlab = 'Time', ylab = 'X')
// [[Rcpp::export]]
Rcpp::NumericVector min_sampler(const double &x,
                                const double &y,
                                const double &s,
                                const double &t,
                                const double &low_bound,
                                const double &up_bound,
                                const bool &checks = true)
{
  if (checks) {
    if (low_bound >= up_bound) {
      stop("layeredBB::min_sampler: low_bound >= up_bound. Must have low_bound < up_bound <= min(x,y)");
    } else if (up_bound > std::min(x, y)) {
      stop("layeredBB::min_sampler: up_bound > min(x,y). Must have low_bound < up_bound <= min(x,y)");
    } else if (t <= s) {
      stop("layeredBB::min_sampler: t <= s. Must have s < t");
    }
  }
  const double u1 = R::runif(M_func(low_bound,x,y,s,t), M_func(up_bound,x,y,s,t));
  // if u1 == 1, then m = min(x,y) as log(u1) = 0
  // when u1 very close to 1, causes numerical precision issues
  if (almostEqual(u1, 1, 1e-07)) {
    if (x <= y) {
      // if x <= y, return (x, s)
      return Rcpp::NumericVector::create(Named("min", x), Named("tau", s));
    } else if (y < x) {
      // if x > y, return (y, t)
      return Rcpp::NumericVector::create(Named("min", y), Named("tau", t));
    }
  }
  // set simulated minimum value
  const double m = x - (0.5*(sqrt((y-x)*(y-x) - 2.0*(t-s)*log(u1)) - y + x));
  // simulating from Inverse Gaussian to set V and tau
  double mu, lambda, V;
  const double condition = (x - m) / (x + y - 2.0*m);
  const double u2 = R::runif(0, 1);
  if (u2 < condition) {
    mu = (y-m)/(x-m);
    lambda = (y-m)*(y-m)/(t-s);
    V = inv_gauss_sampler(mu, lambda);
  } else {
    mu = (x-m)/(y-m);
    lambda = (x-m)*(x-m)/(t-s);
    V = 1.0 / inv_gauss_sampler(mu, lambda);
  }
  const double tau = (s*V+t)/(1.0+V);
  if (tau == s || (std::isnan(tau))) {
    // if s == tau, just return (x, m)
    // tau is NaN if V is inf
    return Rcpp::NumericVector::create(Named("min", x), Named("tau", s));
  } else if (tau == t) {
    // if t == tau, just return (y, m)
    return Rcpp::NumericVector::create(Named("min", y), Named("tau", t));
  } else {
    return Rcpp::NumericVector::create(Named("min", m), Named("tau", tau));
  }
}

//' Bessel Bridge point sampler given minimum (Algorithm 15 in ST329)
//'
//' Simulation of a point of a Bessel bridge at time q, given minimum occurs at time tau
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start time of Bessel bridge
//' @param t end time of Bessel bridge
//' @param m minimum point
//' @param tau time of minimum point
//' @param q time of simulation
//' @param checks logical value to determine if arguments that are passed into 
//'        the function are checked. Things that are checked include that
//'        s < t, that q is in [s,t], that tau is in [s,t], that m <= min(x,y) 
//'        and that if tau == s or tau == t, then m == x or m == y, respectively
//'
//' @return simulated point of the Bessel bridge at time q
//'
//' @examples
//' # simulating a point at q=0.2 for a Bessel bridge starting at 0 and ending 
//' # at 0 in time [0,1] given minimum is at -0.4 at time 0.6
//' min_Bessel_bridge_sampler(x = 0,
//'                           y = 0,
//'                           s = 0,
//'                           t = 1,
//'                           m = -0.4,
//'                           tau = 0.6,
//'                           q = 0.2)
// [[Rcpp::export]]
double min_Bessel_bridge_sampler(const double &x,
                                 const double &y,
                                 const double &s,
                                 const double &t,
                                 const double &m,
                                 const double &tau,
                                 const double &q,
                                 const bool &checks = true)
{
  if (checks) {
    // check that s < t
    if (t <= s) {
      stop("layeredBB::min_Bessel_bridge_sampler: t <= s. Must have s < t");
    }
    // check requested simulation time q is in [s,t]
    if (q < s) {
      stop("layeredBB::min_Bessel_bridge_sampler: requested simulation time q < s");
    } else if (q > t) {
      stop("layeredBB::min_Bessel_bridge_sampler: requested simulation time q > t");
    }
    // check tau is between [s,t]
    if (tau < s) {
      stop("layeredBB::min_Bessel_bridge_sampler: time of minimum tau < s");
    } else if (tau > t) {
      stop("layeredBB::min_Bessel_bridge_sampler: time of minimum tau > t");
    }
    // check m <= min(x,y)
    if (m > std::min(x,y)) {
      stop("layeredBB::min_Bessel_bridge_sampler: m > min(x,y). Must have m <= min(x,y)");
    }
    // if tau == s or tau == t
    // check that they are consistent with the given points for m, x, y
    if (tau == s) {
      if (m != x) {
        stop("layeredBB::min_Bessel_bridge_sampler: tau == s and minimum point m != x (within reasonable precision)");
      }
    } else if (tau == t) {
      if (m != y) {
        stop("layeredBB::min_Bessel_bridge_sampler: tau == t and minimum point m != y (within reasonable precision)");
      }
    }
  }
  // function simulates a Bessel bridge at a given time (q) with minimum (min) at time (tau)
  // cases where we already know the location at time q
  if (q == s) {
    return x;
  } else if (q == t) {
    return y;
  } else if (q == tau) {
    return m;
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
  const Rcpp::NumericVector b = rnorm(3, 0.0, (sqrt(fabs(tau-q) * fabs(q-r)) / fabs(tau-r)));
  // set simulated value and return
  const double term1 = ((Wr-m)*fabs(tau-q)/pow(fabs(tau-r), 1.5)) + b.at(0);
  return m + sqrt(fabs(tau-r) * (term1*term1 + b.at(1)*b.at(1) + b.at(2)*b.at(2)));
}

//' Bessel Bridge path sampler given minimum
//'
//' Simulation of a path of a Bessel bridge at given times, given minimum occurs at time tau
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start time of Bessel bridge
//' @param t end time of Bessel bridge
//' @param m minimum point
//' @param tau time of minimum point
//' @param times vector of real numbers to simulate Bessel bridge
//' @param checks logical value to determine if arguments that are passed into 
//'        the function are checked. Things that are checked include that 
//'        s < t, that requested simulation times are in [s,t], that m <= min(x,y)
//'        and that if tau == s or tau == t, then m == x or m == y, respectively
//'
//' @return A list with the following components
//' \describe{
//'   \item{full_path}{Matrix of the simulated Bessel bridge path at all 
//'                    included time points, i.e. s, t and times. The times
//'                    are sorted and duplicates are removed. The first row
//'                    are the points of the Brownian bridge (named 'X') 
//'                    second row are corresponding times (named 'time')}
//'   \item{simulated_path}{Matrix of the simulated Bessel bridge path only at 
//'                         the specified times passed into the function, i.e. 
//'                         the times vector. The times are not sorted and
//'                         duplicates are not removed. The first row
//'                         are the points of the Bessel bridge (named 'X') 
//'                         second row are corresponding times (named 'time')}
//'   \item{remove_m_path}{Matrix of the simulated Bessel bridge path only at 
//'                        all included times points excluding tau. These times
//'                        are sorted and duplicates are removed. The first row
//'                        are the points of the Bessel bridge (named 'X') 
//'                        second row are corresponding times (named 'time'). 
//'                        Note that the minimum point is included if it is 
//'                        passed into the times vector}
//' }
//'
//' @examples
//' # simulating a path at times=c(0.2, 0.4, 0.8) for a Bessel bridge starting 
//' # at 0 and ending at 0 in time [0,1] given minimum is at -0.4 at time 0.6
//' min_Bessel_bridge_path_sampler(x = 0,
//'                                y = 0,
//'                                s = 0,
//'                                t = 1,
//'                                m = -0.4,
//'                                tau = 0.6,
//'                                times = c(0.2, 0.4, 0.8))
//' 
//' # note that remove_m_path will still include the minimum if passed into times
//' # also note that simulated_path does not remove duplicates passed into times
//' min_Bessel_bridge_path_sampler(x = 0,
//'                                y = 0,
//'                                s = 0,
//'                                t = 1,
//'                                m = -0.4,
//'                                tau = 0.6,
//'                                times = c(0.2, 0.4, 0.6, 0.8, 0.6))
//' 
//' # another example
//' start <- runif(1, -1, 1)
//' end <- runif(1, -1, 1)
//' min <- min_sampler(x = start,
//'                    y = end,
//'                    s = 0,
//'                    t = 1,
//'                    low_bound = min(start, end)-0.4,
//'                    up_bound = min(start, end)-0.2)
//' path <- min_Bessel_bridge_path_sampler(x = start,
//'                                        y = end,
//'                                        s = 0,
//'                                        t = 1,
//'                                        m = min['min'],
//'                                        tau = min['tau'],
//'                                        times = seq(0, 1, 0.01))$full_path
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X')
//' lines(x = path['time',], y = path['X',])
//' points(x = min['tau'], y = min['min'], col = 'red', pch = 20)
// [[Rcpp::export]]
Rcpp::List min_Bessel_bridge_path_sampler(const double &x,
                                          const double &y,
                                          const double &s,
                                          const double &t,
                                          const double &m,
                                          const double &tau,
                                          const Rcpp::NumericVector &times,
                                          const bool &checks = true)
{
  if (checks) {
    // check that s < t
    if (t <= s) {
      stop("layeredBB::min_Bessel_bridge_path_sampler: t <= s. Must have s < t");
    }
    // check requested simulation times are in [s,t]
    if (Rcpp::min(times) < s) {
      stop("layeredBB::min_Bessel_bridge_path_sampler: minimum of specified times is less than s");
    } else if (Rcpp::max(times) > t) {
      stop("layeredBB::min_Bessel_bridge_path_sampler: maximum of specified times is greater than t");
    }
    // check tau is between [s,t]
    if (tau < s) {
      stop("layeredBB::min_Bessel_bridge_path_sampler: time of minimum tau < s");
    } else if (tau > t) {
      stop("layeredBB::min_Bessel_bridge_path_sampler: time of minimum tau > t");
    }
    // check m <= min(x,y)
    if (m > std::min(x,y)) {
      stop("layeredBB::min_Bessel_bridge_path_sampler: m > min(x,y). Must have m <= min(x,y)");
    }
    // if tau == s or tau == t
    // check that they are consistent with the given points for m, x, y
    if (tau == s) {
      if (m != x) {
        stop("layeredBB::min_Bessel_bridge_path_sampler: tau == s and minimum point m != x");
      }
    } else if (tau == t) {
      if (y != m) {
        stop("layeredBB::min_Bessel_bridge_path_sampler:: tau == t and minimum point m != y");
      }
    }
  }
  // ----- collect all times into one vector
  // and remove duplicates and sort full_times vector
  Rcpp::NumericVector full_times = times;
  full_times.push_front(s);
  full_times.push_back(t);
  full_times.push_back(tau);
  full_times = Rcpp::sort_unique(full_times);
  // ----- create variable for all times besides the time of the minimum tau
  // and remove duplicates and sort full_times vector
  Rcpp::NumericVector remove_m_times = times;
  remove_m_times.push_front(s);
  remove_m_times.push_back(t);
  remove_m_times = Rcpp::sort_unique(remove_m_times);
  // ----- create three paths:
  // BB at all times included
  Rcpp::NumericMatrix full_bb(2, full_times.size());
  full_bb.row(1) = full_times;
  rownames(full_bb) = CharacterVector::create("X", "time");
  // BB at specified times
  Rcpp::NumericMatrix simulated_bb(2, times.size());
  simulated_bb.row(1) = times;
  rownames(simulated_bb) = CharacterVector::create("X", "time");
  // BB at times minus the minimum point
  Rcpp::NumericMatrix remove_m_bb(2, remove_m_times.size());
  remove_m_bb.row(1) = remove_m_times;
  rownames(remove_m_bb) = CharacterVector::create("X", "time");
  // when simulating the path, want to work from left to right when left of the min
  // and work from right to left when right of the min
  // this is so we can use the Markov property
  for (int i = 0; i < full_times.size(); ++i) {
    if (full_times.at(i) <= tau) {
      // simulate the point at each time
      if (full_times.at(i) == s) {
        full_bb(0, i) = x;
      } else if (full_times.at(i) == tau) {
        full_bb(0, i) = m;
      } else if (full_times.at(i) == t) {
        full_bb(0, i) = y;
      }  else {
        // if left of tau, then we simulate the next point normally, 
        // starting from the previous point to the end
        full_bb(0, i) = min_Bessel_bridge_sampler(full_bb(0, i-1),
                y,
                full_times.at(i-1),
                t,
                m,
                tau,
                full_times.at(i),
                false);
      }
      // loop through times vector, if times.at(j) == full_times.at(i), 
      // then we know that we just simulated a point that was specified in the
      // times vector and add it to the simulated_bb path
      for (int j = 0; j < times.size(); ++j) {
        if (times.at(j) == full_times.at(i)) {
          simulated_bb(0, j) = full_bb(0, i);
        }
      }
      // do same for remove_m path
      for (int j = 0; j < remove_m_times.size(); ++j) {
        if (remove_m_times.at(j) == full_times.at(i)) {
          remove_m_bb(0, j) = full_bb(0, i);
        }
      }
    } else {
      break;
    }
  }
  for (int i = full_times.size()-1; i > 0; --i) {
    if (full_times.at(i) > tau) {
      // simulate the point at each time
      if (full_times.at(i) == t) {
        full_bb(0, i) = y;
      } else {
        // reflect it by taking the absolute value of the (time wanted - t)
        full_bb(0, i) = min_Bessel_bridge_sampler(full_bb(0, i+1),
                x,
                fabs(full_times.at(i+1)-t),
                fabs(s-t),
                m,
                fabs(tau-t),
                fabs(full_times.at(i)-t),
                false);
      }
      // loop through times vector, if times.at(j) == full_times.at(i), 
      // then we know that we just simulated a point that was specified in the
      // times vector and add it to the simulated_bb path
      for (int j = 0; j < times.size(); ++j) {
        if (times.at(j) == full_times.at(i)) {
          simulated_bb(0, j) = full_bb(0, i);
        }
      }
      // do same for remove_m path
      for (int j = 0; j < remove_m_times.size(); ++j) {
        if (remove_m_times.at(j) == full_times.at(i)) {
          remove_m_bb(0, j) = full_bb(0, i);
        }
      }
    } else {
      break;
    }
  }
  return Rcpp::List::create(Named("full_path", full_bb),
                            Named("simulated_path", simulated_bb),
                            Named("remove_m_path", remove_m_bb));
}

//' Brownian Bridge maximum point sampler
//'
//' Simulation of a maximum point of a Brownian bridge
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start time of Brownian bridge
//' @param t end time of Brownian bridge
//' @param low_bound Lower bound of maximum point (max(x,y) <= low_bound < up_bound)
//' @param up_bound Upper bound of maximum point (max(x,y) <= low_bound < up_bound)
//' @param checks logical value to determine if arguments that are passed into 
//'        the function are checked. Things that are checked include that 
//'        max(x,y) <= low_bound < up_bound and that s < t
//'
//' @return vector: the simulated maximum, 'max', and time where maximum occurs, 'tau'
//'
//' @examples
//' # simulate a maximum between 0 and 1 of a Brownian bridge starting at 
//' # 0 and ending at 0 in time [0,1]
//' max_sampler(x = 0, y = 0, s = 0, t = 1, low_bound = 0, up_bound = 1)
//'
//' # plotting multiple simulated maximums and their times
//' maximums <- sapply(1:5000, function(i) {
//'   max_sampler(x = 0, y = 0, s = 0, t = 1, low_bound = 0 , up_bound = 10)
//' })
//' plot(x = maximums[2,], y = maximums[1,], pch = 20, lwd = 0.1,
//'      xlab = 'Time', ylab = 'X')
// [[Rcpp::export]]
Rcpp::NumericVector max_sampler(const double &x,
                                const double &y,
                                const double &s,
                                const double &t,
                                const double &low_bound,
                                const double &up_bound,
                                const bool &checks = true)
{
  if (checks) {
    if (low_bound >= up_bound) {
      stop("layeredBB::max_sampler: low_bound >= up_bound. Must have max(x,y) <= low_bound < up_bound");
    } else if (low_bound < std::max(x,y)) {
      stop("layeredBB::max_sampler: low_bound < max(x,y). Must have max(x,y) <= low_bound < up_bound");
    } else if (t <= s) {
      stop("layeredBB::max_sampler: t <= s. Must have s < t");
    }
  }
  // reflect the problem to simulate a minimum
  Rcpp::NumericVector sim_min = min_sampler(-x, -y, s, t, -up_bound, -low_bound, false);
  // reflect on x-axis
  return Rcpp::NumericVector::create(Named("max", -sim_min["min"]), 
                                     Named("tau", sim_min["tau"]));
}

//' Bessel Bridge point sampler given maximum
//'
//' Simulation of a point of a Bessel bridge at time q, given maximum occurs at time tuu
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start time of Bessel bridge
//' @param t end time of Bessel bridge
//' @param m maximum point 
//' @param tau time of maximum point
//' @param q time of simulation
//' @param checks logical value to determine if arguments that are passed into 
//'        the function are checked. Things that are checked include that
//'        s < t, that q is in [s,t], that tau is in [s,t], that m >= min(x,y) 
//'        and that if tau == s or tau == t, then m == x or m == y, respectively
//' 
//' @return simulated point of the Bessel bridge at time q
//'
//' @examples
//' # simulating a point at q=0.2 for a Bessel bridge starting at 0 and ending 
//' # at 0 in time [0,1] given maximum is at 0.4 at time 0.6
//' max_Bessel_bridge_sampler(x = 0,
//'                           y = 0,
//'                           s = 0,
//'                           t = 1,
//'                           m = 0.4,
//'                           tau = 0.6,
//'                           q = 0.2)
// [[Rcpp::export]]
double max_Bessel_bridge_sampler(const double &x,
                                 const double &y,
                                 const double &s,
                                 const double &t,
                                 const double &m,
                                 const double &tau,
                                 const double &q,
                                 const bool &checks = true)
{
  if (checks) {
    // check that s < t
    if (t <= s) {
      stop("layeredBB::max_Bessel_bridge_sampler: t <= s. Must have s < t");
    }
    // check requested simulation time q is in [s,t]
    if (q < s) {
      stop("layeredBB::max_Bessel_bridge_sampler: requested simulation time q < s");
    } else if (q > t) {
      stop("layeredBB::max_Bessel_bridge_sampler: requested simulation time q > t");
    }
    // check tau is between [s,t]
    if (tau < s) {
      stop("layeredBB::max_Bessel_bridge_sampler: time of maximum tau < s");
    } else if (tau > t) {
      stop("layeredBB::max_Bessel_bridge_sampler: time of maximum tau > t");
    }
    // check m >= max(x,y)
    if (m < std::max(x,y)) {
      stop("layeredBB::max_Bessel_bridge_sampler: m < max(x,y). Must have m >= max(x,y)");
    }
    // if tau == s or tau == t
    // check that they are consistent with the given points for m, x, y
    if (tau == s) {
      if (m != x) {
        stop("layeredBB::max_Bessel_bridge_sampler: tau == s and maximum point m != x");
      }
    } else if (tau == t) {
      if (m != y) {
        stop("layeredBB::max_Bessel_bridge_sampler:: tau == t and maximum point m != y");
      }
    }
  }
  // reflect the problem to simulate a Bessel bridge with a given minimum point
  // reflect on x-axis
  return -min_Bessel_bridge_sampler(-x, -y, s, t, -m, tau, q, false);
}

//' Bessel Bridge path sampler given maximum
//'
//' Simulation of a path of a Bessel bridge at given times, given maximum occurs at time tuu
//'
//' @param x start value of Bessel bridge
//' @param y end value of Bessel bridge
//' @param s start time of Bessel bridge
//' @param t end time of Bessel bridge
//' @param m maximum point 
//' @param tau time of maximum point
//' @param times vector of real numbers to simulate Bessel bridge
//' @param checks logical value to determine if arguments that are passed into 
//'        the function are checked. Things that are checked include that 
//'        s < t, that requested simulation times are in [s,t], that m >= max(x,y)
//'        and that if tau == s or tau == t, then m == x or m == y, respectivelys
//'
//' @return A list with the following components
//' \describe{
//'   \item{full_path}{Matrix of the simulated Bessel bridge path at all 
//'                    included time points, i.e. s, t and times. The times
//'                    are sorted and duplicates are removed. The first row
//'                    are the points of the Brownian bridge (named 'X') 
//'                    second row are corresponding times (named 'time')}
//'   \item{simulated_path}{Matrix of the simulated Bessel bridge path only at 
//'                         the specified times passed into the function, i.e. 
//'                         the times vector. The times are not sorted and
//'                         duplicates are not removed. The first row
//'                         are the points of the Bessel bridge (named 'X') 
//'                         second row are corresponding times (named 'time')}
//'   \item{remove_m_path}{Matrix of the simulated Bessel bridge path only at 
//'                        all included times points excluding tau. These times
//'                        are sorted and duplicates are removed. The first row
//'                        are the points of the Bessel bridge (named 'X') 
//'                        second row are corresponding times (named 'time'). 
//'                        Note that the maximum point is included if it is 
//'                        passed into the times vector}
//' }
//'
//' @examples
//' # simulating a path at times = c(0.2, 0.4, 0.8) for a Bessel bridge starting 
//' # at 0 and ending at 0 in time [0,1] given maximum is at 0.4 at time 0.6
//' max_Bessel_bridge_path_sampler(x = 0,
//'                                y = 0,
//'                                s = 0,
//'                                t = 1,
//'                                m = 0.4,
//'                                tau = 0.6,
//'                                times = c(0.2, 0.4, 0.8))
//'
//' # note that remove_m_path will still include the minimum if passed into times
//' # also note that simulated_path does not remove duplicates passed into times
//' max_Bessel_bridge_path_sampler(x = 0,
//'                                y = 0,
//'                                s = 0,
//'                                t = 1,
//'                                m = 0.4,
//'                                tau = 0.6,
//'                                times = c(0.2, 0.4, 0.6, 0.8, 0.6))
//'
//' # another example
//' start <- runif(1, -1, 1)
//' end <- runif(1, -1, 1)
//' max <- max_sampler(x = start,
//'                    y = end,
//'                    s = 0,
//'                    t = 1,
//'                    low_bound = max(start, end)+0.2,
//'                    up_bound = max(start, end)+0.4)
//' path <- max_Bessel_bridge_path_sampler(x = start,
//'                                        y = end,
//'                                        s = 0,
//'                                        t = 1,
//'                                        m = max['max'],
//'                                        tau = max['tau'],
//'                                        times = seq(0, 1, 0.01))$full_path
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X')
//' lines(x = path['time',], y = path['X',])
//' points(x = max['tau'], y = max['max'], col = 'red', pch = 20)
// [[Rcpp::export]]
Rcpp::List max_Bessel_bridge_path_sampler(const double &x,
                                          const double &y,
                                          const double &s,
                                          const double &t,
                                          const double &m,
                                          const double &tau,
                                          const Rcpp::NumericVector &times,
                                          const bool &checks = true)
{
  if (checks) {
    // check that s < t
    if (t <= s) {
      stop("layeredBB::max_Bessel_bridge_path_sampler: t <= s. Must have s < t");
    }
    // check requested simulation times are in [s,t]
    if (Rcpp::min(times) < s) {
      stop("layeredBB::max_Bessel_bridge_path_sampler: minimum of specified times is less than s");
    } else if (Rcpp::max(times) > t) {
      stop("layeredBB::max_Bessel_bridge_path_sampler: maximum of specified times is greater than t");
    }
    // check tau is between [s,t]
    if (tau < s) {
      stop("layeredBB::max_Bessel_bridge_path_sampler: time of maximum tau < s");
    } else if (tau > t) {
      stop("layeredBB::max_Bessel_bridge_path_sampler: time of maximum tau > t");
    }
    // check m >= max(x,y)
    if (m < std::max(x,y)) {
      stop("layeredBB::max_Bessel_bridge_path_sampler: m < max(x,y). Must have m >= max(x,y)");
    }
    // if tau == s or tau == t
    // check that they are consistent with the given points for m, x, y
    if (tau == s) {
      if (m != x) {
        stop("layeredBB::max_Bessel_bridge_path_sampler: tau == s and maximum point m != x");
      }
    } else if (tau == t) {
      if (m != y) {
        stop("layeredBB::max_Bessel_bridge_path_sampler:: tau == t and maximum point m != y");
      }
    }
  }
  // reflect the problem to simulate a Bessel bright path with a given minimum point
  Rcpp::List minimim_BB = min_Bessel_bridge_path_sampler(-x, -y, s, t, -m, tau, times, false);
  Rcpp::NumericMatrix full_path = minimim_BB["full_path"];
  Rcpp::NumericMatrix simulated_path = minimim_BB["simulated_path"];
  Rcpp::NumericMatrix remove_m_path = minimim_BB["remove_m_path"];
  // reflect on x-axis for full_path
  for (int i=0; i < full_path.ncol(); ++i) {
    full_path(0, i) = -full_path(0, i);
  }
  // reflect on x-axis for simulated_path
  for (int i=0; i < simulated_path.ncol(); ++i) {
    simulated_path(0, i) = -simulated_path(0, i);
  }
  // reflect on x-axis for remove_m_path
  for (int i=0; i < remove_m_path.ncol(); ++i) {
    remove_m_path(0, i) = -remove_m_path(0, i);
  }
  return Rcpp::List::create(Named("full_path", full_path),
                            Named("simulated_path", simulated_path),
                            Named("remove_m_path", remove_m_path));
}
