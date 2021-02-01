#include "../inc/layered_brownian_bridge.hpp"
#include "../inc/brownian_bridge.hpp"
#include "../inc/inverse_gauss.hpp"
#include "../inc/cauchy_sums.hpp"
#include "../inc/coins_cauchy.hpp"
#include "../inc/bessel_layer_sim.hpp"

using namespace Rcpp;

//' Layered Brownian Bridge sampler (Algorithm 33 in ST329)
//'
//' Simulation of a layered Brownian Bridge given a Bessel layer at user specified times
//'
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param a vector/sequence of numbers
//' @param l integer number denoting Bessel layer, i.e. Brownian bridge 
//'        is contained in [min(x,y)-a[l], max(x,y)+a[l]]
//' @param times vector of real numbers to simulate Bessel bridge
//' 
//' @return A list with the following components
//' \describe{
//'   \item{full_path}{Matrix of the simulated layered Brownian bridge path at 
//'                    all included time points, i.e. s, t and times. The times
//'                    are sorted and duplicates are removed. The first row
//'                    are the points of the Brownian bridge (named 'X') 
//'                    second row are corresponding times (named 'times')}
//'   \item{simulated_path}{Matrix of the simulated layered Brownian bridge path 
//'                         only at the specified times passed into the function, 
//'                         i.e. the times vector. The times are not sorted and 
//'                         duplicates are not removed. The first row 
//'                         are the points of the layered Brownian bridge (named 
//'                         'X') second row are corresponding times (named 
//'                         'times')}
//'   \item{remove_m_path}{Matrix of the simulated layered Brownian bridge path 
//'                        only at all included times points excluding tau. These 
//'                        times are sorted and duplicates are removed. The first
//'                        row are the points of the layered Brownian bridge 
//'                        (named 'X') second row are corresponding times 
//'                        (named 'times')}
//' }
//'
//' @examples
//' # simulate Bessel layer
//' bes_layer <- bessel_layer_simulation(x = 0,
//'                                      y = 0,
//'                                      s = 0,
//'                                      t = 1,
//'                                      mult = 0.2)
//' # simulate layered Brownian bridge
//' # notice full_path has all included times and are sorted and have no duplicates
//' # simulated_path only returns points that are passed into times
//' # remove_m does not include the simulated minimum or maximum point
//' layered_brownian_bridge(x = 0,
//'                         y = 0,
//'                         s = 0,
//'                         t = 1,
//'                         bessel_layer = bes_layer,
//'                         times = c(0.2, 0.4, 0.6, 0.8))
//'
//' # note that simulated_path does not remove duplicates passed into times
//' layered_brownian_bridge(x = 0,
//'                         y = 0,
//'                         s = 0,
//'                         t = 1,
//'                         bessel_layer = bes_layer,
//'                         times = c(0.2, 0.4, 0.6, 0.8, 0.4, 0.6))
//'
//' # another example
//' start <- runif(1, -1, 1)
//' end <- runif(1, -1, 1)
//' bes_layer <- bessel_layer_simulation(x = start, y = end, s = 0, t = 1, mult = 0.2)
//' path <- layered_brownian_bridge(x = start,
//'                                 y = end,
//'                                 s = 0,
//'                                 t = 1,
//'                                 bessel_layer = bes_layer,
//'                                 times = seq(0, 1, 0.01))$full_path
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X',
//'      ylim = c(bes_layer$L, bes_layer$U))
//' lines(x = path['time',], y = path['X',])
//' abline(h=c(bes_layer$L, bes_layer$U), col = 'red')
//' abline(h=c(bes_layer$l, bes_layer$u), col = 'red', lty = 2)
//' 
//' # compare the simulated distribution of simulated points to the
//' # theoretical distribution of simulated points
//' # for large Bessel layers, it should look like a unconditional Brownian bridge
//' x <- 0.53
//' y <- 4.32
//' s <- 0.53
//' t <- 2.91
//' q <- 1.72
//' replicates <- 10000
//' paths <- list()
//' large_bessel_layer <- bessel_layer_simulation(x = x,
//'                                               y = y,
//'                                               s = s,
//'                                               t = t,
//'                                               mult = 100)
//' # repeatedly simulate Brownian bridge 
//' for (i in 1:replicates) {
//'   paths[[i]] <- layered_brownian_bridge(x = x,
//'                                         y = y,
//'                                         s = s,
//'                                         t = t,
//'                                         bessel_layer = large_bessel_layer,
//'                                         times = seq(s, t, 0.01))
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
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List layered_brownian_bridge(const double &x,
                                   const double &y,
                                   const double &s,
                                   const double &t,
                                   const Rcpp::List &bessel_layer,
                                   const Rcpp::NumericVector &times)
{
  if (t <= s) {
    stop("layeredBB::layered_brownian_bridge: t <= s. Must have s < t");
  } else if (Rcpp::min(times) < s) {
    stop("layeredBB::layered_brownian_bridge: minimum of specified times is less than s");
  } else if (Rcpp::max(times) > t) {
    stop("layeredBB::layered_brownian_bridge: maximum of specified times is greater than t");
  }
  while (true) {
    // add in line to abort C++ if user has pressed Ctrl/Cmd+C or Escape in R
    Rcpp::checkUserInterrupt();
    const Rcpp::NumericVector u = Rcpp::runif(2, 0.0, 1.0);
    double l1, l2, v1, v2;
    Rcpp::List simulated_BB;
    if (u[0] < 0.5) {
      Rcpp::NumericVector sim_m = min_sampler(x, y, s, t, bessel_layer["L"], bessel_layer["l"]);
      l1 = sim_m["min"], l2 = sim_m["min"];
      v1 = bessel_layer["u"], v2 = bessel_layer["U"];
      // simulate Bessel bridge conditional on the minimum point at intermediate points
      simulated_BB = min_Bessel_bridge_path_sampler(x, y, s, t, sim_m["min"], sim_m["tau"], times);
      const Rcpp::NumericMatrix &full_path = simulated_BB["full_path"];
      // check that none of the simulated points are outside of the layer: 
      // if so, resample Bessel bridge
      if (Rcpp::max(full_path.row(0)) > v2) {
        while (true) {
          simulated_BB = min_Bessel_bridge_path_sampler(x, y, s, t, sim_m["min"], sim_m["tau"], times);
          const Rcpp::NumericMatrix &full_path = simulated_BB["full_path"];
          if (Rcpp::max(full_path.row(0)) <= v2) {
            break;
          }
        }
      }
    } else {
      Rcpp::NumericVector sim_m = max_sampler(x, y, s, t, bessel_layer["u"], bessel_layer["U"]);
      l1 = bessel_layer["l"], l2 = bessel_layer["L"];
      v1 = sim_m["max"], v2 = sim_m["max"];
      // simulate Bessel bridge conditional on the maximum at intermediate points
      simulated_BB = max_Bessel_bridge_path_sampler(x, y, s, t, sim_m["max"], sim_m["tau"], times);
      const Rcpp::NumericMatrix &full_path = simulated_BB["full_path"];
      // check that none of the simulated points are outside of the layer: 
      // if so, resample Bessel bridge
      if (Rcpp::min(full_path.row(0)) < l2) {
        while (true) {
          simulated_BB = max_Bessel_bridge_path_sampler(x, y, s, t, sim_m["max"], sim_m["tau"], times);
          const Rcpp::NumericMatrix &full_path = simulated_BB["full_path"];
          if (Rcpp::min(full_path.row(0)) >= l2) {
            break;
          }
        }
      }
    }
    // checking if BB remains in [l1, v1]
    const Rcpp::NumericMatrix &full_path = simulated_BB["full_path"];
    double d1 = fabs(v1-l1);
    int j = ceil(sqrt(t-s + d1*d1) / (2*d1));
    // decide whether or not to accept simmulated Brownian bridge
    if (u[0] < 0.5) {
      // flip delta coin #1
      if (delta_coin_intervals(u[1], j, full_path.row(0), full_path.row(1), l1, v1)) {
        return simulated_BB;
      } else {
        // checking if BB remains in [l2, v2]
        double d2 = fabs(v2-l2);
        int k = ceil(sqrt(t-s + d2*d2) / (2*d2));
        // flip detla coin #2
        if (delta_coin_intervals(u[1], k, full_path.row(0), full_path.row(1), l2, v2)) {
          if (Rcpp::runif(1, 0.0, 1.0)[0] < 0.5) {
            return simulated_BB;
          }
        }
      }
    } else {
      // flip delta coin #1
      if (delta_coin_intervals(u[1], j, -full_path.row(0), full_path.row(1), -v1, -l1)) {
        return simulated_BB;
      } else {
        // checking if BB remains in [l2, v2]
        double d2 = fabs(v2-l2);
        int k = ceil(sqrt(t-s + d2*d2) / (2*d2));
        // flip detla coin #2
        if (delta_coin_intervals(u[1], k, -full_path.row(0), full_path.row(1), -v2, -l2)) {
          if (Rcpp::runif(1, 0.0, 1.0)[0] < 0.5) {
            return simulated_BB;
          }
        }
      }
    }
  }
}

//' Multi-dimensional Layered Brownian Bridge sampler
//'
//' Simulation of a multi-dimensional layered Brownian Bridge given Bessel layers
//' at user-specified times
//'
//' @param dim dimension of Brownian bridge
//' @param x start value of Brownian bridge
//' @param y end value of Brownian bridge
//' @param s start value of Brownian bridge
//' @param t end value of Brownian bridge
//' @param bessel_layers a list of length dim where list[i] is the Bessel layer for component i
//' @param times vector of real numbers to simulate Bessel bridge
//' 
//' @return A list with the following components
//' \describe{
//'   \item{full_path}{Matrix of the simulated layered Brownian bridge path at all 
//'                    included time points, i.e. s, t and times. The times
//'                    are sorted and duplicates are removed. The first dim rows
//'                    are the points of the Brownian bridge in each component,
//'                    last row gives the corresponding times}
//'   \item{simulated_path}{Matrix of the simulated layered Brownian bridge path 
//'                         only at the specified times passed into the function, 
//'                         i.e. the times vector. The times are not sorted and
//'                         duplicates are not removed. The first dim rows
//'                         are the points of the layered Brownian bridge in 
//'                         each component, last row gives the corresponding times}
//' }
//'
//' @examples
//' # simulate Bessel layer for two-dimensional Brownian bridge starting 
//' # and ending at (0,0) in time [0,1]
//' bes_layers <- multi_bessel_layer_simulation(dim = 2,
//'                                             x = c(0, 0),
//'                                             y = c(0, 0),
//'                                             s = 0,
//'                                             t = 1,
//'                                             mult = 0.5)
//' # simulate two-dimensional Brownian bridge starting 
//' # and ending at (0,0) in time [0,1]
//' multi_layered_brownian_bridge(dim = 2,
//'                               x = c(0,0),
//'                               y = c(0,0),
//'                               s = 0,
//'                               t = 1,
//'                               bessel_layers = bes_layers,
//'                               times = c(0.2, 0.4, 0.6, 0.8))
//'
//' # note that simulated_path does not remove duplicates passed into times
//' multi_layered_brownian_bridge(dim = 2,
//'                               x = c(0,0),
//'                               y = c(0,0),
//'                               s = 0,
//'                               t = 1,
//'                               bessel_layers = bes_layers,
//'                               times = c(0.2, 0.4, 0.6, 0.8, 0.4, 0.6))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List multi_layered_brownian_bridge(const int &dim,
                                         const arma::vec &x,
                                         const arma::vec &y,
                                         const double &s,
                                         const double &t,
                                         const Rcpp::List &bessel_layers,
                                         const Rcpp::NumericVector &times)
{
  // check that x and y match the dimensions of dim
  if (t <= s) {
    stop("layeredBB::multi_layered_brownian_bridge: t <= s. Must have s < t");
  } if (x.size() != dim) {
    stop("layeredBB::multi_layered_brownian_bridge: size of x is not equal to dim");
  } else if (y.size() != dim) {
    stop("layeredBB::multi_layered_brownian_bridge: size of y is not equal to dim");
  } else if (bessel_layers.size() != dim) {
    stop("layeredBB::multi_layered_brownian_bridge: size of bessel_layers is not equal to dim");
  } else if (Rcpp::min(times) < s) {
    stop("layeredBB::multi_layered_brownian_bridge: minimum of specified times is less than s");
  } else if (Rcpp::max(times) > t) {
    stop("layeredBB::multi_layered_brownian_bridge: maximum of specified times is greater than t");
  } 
  // ----- collect all times into one vector
  // and remove duplicates and sort full_times vector
  Rcpp::NumericVector full_times = times;
  full_times.push_front(s);
  full_times.push_back(t);
  full_times = Rcpp::sort_unique(full_times);
  full_times = Rcpp::unique(full_times);
  // for component, we simulate a layered Brownian bridge
  // ----- for each component, we simulate a Brownian bridge
  // BB at all times included
  Rcpp::NumericMatrix multi_full_bb(dim+1, full_times.size());
  multi_full_bb.row(dim) = full_times;
  // BB at specified times
  Rcpp::NumericMatrix multi_simulated_bb(dim+1, times.size());
  multi_simulated_bb.row(dim) = times;
  // loop through the components and simulate a layered Brownian bridge
  // we keep the simulated values at the times we want
  // we 'throw away' auxiliary information about the maximum and minimum simulated
  for (int d=0; d < dim; ++d) {
    // getting layer information for component d
    const Rcpp::List &sublist = bessel_layers[d];
    // simulate layered Brownian bridge for component d
    Rcpp::List component_BB = layered_brownian_bridge(x.at(d), y.at(d), s, t, sublist, times);
    // fill in full path BB with the remove_m_path (i.e. path without simulated minimum or maximum point)
    const Rcpp::NumericMatrix &remove_m_path = component_BB["remove_m_path"];
    multi_full_bb.row(d) = remove_m_path.row(0);
    // fill in simulated BB
    const Rcpp::NumericMatrix &simulated_path = component_BB["simulated_path"];
    multi_simulated_bb.row(d) = simulated_path.row(0);
  }
  return Rcpp::List::create(Named("full_path", multi_full_bb),
                            Named("simulated_path", multi_simulated_bb));
}
