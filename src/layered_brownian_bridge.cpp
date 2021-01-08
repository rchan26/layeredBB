#include "../inc/layered_brownian_bridge.hpp"
#include "../inc/brownian_bridge_min.hpp"
#include "../inc/inverse_gauss.hpp"
#include "../inc/cauchy_sums.hpp"
#include "../inc/coins_cauchy.hpp"
#include "../inc/bessel_layer_sim.hpp"

using namespace Rcpp;

//' Find maximum of vector
//'
//' This function finds the maximum value of a numerical vector
//'
//' @param vect numerical vector
//'
//' @return maximum of vector given
//'
//' @examples
//' # returns 0.8
//' find_max(c(-0.4, 0.23, -0.4321, 0.6, 0.3, 0.8, 0.54)) 
//'
//' @export
// [[Rcpp::export]]
double find_max(const Rcpp::NumericVector vect) {
  // finds the maximum value in a vector
  double current_max = vect.at(0); // first element
  for (const auto &element: vect) {
    if (element > current_max) {
      current_max = element;
    }
  }
  return current_max;
}

//' Find minimum of vector
//'
//' This function finds the minimum value of a numerical vector
//'
//' @param vect numerical vector
//'
//' @return minimum of vector given
//' 
//' @examples
//' # returns -0.4321
//' find_min(c(-0.4, 0.23, -0.4321, 0.6, 0.3, 0.8, 0.54)) 
//'
//' @export
// [[Rcpp::export]]
double find_min(const Rcpp::NumericVector vect) {
  // finds the minimum value in a vector
  double current_min = vect.at(0); // first element
  for (const auto &element: vect) {
    if (element < current_min) {
      current_min = element;
    }
  }
  return current_min;
}

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
//' @param remove_m boolean to indicate whether or not simulated minimum or 
//'        maximum is removed (default is TRUE)
//' 
//' @return matrix of the simulated layered Brownian bridge path, 
//'         first row is points X, second row are corresponding times
//'
//' @examples
//' # simulate Bessel layer
//' bes_layer <- bessel_layer_simulation(x = 0,
//'                                      y = 0,
//'                                      s = 0,
//'                                      t = 1,
//'                                      mult = 0.2)
//' # simulate layered Brownian bridge
//' # simulated minimum or maximum is removed
//' layered_brownian_bridge(x = 0,
//'                         y = 0,
//'                         s = 0,
//'                         t = 1,
//'                         bessel_layer = bes_layer,
//'                         times = seq(0.2, 0.8, 0.2))
//' # simulated minimum or maximum is kept
//' layered_brownian_bridge(x = 0,
//'                         y = 0,
//'                         s = 0,
//'                         t = 1,
//'                         bessel_layer = bes_layer,
//'                         times = seq(0.2, 0.8, 0.2),
//'                         remove_m = FALSE)
//' 
//' # another example
//' # simulated minimum or maximum is removed
//' start <- runif(1, -1, 1)
//' end <- runif(1, -1, 1)
//' bes_layer <- bessel_layer_simulation(x = start, y = end, s = 0, t = 1, mult = 0.2)
//' path <- layered_brownian_bridge(x = start,
//'                                 y = end,
//'                                 s = 0,
//'                                 t = 1,
//'                                 bessel_layer = bes_layer,
//'                                 times = seq(0, 1, 0.01))
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X',
//'      ylim = c(bes_layer$L, bes_layer$U))
//' lines(x = path['time',], y = path['X',])
//' abline(h=c(bes_layer$L, bes_layer$U), col = 'red')
//' abline(h=c(bes_layer$l, bes_layer$u), col = 'red', lty = 2)
//' 
//' # simulated miniumum or maximum is kept
//' start <- runif(1, -1, 1)
//' end <- runif(1, -1, 1)
//' bes_layer <- bessel_layer_simulation(x = start, y = end, s = 0, t = 1, mult = 0.2)
//' path <- layered_brownian_bridge(x = start,
//'                                 y = end,
//'                                 s = 0,
//'                                 t = 1,
//'                                 bessel_layer = bes_layer,
//'                                 times = seq(0, 1, 0.01),
//'                                 remove_m = FALSE)
//' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X',
//'      ylim = c(bes_layer$L, bes_layer$U))
//' lines(x = path['time',], y = path['X',])
//' abline(h=c(bes_layer$L, bes_layer$U), col = 'red')
//' abline(h=c(bes_layer$l, bes_layer$u), col = 'red', lty = 2)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix layered_brownian_bridge(const double &x, 
                                            const double &y,
                                            const double &s, 
                                            const double &t,
                                            const Rcpp::List &bessel_layer,
                                            const Rcpp::NumericVector &times,
                                            const bool &remove_m = true)
{
  while (true) {
    // add in line to abort C++ if user has pressed Ctrl/Cmd+C or Escape in R
    Rcpp::checkUserInterrupt();
    double u = Rcpp::runif(1, 0.0, 1.0)[0];
    double l1, l2, v1, v2;
    Rcpp::NumericVector sim_m;
    Rcpp::NumericMatrix simulated_BB;
    if (u < 0.5) {
      sim_m = min_sampler(x, y, s, t, bessel_layer["L"], bessel_layer["l"]);
      l1 = sim_m["min"], l2 = sim_m["min"];
      v1 = bessel_layer["u"], v2 = bessel_layer["U"];
      // simulate Bessel bridge conditional on the minimum point at intermediate points
      simulated_BB = min_Bessel_bridge_path_sampler(x, y, s, t, sim_m["min"], sim_m["tau"], times);
      // check that none of the simulated points are outside of the layer: 
      // if so, resample Bessel bridge
      while (find_max(simulated_BB(0, _)) > v2) { 
        simulated_BB = min_Bessel_bridge_path_sampler(x, y, s, t, sim_m["min"], sim_m["tau"], times);
      }
    } else {
      sim_m = max_sampler(x, y, s, t, bessel_layer["u"], bessel_layer["U"]);
      l1 = bessel_layer["l"], l2 = bessel_layer["L"];
      v1 = sim_m["max"], v2 = sim_m["max"];
      // simulate Bessel bridge conditional on the maximum at intermediate points
      simulated_BB = max_Bessel_bridge_path_sampler(x, y, s, t, sim_m["max"], sim_m["tau"], times);
      // check that none of the simulated points are outside of the layer: 
      // if so, resample Bessel bridge
      while (find_min(simulated_BB(0, _)) < l2) {
        simulated_BB = max_Bessel_bridge_path_sampler(x, y, s, t, sim_m["max"], sim_m["tau"], times);
      }
    }
    // checking if BB remains in [l1, v1]
    double d1 = fabs(v1-l1);
    int j = ceil(sqrt(t-s + d1*d1) / (2*d1));
    // decide whether or not to accept simmulated Brownian bridge
    bool success = false;
    if (u < 0.5) {
      // flip delta coin #1
      if (delta_coin_intervals(j, simulated_BB(0,_), simulated_BB(1,_), l1, v1)) {
        success = true;
      } else {
        // checking if BB remains in [l2, v2]
        double d2 = fabs(v2-l2);
        int k = ceil(sqrt(t-s + d2*d2) / (2*d2));
        // flip detla coin #2
        if (delta_coin_intervals(k, simulated_BB(0,_), simulated_BB(1,_), l2, v2)) {
          if (Rcpp::runif(1, 0.0, 1.0)[0] < 0.5) {
            success = true;
          }
        } else {
          continue;
        }
      } 
    } else {
      // flip delta coin #1
      if (delta_coin_intervals(j, -simulated_BB(0,_), simulated_BB(1,_), -v1, -l1)) {
        success = true;
      } else {
        // checking if BB remains in [l2, v2]
        double d2 = fabs(v2-l2);
        int k = ceil(sqrt(t-s + d2*d2) / (2*d2));
        // flip detla coin #2
        if (delta_coin_intervals(k, -simulated_BB(0,_), simulated_BB(1,_), -v2, -l2)) {
          if (Rcpp::runif(1, 0.0, 1.0)[0] < 0.5) {
            success = true;
          } 
        } else {
          continue;
        }
      }
    }
    // if accepted, remove simulated minimum or maximum if remove_m == TRUE
    if (success) {
      if (remove_m) {
        Rcpp::NumericMatrix simulated_BB_removed_m(simulated_BB.nrow(), simulated_BB.ncol()-1);
        int col_index = 0;
        // want to only keep the user-specified times and the start and end time points
        Rcpp::NumericVector keep_times = times;
        keep_times.push_back(s);
        keep_times.push_back(t);
        for (int j=0; j < simulated_BB.ncol(); ++j) {
          if (std::find(keep_times.begin(), keep_times.end(), simulated_BB(1,j)) != keep_times.end()) {
            // check that the time is in the keep_times vector
            for (int i=0; i < simulated_BB_removed_m.nrow(); ++i) {
              simulated_BB_removed_m(i, col_index) = simulated_BB(i,j);
            }
            col_index += 1;
          }
        }
        rownames(simulated_BB_removed_m) = CharacterVector::create("X", "time");
        return simulated_BB_removed_m;
      } else {
        return simulated_BB;
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
//' @return matrix of the simulated layered Brownian bridge path, 
//'         first dim rows are points for X in each component, 
//'         last row are corresponding times
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
//'                               times = seq(0.2, 0.8, 0.2))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix multi_layered_brownian_bridge(const int &dim,
                                                  const arma::vec &x,
                                                  const arma::vec &y,
                                                  const double &s,
                                                  const double &t,
                                                  const Rcpp::List &bessel_layers,
                                                  Rcpp::NumericVector times)
{
  // check that x and y match the dimensions of dim
  if (x.size() != dim) {
    stop("multi_layered_brownian_bridge: size of x is not equal to dim");
  } else if (y.size() != dim) {
    stop("multi_layered_brownian_bridge: size of y is not equal to dim");
  } else if (bessel_layers.size() != dim) {
    stop("multi_layered_brownian_bridge: size of bessel_layers is not equal to dim");
  }
  // collect all times into one vector
  times.insert(times.end(), s);
  times.insert(times.end(), t);
  // sort the vector 'times' forward in time
  times.sort();
  // delete any duplicates
  times.erase(std::unique(times.begin(), times.end()), times.end());
  // for component, we simulate a layered Brownian bridge
  // multi_BB is a matrix with dimensions (dim+1) x times.size()
  Rcpp::NumericMatrix multi_BB(dim+1, times.size());
  multi_BB.row(dim) = times;
  // loop through the components and simulate a layered Brownian bridge
  // we keep the simulated values at the times we want
  // we 'throw away' auxiliary information about the maximum and minimum simualted
  for (int d=0; d < dim; ++d) {
    // getting layer information for component d
    Rcpp::List sublist = bessel_layers[d];
    // simulate layered Brownian bridge for component d
    Rcpp::NumericMatrix component_BB = layered_brownian_bridge(x.at(d), y.at(d), s, t, sublist, times, true);
    multi_BB(d, _) = component_BB(0, _);
  }
  return(multi_BB);
}

// // function to pick out the associated values of the layered Brownian bridge at the Poisson points
// // [[Rcpp::export]]
// arma::mat get_pois_points_multi(const Rcpp::NumericMatrix &skeleton,
//                                 const Rcpp::NumericVector &pois_times) {
//   arma::mat pois_points(skeleton.nrow()-1, pois_times.size());
//   int col_index = 0;
//   // loop through the columns of the skeleton
//   for (int j=0; j < skeleton.ncol(); ++j) {
//     if (std::find(pois_times.begin(), pois_times.end(), skeleton((skeleton.nrow()-1),j)) != pois_times.end()) {
//       // if the current time is in pois_times, then we store the corresponding values of the Brownian bridge into another matrix
//       // loop through the rows and storing the values into the matrix pois_points
//       for (int i=0; i < skeleton.nrow()-1; ++i) {
//         pois_points(i, col_index) = skeleton(i, j);
//       }
//       col_index = col_index + 1;
//     }
//   }
//   return(pois_points);
// }
// 
// // function to pick out the associated values of the layered Brownian bridge at Poisson points by
// // accessing elements quickly if the times for the skeleton are sorted
// // [[Rcpp::export]]
// Rcpp::NumericMatrix get_pois_points_fast(Rcpp::NumericMatrix &skeleton) {
//   Rcpp::NumericMatrix pois_points = skeleton(Rcpp::Range(0, skeleton.nrow()-2),
//                                              Rcpp::Range(1, skeleton.ncol()-2));
//   return(pois_points);
// }
