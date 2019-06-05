#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

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
  double min = x - (0.5*(sqrt((y-x)*(y-x)-2.0*(t-s)*log(u1)) - y + x));
  
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
    V = (1.0 / inv_gauss_sampler(mu, lambda));
  }
  
  // set tau (time of simualted minimum)
  double tau = ((s*V)+t)/(1.0+V);
  
  // setting simulated minimum and tau in array
  Rcpp::NumericVector simulated_min = Rcpp::NumericVector::create(Named("min", min), Named("tau", tau));
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
  
  // set simulated value
  double term1 = ((Wr-min)*fabs(tau-q)/(pow(fabs(tau-r), 1.5))) + b.at(0);
  double W = min + sqrt(fabs(tau-r)*(term1*term1 + b.at(1)*b.at(1) + b.at(2)*b.at(2)));
  
  return W;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix min_Bessel_bridge_path_sampler(const double &x, const double &y,
                                                   const double &s, const double &t,
                                                   const double &min, const double &tau,
                                                   Rcpp::NumericVector times)
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
  Rcpp::NumericVector simulated_bb(times.size());
  
  // when simulating the path, want to work from left to right when left of the min
  // and work from right to left when right of the min
  // this is so we can use the Markov property
  
  for (int i = 0; times[i] <= tau; ++i) {
    // simulate the point at each time
    if (times[i] == s) {
      simulated_bb[i] = x;
    } else if (times[i] == tau) {
      simulated_bb[i] = min;
    } else {
      // if left of tau, then we simulate the next point normally, starting from the previous point to the end
      simulated_bb[i] = min_Bessel_bridge_sampler(simulated_bb[i-1], y, times[i-1], t, min, tau, times[i]);
    }
  }
  
  for (int i = times.size()-1; times[i] > tau; --i) {
    // simulate the point at each time
    if (times[i] == t) {
      simulated_bb[i] = y;
    } else {
      // reflect it by taking the absolute value of the (time wanted - t)
      simulated_bb[i] = min_Bessel_bridge_sampler(simulated_bb[i+1], x, fabs(times[i+1]-t), fabs(s-t), 
                                                  min, fabs(tau-t), fabs(times[i]-t));
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
  return Rcpp::NumericVector::create(Named("max", -sim_min["min"]), Named("tau", sim_min["tau"]));
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
                                                   Rcpp::NumericVector times)
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

// [[Rcpp::export]]
double sigma_bar(const double &j, const double &x, const double &y,
                 const double &s, const double &t,
                 const double &l, const double &v)
{
  return exp((-2.0/(t-s))*((fabs(v-l)*j)+std::min(l,v)-x)*((fabs(v-l)*j)+std::min(l,v)-y));
}

// [[Rcpp::export]]
double sigma(const double &j, const double &x, const double &y, 
             const double &s, const double &t,
             const double &l, const double &v)
{
  return sigma_bar(j,x,y,s,t,l,v) + sigma_bar(j,-x,-y,s,t,-l,-v);
}

// [[Rcpp::export]]
double phi_bar(const double &j, const double &x, const double &y,
               const double &s, const double &t,
               const double &l, const double &v)
{
  return exp(-(2.0*j/(t-s)) * (fabs(v-l)*fabs(v-l)*j + fabs(v-l)*(x-y))); 
}

// [[Rcpp::export]]
double phi(const double &j, const double &x, const double &y, 
           const double &s, const double &t,
           const double &l, const double &v)
{
  return phi_bar(j,x,y,s,t,l,v) + phi_bar(j,-x,-y,s,t,-l,-v);
}

// [[Rcpp::export]]
double psi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v)
{
  return (((2*fabs(v-min)*j)-(std::max(x,y)-min))*(exp(-(2.0*fabs(v-min)*j/(t-s))*((fabs(v-min)*j)-(std::max(x,y)-min)))));
}

// [[Rcpp::export]]
double chi(const double &j, const double &x, const double &y,
           const double &s, const double &t,
           const double &min, const double &v)
{
  return (((2*fabs(v-min)*j)+(std::max(x,y)-min))*(exp(-(2.0*fabs(v-min)*j/(t-s))*((fabs(v-min)*j)+(std::max(x,y)-min)))));
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SgammaK_intervals(const int &k, const double &x, const double &y,
                                           const double &s, const double &t,
                                           const double &l, const double &v)
{
  // function calculates (S_{2k+1}^{gamma}, S_{2k}^{gamma}) for a given k
  
  // if k==0, then we need to calculate (S_{1}^{gamma}, S_{0}^{gamma}) = (1-sigma(1), 1)
  if (k == 0) {
    return Rcpp::NumericVector::create(1-sigma(1.0,x,y,s,t,l,v), 1);
  } 
  
  // calculating S_{2k}^{gamma} = 1 - sum_{j=1}^{k} sigma(j)-phi(j)
  double S_2k = 1;
  for (int j=1; j <= k; ++j) {
    S_2k -= (sigma(j,x,y,s,t,l,v) - phi(j,x,y,s,t,l,v));
  }
  // calculating S_{2k+1}^{gamma} = S_{2k} - sigma(k+1)
  double S_2k_plus_1 = S_2k - sigma(k+1,x,y,s,t,l,v);
  
  return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_1_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v)
{
  // function calculates (S_{2k+1}^{delta,1}, S_{2k}^{delta,1}) for a given k
  // S_{k}^{delta,1} = S_{k}^{gamma}/denom
  double denom = 1 - exp(-2.0*(x-min)*(y-min)/(t-s));
  // use calc_SgammaK_intervals and divide the result by denom
  Rcpp::NumericVector SgammaK = calc_SgammaK_intervals(k,x,y,s,t,min,v);
  
  return Rcpp::NumericVector::create(SgammaK.at(0)/denom, SgammaK.at(1)/denom);
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_2_intervals(const int &k, const double &x, const double &y,
                                             const double &s, const double &t,
                                             const double &min, const double &v)
{
  // function calculates (S_{2k+1}^{delta,2}, S_{2k}^{delta,2}) for a given k
  // checking k is large enough for this to be valid
  double K = sqrt((t-s)+(fabs(v-min)*fabs(v-min)))/(2*fabs(v-min));
  if (k < K) {
    stop("error in calc_SdeltaK_2_intervals: given k is too small");
  }
  
  // calculating S_{2k+1}^{delta,2} = 1 - (psi(j)-chi(j)) / abs(x-y)
  double S_2k = 1;
  for (int j=1; j <= k; ++j) {
    S_2k -= ((psi(j,x,y,s,t,min,v)-chi(j,x,y,s,t,min,v))/fabs(x-y));
  }
  double S_2k_plus_1 = S_2k - (psi(k+1,x,y,s,t,min,v)/fabs(x-y));
  
  return Rcpp::NumericVector::create(S_2k_plus_1, S_2k);
}

// [[Rcpp::export]]
Rcpp::NumericVector calc_SdeltaK_intervals(const int &k, const double &x, const double &y,
                                           const double &s, const double &t,
                                           const double &min, const double &v)
{
  if (std::min(x,y) > min) {
    return calc_SdeltaK_1_intervals(k,x,y,s,t,min,v);
  } else if (std::min(x,y) == min) {
    return calc_SdeltaK_2_intervals(k,x,y,s,t,min,v);
  } else {
    stop("error in calc_SdeltaK_intervals: min(x,y) < min - given minimum point is not the minimum of the BB");
  }
}

// [[Rcpp::export]]
double product_vector_elements(const Rcpp::NumericVector &vect) {
  double prod = 1;
  for (const auto &element: vect) {
    prod *= element;
  }
  return prod;
}

// [[Rcpp::export]]
bool gamma_coin(const double &x, const double &y, 
                const double &s, const double &t,
                const double &l, const double &v,
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

// [[Rcpp::export]]
bool gamma_coin_intervals(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, const Rcpp::NumericVector &t,
                          const double &l, const double &v, int k)
{
  // check if vector lengths are all the same
  if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
    stop("error in gamma_coin_intervals: vector lengths are not equal");
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

// [[Rcpp::export]]
bool delta_coin(const double &x, const double &y, 
                const double &s, const double &t,
                const double &min, const double &v,
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

// [[Rcpp::export]]
bool delta_coin_intervals(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y,
                          const Rcpp::NumericVector &s, const Rcpp::NumericVector &t,
                          const double &min, const double &v, int k)
{
  // check if vector lengths are all the same
  if (x.size()!=y.size() || x.size()!=s.size() || x.size()!=t.size()) {
    stop("error in gamma_coin_intervals: vector lengths are not equal");
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

// [[Rcpp::export]]
List bessel_layer_simulation(const double &x, const double &y,
                             const double &s, const double &t,
                             Rcpp::NumericVector &a)
{
  // initialise Bessel layer
  int l = 0;
  
  while (true) {
    // flip gamma coin to determine if BB stays within (min(x,y)-a[l], max(x,y)+a[l])
    if (!(x==y && a.at(l)==0)) {
      if (gamma_coin(x,y,s,t,std::min(x,y)-a.at(l),std::max(x,y)+a.at(l),0)) {
        // if true, then return current Bessel layer (l) and current vector (a) as a list
        // we return (l+1), since indicies start from 0 in C++ but start from 1 in R
        return List::create(_["a"] = a, _["l"] = l+1);
      }
    }
    
    // if false, l = l+1
    l += 1;
    // need to check that we are not at the end of the sequence (a)
    // in this case, we extend the sequence of (a), and try carry on to find a layer
    if (l >= a.size()) {
      // find the size of the vector (a) and set it to (size)
      auto size = a.size();
      // find the last element of the vector (a) and set it to (last)
      // since indices start from 0 in C++, this is the (size-1)-th element
      double last = a[(size-1)];
      for (int i = 0; i < size; ++i) {
        a.push_back(a[i] + last);
      }
      // remove any duplicates
      a.erase(std::unique(a.begin(), a.end()), a.end());
    }
  }
}

// [[Rcpp::export]]
double find_max(const Rcpp::NumericVector vect) 
{
  // finds the maximum value in a vector
  double current_max = vect.at(0); // first element
  for (const auto &element: vect) {
    if (element > current_max) {
      current_max = element;
    }
  }
  return current_max;
}

// [[Rcpp::export]]
double find_min(const Rcpp::NumericVector vect)
{
  // finds the minimum value in a vector
  double current_min = vect.at(0); // first element
  for (const auto &element: vect) {
    if (element < current_min) {
      current_min = element;
    }
  }
  return current_min;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix layered_brownian_bridge(const double &x, const double &y,
                                            const double &s, const double &t,
                                            const Rcpp::NumericVector &a, int l,
                                            const Rcpp::NumericVector &times)
{
  // layered Brownian bridge sampler
  // simulate from uniform random variable
  double u = Rcpp::runif(1, 0.0, 1.0)[0];
  // we minus 1 from the Bessel layer (l), since indicies start from 0 in C++ but start from 1 in R
  l -= 1;
  
  while (true) {
    // add in line to abort C++ if user has pressed Ctrl/Cmd+C or Escape in R
    Rcpp::checkUserInterrupt();
    
    // declare variables required
    double l1, l2, v1, v2;
    Rcpp::NumericMatrix simulated_BB;
    // Rcpp::NumericVector bridge_times = times;
    
    if (u < 0.5) {
      // simulate minimum point
      Rcpp::NumericVector simulated_min;
      
      // setting l1, l2, v1, v2 values in the different cases
      // simulating minimum points if needed
      if (l==0 && a.at(0)==0) {
        // simulated Bessel layer was [min(x,y), max(x,y)]
        // set minimum to occur at min(x,y)
        if (x < y) {
          // mimimum occurs at (x) at time (s)
          l1 = x, l2 = x;
          v1 = y, v2 = y;
          simulated_min = Rcpp::NumericVector::create(Named("min", x), Named("tau", s));
        } else {
          // minimum occurs at (y) at time (t)
          l1 = y, l2 = y;
          v1 = x, v2 = x;
          simulated_min = Rcpp::NumericVector::create(Named("min", y), Named("tau", t));
        }
      } else if (l==0 && a.at(0)!=0) {
        // simulated Bessel layer was [min(x,y)-a[0], max(x,y)+a[0]]
        // but this is the first layer, don't have a[l-1]
        simulated_min = min_sampler(x, y, s, t, std::min(x,y)-a.at(l), std::min(x,y));
        l1 = simulated_min["min"], l2 = simulated_min["min"];
        v1 = std::max(x,y)+a.at(l), v2 = std::max(x,y)+a.at(l);
      } else {
        // simulated Bessel layer was [min(x,y)-a[l], max(x,y)+a[l]]
        simulated_min = min_sampler(x, y, s, t, std::min(x,y)-a.at(l), std::min(x,y)-a.at(l-1));
        l1 = simulated_min["min"], l2 = simulated_min["min"];
        v1 = std::max(x,y)+a.at(l-1), v2 = std::max(x,y)+a.at(l);
      }
      
      // simulate Bessel bridge conditional on the minimum point at intermediate points
      simulated_BB = min_Bessel_bridge_path_sampler(x, y, s, t, simulated_min["min"], simulated_min["tau"], times);
      
      // check that none of the simulated points are outside of the layer: if so, resample Bessel bridge
      while (find_max(simulated_BB(0, _)) > v2) { 
        simulated_BB = min_Bessel_bridge_path_sampler(x, y, s, t, simulated_min["min"], simulated_min["tau"], times);
      }
    } else {
      // simulate maximum point
      Rcpp::NumericVector simulated_max;
      if (l==0 && a.at(0)==0) {
        // simulated Bessel layer was [min(x,y), max(x,y)]
        // set maximum to occur at max(x,y)
        if (x < y) {
          // maximum occurs at (y) at time (t)
          l1 = x, l2 = x;
          v1 = y, v2 = y;
          simulated_max = Rcpp::NumericVector::create(Named("max", y), Named("tau", t));
        } else {
          // maximum occurs at (x) at time (s)
          l1 = y, l2 = y;
          v1 = x, v2 = x;
          simulated_max = Rcpp::NumericVector::create(Named("max", x), Named("tau", s));
        }
      } else if (l==0 && a.at(0)!=0) {
        // simulated Bessel layer was [min(x,y)-a[0], max(x,y)+a[0]]
        // but this is the first layer, don't have a[l-1]
        simulated_max = max_sampler(x, y, s, t, std::max(x,y), std::max(x,y)+a.at(l));
        l1 = std::min(x,y)-a.at(l), l2 = std::min(x,y)-a.at(l);
        v1 = simulated_max["max"], v2 = simulated_max["max"];
      } else {
        // simulated Bessel layer was [min(x,y)-a[l], max(x,y)+a[l]]
        simulated_max = max_sampler(x, y, s, t, std::max(x,y)+a.at(l-1), std::max(x,y)+a.at(l));
        l1 = std::min(x,y)-a.at(l-1), l2 = std::min(x,y)-a.at(l);
        v1 = simulated_max["max"], v2 = simulated_max["max"];
      }
      
      // simulate Bessel bridge conditional on the maximum at intermediate points
      simulated_BB = max_Bessel_bridge_path_sampler(x, y, s, t, simulated_max["max"], simulated_max["tau"], times);
      
      // check that none of the simulated points are outside of the layer: if so, resample Bessel bridge
      while (find_min(simulated_BB(0, _)) < l2) {
        simulated_BB = max_Bessel_bridge_path_sampler(x, y, s, t, simulated_max["max"], simulated_max["tau"], times);
      }
    }
    
    // checking if BB remains in [l1, v1]
    double starting_index_1 = sqrt((t-s)+fabs(v1-l1)*fabs(v1-l1))/(2*fabs(v1-l1));
    int j = ceil(starting_index_1);
    
    Rcpp::NumericVector x_values, y_values, s_values, t_values;
    for (int i=0; i < (simulated_BB.ncol()-1); ++i) {
      x_values.push_back(simulated_BB(0, i));
      y_values.push_back(simulated_BB(0, i+1));
      s_values.push_back(simulated_BB(1, i));
      t_values.push_back(simulated_BB(1, i+1));
    }
    
    if (u < 0.5) {
      // flip delta coin #1
      if (delta_coin_intervals(x_values, y_values, s_values, t_values, l1, v1, j)) {
        return simulated_BB;
      } else {
        // checking if BB remains in [l1, v1]
        double starting_index_2 = sqrt((t-s)+fabs(v2-l2)*fabs(v2-l2))/(2*fabs(v2-l2));
        int k = ceil(starting_index_2);
        
        // flip detla coin #2
        if (delta_coin_intervals(x_values, y_values, s_values, t_values, l2, v2, k)) {
          return simulated_BB;
        } else {
          continue;
        }
      } 
    } else {
      // flip delta coin #1
      if (delta_coin_intervals(-x_values, -y_values, s_values, t_values, -v1, -l1, j)) {
        return simulated_BB;
      } else {
        // checking if BB remains in [l1, v1]
        double starting_index_2 = sqrt((t-s)+fabs(v2-l2)*fabs(v2-l2))/(2*fabs(v2-l2));
        int k = ceil(starting_index_2);
        
        // flip detla coin #2
        if (delta_coin_intervals(-x_values, -y_values, s_values, t_values, -v2, -l2, k)) {
          return simulated_BB;
        } else {
          continue;
        }
      }
    }
  }	
}
