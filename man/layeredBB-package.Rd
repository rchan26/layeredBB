\name{layeredBB-package}
\alias{layeredBB-package}
\alias{layeredBB}
\docType{package}
\title{
  Simulates layered Brownian bridges using C++ and Rcpp
}
\description{
  R access to C++ implementation of layered Brownian Bridge simulation and Bessel layer simulation.
}
\details{
  This package includes functions to simulated layered Brownian bridges using Bessel layers. See Pollock et al. (2016) and Pollock (2013) for more details on simulating Bessel layers using Cauchy sums. 
}
\author{
  Ryan Chan, rchan@turing.ac.uk.
}
\references{
  M. Pollock. Some Monte Carlo Methods for Jump Diffusions. PhD thesis, University of Warwick, September 2013. 
  
  M.Pollock, A.M. Johansen, G.O. Roberts. 2016. On the Exact and epsilon-Strong Simulation of (Jump) Diffusions. Bernoulli. 22(2), pp. 794-856. 
  
  Beskos, A., Papaspiliopoulos, O., Roberts, G.O., and Fearnhead, P. 2006. Exact and computationally effcient likelihood-based estimation for discretely observed diffusion processes (with discussion). Journal of Royal Statistical Society, B. 68, pp. 333-382.
}
\keyword{ package }
\examples{
# simulate Bessel layer
bes_layer <- bessel_layer_simulation(x = 0, y = 0, s = 0, t = 1, a = seq(0.1, 1.0, 0.1))
# simulate layered Brownian bridge
layered_brownian_bridge(x = 0, y = 0, s = 0, t = 1, a = bes_layer$a, l = bes_layer$l, times = seq(0.2, 0.8, 0.2))
  
}
