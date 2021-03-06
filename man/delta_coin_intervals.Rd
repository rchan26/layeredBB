% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{delta_coin_intervals}
\alias{delta_coin_intervals}
\title{Delta coin flipper for intervals (used for Algorithm 33 in ST329)}
\usage{
delta_coin_intervals(u, k, X, times, min, v)
}
\arguments{
\item{u}{simulated value from random U[0,1]}

\item{k}{integer value starting index for calculating the intervals}

\item{X}{vector of values of Brownian bridge}

\item{times}{vector of times}

\item{min}{minimum of Brownian bridge}

\item{v}{upper bound of Brownian bridge}
}
\value{
boolean value: if T, accept probability that Brownian bridge with 
        minimum, min, remains in [min,v], otherwise reject
}
\description{
Flips 'Delta coin' for intervals; takes the product of the Cauchy sequence S^{delta}_{k} to 
determine whether or not the Brownian bridge with minimum, min, remains in the interval [l,v]
}
\examples{
# setting up Brownian bridge variable
brownian_bridge <- matrix(c(0, 0, -0.2, 0.4, 0.3, 0.5, 1, 1),
                          ncol = 4, nrow = 2)

# flip delta coin whether or not Brownian bridge remains in [-0.2, 1.5]
d <- abs(1.5 - -0.2)
k <- ceiling(sqrt(1 + d^2)/(2*d))
delta_coin_intervals(u = runif(1, 0, 1),
                     k = k,
                     X = brownian_bridge[1,],
                     times = brownian_bridge[2,],
                     min = -0.2,
                     v = 1.5)

}
