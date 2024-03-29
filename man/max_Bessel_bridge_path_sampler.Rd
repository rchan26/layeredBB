% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{max_Bessel_bridge_path_sampler}
\alias{max_Bessel_bridge_path_sampler}
\title{Bessel Bridge path sampler given maximum}
\usage{
max_Bessel_bridge_path_sampler(x, y, s, t, m, tau, times, checks = TRUE)
}
\arguments{
\item{x}{start value of Bessel bridge}

\item{y}{end value of Bessel bridge}

\item{s}{start time of Bessel bridge}

\item{t}{end time of Bessel bridge}

\item{m}{maximum point}

\item{tau}{time of maximum point}

\item{times}{vector of real numbers to simulate Bessel bridge}

\item{checks}{logical value to determine if arguments that are passed into 
the function are checked. Things that are checked include that 
s < t, that requested simulation times are in [s,t], that m >= max(x,y)
and that if tau == s or tau == t, then m == x or m == y, respectivelys}
}
\value{
A list with the following components
\describe{
  \item{full_path}{Matrix of the simulated Bessel bridge path at all 
                   included time points, i.e. s, t and times. The times
                   are sorted and duplicates are removed. The first row
                   are the points of the Brownian bridge (named 'X') 
                   second row are corresponding times (named 'time')}
  \item{simulated_path}{Matrix of the simulated Bessel bridge path only at 
                        the specified times passed into the function, i.e. 
                        the times vector. The times are not sorted and
                        duplicates are not removed. The first row
                        are the points of the Bessel bridge (named 'X') 
                        second row are corresponding times (named 'time')}
  \item{remove_m_path}{Matrix of the simulated Bessel bridge path only at 
                       all included times points excluding tau. These times
                       are sorted and duplicates are removed. The first row
                       are the points of the Bessel bridge (named 'X') 
                       second row are corresponding times (named 'time'). 
                       Note that the maximum point is included if it is 
                       passed into the times vector}
}
}
\description{
Simulation of a path of a Bessel bridge at given times, given maximum occurs at time tuu
}
\examples{
# simulating a path at times = c(0.2, 0.4, 0.8) for a Bessel bridge starting 
# at 0 and ending at 0 in time [0,1] given maximum is at 0.4 at time 0.6
max_Bessel_bridge_path_sampler(x = 0,
                               y = 0,
                               s = 0,
                               t = 1,
                               m = 0.4,
                               tau = 0.6,
                               times = c(0.2, 0.4, 0.8))

# note that remove_m_path will still include the minimum if passed into times
# also note that simulated_path does not remove duplicates passed into times
max_Bessel_bridge_path_sampler(x = 0,
                               y = 0,
                               s = 0,
                               t = 1,
                               m = 0.4,
                               tau = 0.6,
                               times = c(0.2, 0.4, 0.6, 0.8, 0.6))

# another example
start <- runif(1, -1, 1)
end <- runif(1, -1, 1)
max <- max_sampler(x = start,
                   y = end,
                   s = 0,
                   t = 1,
                   low_bound = max(start, end)+0.2,
                   up_bound = max(start, end)+0.4)
path <- max_Bessel_bridge_path_sampler(x = start,
                                       y = end,
                                       s = 0,
                                       t = 1,
                                       m = max['max'],
                                       tau = max['tau'],
                                       times = seq(0, 1, 0.01))$full_path
plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X')
lines(x = path['time',], y = path['X',])
points(x = max['tau'], y = max['max'], col = 'red', pch = 20)
}
