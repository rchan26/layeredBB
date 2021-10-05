#' Brownian Motion path sampler 
#'
#' Simulation of a path of a Brownian motion at given times
#'
#' @param x start value of Brownian motion
#' @param times vector of real numbers to simulate Brownian motion
#'
#' @return Matrix of the simulated Brownian motion path at all
#'         included time points. The times are sorted. 
#'         The first row are the points of the Brownian motion (named 'X')
#'         second row are corresponding times (named 'times')
#'
#' @examples
#' # simulating path for Brownian motion starting at (0,0) between [0,1]
#' path <- Brownian_motion_path_sampler(x = 0, times = seq(0, 1, 0.01))
#' plot(x = path['time',], y = path['X',], pch = 20, xlab = 'Time', ylab = 'X')
#' lines(x = path['time',], y = path['X',])
#' 
#' # comparing the simulated distribution of simulated points to the
#' # theoretical distribution of simulated points
#' # set variables
#' x <- 0
#' start_time <- 1.8
#' end_time <- 5
#' replicates <- 10000
#' paths <- list()
#' # repeatedly simulate Brownian bridge 
#' for (i in 1:replicates) {
#'   paths[[i]] <- Brownian_motion_path_sampler(x = x,
#'                                              times = seq(start_time, end_time, 0.01))
#' }
#' # select the points at the specified time q
#' index <- which(seq(start_time, end_time, 0.01)==end_time)
#' simulated_points <- sapply(1:replicates, function(i) paths[[i]]['X', index])
#' # calculate the theoretical mean and standard deviation of the simulated points at time q
#' theoretical_mean <- x
#' theoretical_sd <- sqrt(end_time-start_time)
#' # plot distribution of the simulated points and the theoretical distribution
#' plot(density(simulated_points))
#' curve(dnorm(x, theoretical_mean, theoretical_sd), add = T, col = 'red')
#' print(paste('Theoretical variance is', end_time-start_time, 'and sample variance is', var(simulated_points)))
#' @export
Brownian_motion_path_sampler <- function(x, times) {
  if (length(times) < 1) {
    stop("Brownian_motion_path_sampler: times must be a vector of length at least 1")
  }
  X <- c(x, rep(NA, length(times)-1))
  time <- sort(times)
  for (i in 2:length(time)) {
    X[i] <- rnorm(1, mean = X[i-1], sd = sqrt(time[i]-time[i-1]))
  }
  return(rbind(X, time))
}
