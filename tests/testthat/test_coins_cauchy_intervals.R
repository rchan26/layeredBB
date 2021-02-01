context("Checking functions to flip Gamma and Delta probabilities for intervals using Cauchy sums and retrospective Bernoulli sampling")

test_that("Gamma coin flipper for intervals returns errors if arguments are outside constraints", {
  # setting up Brownian bridge path
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  sim_times <- seq(times[1], times[2], 0.1)
  BB <- Brownian_bridge_path_sampler(x = x, y = y, s = times[1], t = times[2], times = sim_times)$full_path
  # setting up input variables for gamma_coin_intervals
  u <- runif(1, 0, 1)
  X <- BB[1,]
  times <- BB[2,]
  l <- min(BB[1,])-1
  v <- max(BB[1,])+1
  k <- sample(x = 1:10, size = 1)
  # should have error if X and times are not of the same size
  # expect error if length(X) < length(times)
  expect_error(gamma_coin_intervals(u = u,
                                    k = k,
                                    X = X[1:(nrow(BB)-1)],
                                    times = times,
                                    l = l,
                                    v = v),
               "layeredBB::gamma_coin_intervals: vector lengths of X and times are not equal")
  # expect error if length(X) > length(times)
  expect_error(gamma_coin_intervals(u = u,
                                    k = k,
                                    X = X,
                                    times = times[1:(nrow(BB)-1)],
                                    l = l,
                                    v = v),
               "layeredBB::gamma_coin_intervals: vector lengths of X and times are not equal")
  # should have error if u not in [0,1]
  # expect error if u < 0
  expect_error(gamma_coin_intervals(u = -u,
                                    k = k,
                                    X = X,
                                    times = times,
                                    l = l,
                                    v = v),
               "layeredBB::gamma_coin_intervals: u must be in interval [0,1]",
               fixed = TRUE)
  # expect error if u > 1
  expect_error(gamma_coin_intervals(u = u+1,
                                    k = k,
                                    X = X,
                                    times = times,
                                    l = l,
                                    v = v),
               "layeredBB::gamma_coin_intervals: u must be in interval [0,1]",
               fixed = TRUE)
})

test_that("Gamma coin flipper for intervals returns expected results", {
  # setting up Brownian bridge path
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  sim_times <- seq(times[1], times[2], 0.1)
  BB <- Brownian_bridge_path_sampler(x = x, y = y, s = times[1], t = times[2], times = sim_times)$full_path
  # setting up input variables for gamma_coin_intervals
  u <- runif(1, 0, 1)
  X <- BB[1,]
  times <- BB[2,]
  l <- min(BB[1,])-1
  v <- max(BB[1,])+1
  k <- sample(x = 1:10, size = 1)
  returned <- gamma_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   l = l,
                                   v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if l == min(X)
  returned <- gamma_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   l = min(X),
                                   v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if v == max(X)
  returned <- gamma_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   l = l,
                                   v = max(X))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if l == min(X) and v == max(X)
  returned <- gamma_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   l = min(X),
                                   v = max(X))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # if min(X) < l (already outside layer), should be FALSE
  returned <- gamma_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   l = min(X)+1,
                                   v = v)
  expect_equal(returned, FALSE)
  # if max(X) > v (already outside layer), should be TRUE
  returned <- gamma_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   l = l,
                                   v = max(X)-1)
  expect_equal(returned, FALSE)
})

test_that("Delta coin flipper for intervals returns errors if arguments are outside constraints", {
  # setting up Brownian bridge path
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  sim_times <- seq(times[1], times[2], 0.1)
  BB <- Brownian_bridge_path_sampler(x = x, y = y, s = times[1], t = times[2], times = sim_times)$full_path
  # setting up input variables to delta_coin_intervals
  u <- runif(1, 0, 1)
  X <- BB[1,]
  times <- BB[2,]
  min <- min(BB[1,])-1
  v <- max(BB[1,])+1
  k <- sample(x = 1:10, size = 1)
  # should have error if X and times are not of the same size
  # expect error if length(X) < length(times)
  expect_error(delta_coin_intervals(u = u,
                                    k = k,
                                    X = X[1:(nrow(BB)-1)],
                                    times = times,
                                    min = min,
                                    v = v),
               "layeredBB::delta_coin_intervals: vector lengths of X and times are not equal")
  # expect error if length(X) > length(times)
  expect_error(delta_coin_intervals(u = u,
                                    k = k,
                                    X = X,
                                    times = times[1:(nrow(BB)-1)],
                                    min = min,
                                    v = v),
               "layeredBB::delta_coin_intervals: vector lengths of X and times are not equal")
  # should have error if u not in [0,1]
  # expect error if u < 0
  expect_error(delta_coin_intervals(u = -u,
                                    k = k,
                                    X = X,
                                    times = times,
                                    min = min,
                                    v = v),
               "layeredBB::delta_coin_intervals: u must be in interval [0,1]",
               fixed = TRUE)
  # expect error if u > 1
  expect_error(delta_coin_intervals(u = u+1,
                                    k = k,
                                    X = X,
                                    times = times,
                                    min = min,
                                    v = v),
               "layeredBB::delta_coin_intervals: u must be in interval [0,1]",
               fixed = TRUE)
})

test_that("Delta coin flipper for intervals returns expected results", {
  u <- runif(1, 0, 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  sim_times <- seq(times[1], times[2], 0.1)
  BB <- Brownian_bridge_path_sampler(x = x, y = y, s = times[1], t = times[2], times = sim_times)$full_path
  X <- BB[1,]
  times <- BB[2,]
  min <- min(BB[1,])-1
  v <- max(BB[1,])+1
  # calculate k
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- delta_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   min = min,
                                   v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if min == min(X)
  returned <- delta_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   min = min(X),
                                   v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if v =- max(X)
  returned <- delta_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   min = min,
                                   v = max(X))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if min == min(X) and v == max(X)
  returned <- delta_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   min = min(X),
                                   v = max(X))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # if min(X) < min (already outside layer), should be FALSE
  returned <- delta_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   min = min(X)+1,
                                   v = v)
  # if max(X) > v (already outside layer), should be TRUE
  returned <- delta_coin_intervals(u = u,
                                   k = k,
                                   X = X,
                                   times = times,
                                   min = min,
                                   v = max(X)-1)
})
