context("Input checking for Brownian motion path samplers")

test_that("Brownian motion path sampler returns errors if arguments are outside constraints", {
  x <- rnorm(1, 0, 10)
  # should have error if length(times) < 2
  times <- runif(1)
  expect_error(Brownian_motion_path_sampler(x, times),
               "layeredBB::Brownian_motion_path_sampler: times must be a vector of length at least 2")
  # should have error if unique(length(times) < 2)
  times <- rep(runif(1), 2)
  expect_error(Brownian_motion_path_sampler(x, times),
               "layeredBB::Brownian_motion_path_sampler: times must be a vector with at least 2 unique values in")
})
