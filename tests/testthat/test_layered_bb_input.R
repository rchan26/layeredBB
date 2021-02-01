context("Input checking for Layered Brownian bridge path samplers")

test_that("Layered Brownian bridge path sampler returns errors if arguments are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  bes_layer <- bessel_layer_simulation(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[2],
                                       mult = mult)
  sim_times <- seq(times[1], times[2], 0.01)
  # should have error if s >= t
  # expect error if s > t
  expect_error(layered_brownian_bridge(x = x,
                                       y = y,
                                       s = times[2],
                                       t = times[1],
                                       bessel_layer = bes_layer,
                                       times = sim_times),
               "layeredBB::layered_brownian_bridge: t <= s. Must have s < t")
  # expect error if s == t
  expect_error(layered_brownian_bridge(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[1],
                                       bessel_layer = bes_layer,
                                       times = sim_times),
               "layeredBB::layered_brownian_bridge: t <= s. Must have s < t")
  expect_error(layered_brownian_bridge(x = x,
                                       y = y,
                                       s = times[2],
                                       t = times[2],
                                       bessel_layer = bes_layer,
                                       times = sim_times),
               "layeredBB::layered_brownian_bridge: t <= s. Must have s < t")
  # should have errors if any of the times fall outside of [s,t]
  # expect error if min(times) < s
  expect_error(layered_brownian_bridge(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[2],
                                       bessel_layer = bes_layer,
                                       times = c(-10, sim_times)),
               "layeredBB::layered_brownian_bridge: minimum of specified times is less than s")
  # expect error if max(times) > t
  expect_error(layered_brownian_bridge(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[2],
                                       bessel_layer = bes_layer,
                                       times = c(20, sim_times)),
               "layeredBB::layered_brownian_bridge: maximum of specified times is greater than t")
})

test_that("Multivariate Brownian bridge path sampler returns errors if arguments are outside constraints", {
  dim <- sample(x = 1:10, size = 1)
  x <- rnorm(dim, 0, 10)
  y <- rnorm(dim, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  bes_layers <- multi_bessel_layer_simulation(dim = dim,
                                              x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              mult = mult)
  sim_times <- seq(times[1], times[2], 0.01)
  # should have error if s >= t
  # expect error if s > t
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[2],
                                             t = times[1],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: t <= s. Must have s < t")
  # expect error if s == t
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[1],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: t <= s. Must have s < t")
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[2],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: t <= s. Must have s < t")
  # should have errors if the start and end points aren't vectors of length dim
  # expect error if x is not vector of length dim
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = rnorm(dim-1),
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: size of x is not equal to dim")
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = rnorm(dim+1),
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: size of x is not equal to dim")
  # expect error if y is not vector of length dim
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = rnorm(dim+1),
                                             s = times[1],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: size of y is not equal to dim")
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = rnorm(dim-1),
                                             s = times[1],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = sim_times),
               "layeredBB::multi_layered_brownian_bridge: size of y is not equal to dim")
  # should have errors if any of the times fall outside of [s,t]
  # expect error if min(times) < s
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = c(times[1]-1, sim_times)),
               "layeredBB::multi_layered_brownian_bridge: minimum of specified times is less than s")
  # expect error if max(times) > t
  expect_error(multi_layered_brownian_bridge(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             bessel_layers = bes_layers,
                                             times = c(times[2]+1, sim_times)),
               "layeredBB::multi_layered_brownian_bridge: maximum of specified times is greater than t")
})