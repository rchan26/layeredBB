context("Input checking for Brownian bridge maximum and Bessel bridge with maximum path samplers")

test_that("Brownian bridge maximum point sampler returns errors if arguments are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  low_bound <- max(x,y) + 1
  up_bound <- max(x,y) + 2
  times <- sort(runif(2, 0, 10))
  # should have error if low_bound >= up_bound (must have max(x,y) <= low_bound < up_bound)
  # expect error if low_bound > up_bound
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[2],
                           low_bound = up_bound,
                           up_bound = low_bound),
               "layeredBB::max_sampler: low_bound >= up_bound. Must have max(x,y) <= low_bound < up_bound",
               fixed = TRUE)
  # expect error if low_bound == up_bound
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[2],
                           low_bound = low_bound,
                           up_bound = low_bound),
               "layeredBB::max_sampler: low_bound >= up_bound. Must have max(x,y) <= low_bound < up_bound",
               fixed = TRUE)
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[2],
                           low_bound = up_bound,
                           up_bound = up_bound),
               "layeredBB::max_sampler: low_bound >= up_bound. Must have max(x,y) <= low_bound < up_bound",
               fixed = TRUE)
  # should have error if low_bound < max(x,y) (must have max(x,y) <= low_bound < up_bound)
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[2],
                           low_bound = min(x,y),
                           up_bound = up_bound),
               "layeredBB::max_sampler: low_bound < max(x,y). Must have max(x,y) <= low_bound < up_bound",
               fixed = TRUE)
  # should have error if s >= t (must have s < t)
  # expect error if s > t
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[2],
                           t = times[1],
                           low_bound = low_bound,
                           up_bound = up_bound),
               "layeredBB::max_sampler: t <= s. Must have s < t")
  # expect error if t == s
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[1],
                           low_bound = low_bound,
                           up_bound = up_bound),
               "layeredBB::max_sampler: t <= s. Must have s < t")
  expect_error(max_sampler(x = x,
                           y = y,
                           s = times[2],
                           t = times[2],
                           low_bound = low_bound,
                           up_bound = up_bound),
               "layeredBB::max_sampler: t <= s. Must have s < t")
})

test_that("Bessel bridge with maximum sampler returns errors if arguments are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  m <- max(x,y)+1
  tau <- runif(1, times[1], times[2])
  q <- runif(1, times[1], times[2])
  # should have error if s >= t
  # expect error if s > t
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[2],
                                         t = times[1],
                                         m = m,
                                         tau = tau,
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: t <= s. Must have s < t")
  # expect error if s == t
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[1],
                                         m = m,
                                         tau = tau,
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: t <= s. Must have s < t")
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[2],
                                         t = times[2],
                                         m = m,
                                         tau = tau,
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: t <= s. Must have s < t")
  # should have error if q is oustide [s,t]
  # expect error if q < s
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = tau,
                                         q = times[1]-1),
               "layeredBB::max_Bessel_bridge_sampler: requested simulation time q < s")
  # expect error if q > t
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = tau,
                                         q = times[2]+1),
               "layeredBB::max_Bessel_bridge_sampler: requested simulation time q > t")
  # should have error if tau is oustide [s,t]
  # expect error if tau < s
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[1]-1,
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: time of maximum tau < s")
  # expect error if tau > t
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[2]+1,
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: time of maximum tau > t")
  # if tau == s or tau == t, should have error if max != x or max != y
  # expect error if tau == s and max != x
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[1],
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: tau == s and maximum point m != x")
  # expect error if tau == t and max != y
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[2],
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler:: tau == t and maximum point m != y")
  # should have error if max < max(x,y)
  # expect error ix max < max(x,y)
  expect_error(max_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = min(x,y),
                                         tau = tau,
                                         q = q),
               "layeredBB::max_Bessel_bridge_sampler: m < max(x,y). Must have m >= max(x,y)",
               fixed = TRUE)
})

test_that("Bessel bridge with maximum path sampler returns errors if arguents are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  m <- max(x,y)+1
  tau <- runif(1, times[1], times[2])
  sim_times <- seq(times[1], times[2], 0.01)
  # should have error if s >= t
  # expect error if s > t
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[2],
                                              t = times[1],
                                              m = m,
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: t <= s. Must have s < t")
  # expect error if s == t
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[1],
                                              m = m,
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: t <= s. Must have s < t")
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[2],
                                              t = times[2],
                                              m = m,
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: t <= s. Must have s < t")
  # should have error if any of the requested times are outside [s,t]
  # expect error if min(times) < s
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = tau,
                                              times = c(times[1]-1, sim_times)),
               "layeredBB::max_Bessel_bridge_path_sampler: minimum of specified times is less than s")
  # expect error if max(times) > t
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = tau,
                                              times = c(times[2]+1, sim_times)),
               "layeredBB::max_Bessel_bridge_path_sampler: maximum of specified times is greater than t")
  # should have error if tau is outside [s,t]
  # expect error if tau < s
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[1]-1,
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: time of maximum tau < s")
  # expect error if tau > t
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[2]+1,
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: time of maximum tau > t")
  # if tau == s or tau == t, should have error if max != x or max != y
  # expect error if tau == s and max != x
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[1],
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: tau == s and maximum point m != x")
  # expect error if tau == t and max != y
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[2],
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler:: tau == t and maximum point m != y")
  # should have error if max < max(x,y)
  # expect error if max < max(x,y)
  expect_error(max_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = min(x,y),
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::max_Bessel_bridge_path_sampler: m < max(x,y). Must have m >= max(x,y)",
               fixed = TRUE)
})
