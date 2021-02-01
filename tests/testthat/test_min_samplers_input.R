context("Input checking for Brownian bridge minimum and Bessel bridge with minimum path samplers")

test_that("Brownian bridge minimum point sampler returns errors if arguments are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  low_bound <- min(x,y) - 2
  up_bound <- min(x,y) - 1
  times <- sort(runif(2, 0, 10))
  # should have error if low_bound > up_bound (must have low_bound < up_bound <= min(x,y))
  expect_error(min_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[2],
                           low_bound = up_bound,
                           up_bound = low_bound),
               "layeredBB::min_sampler: low_bound > up_bound. Must have low_bound < up_bound <= min(x,y)",
               fixed = TRUE)
  # should have error up_bound > min(x,y) (must have low_bound < up_bound <= min(x,y))
  expect_error(min_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[2],
                           low_bound = low_bound,
                           up_bound = max(x,y)),
               "layeredBB::min_sampler: up_bound > min(x,y). Must have low_bound < up_bound <= min(x,y)",
               fixed = TRUE)
  # should have error if s >= t (must have s < t)
  # expect error if s > t
  expect_error(min_sampler(x = x,
                           y = y,
                           s = times[2],
                           t = times[1],
                           low_bound = low_bound,
                           up_bound = up_bound),
               "layeredBB::min_sampler: t <= s. Must have s < t")
  # expect error if t == s
  expect_error(min_sampler(x = x,
                           y = y,
                           s = times[1],
                           t = times[1],
                           low_bound = low_bound,
                           up_bound = up_bound),
               "layeredBB::min_sampler: t <= s. Must have s < t")
  expect_error(min_sampler(x = x,
                           y = y,
                           s = times[2],
                           t = times[2],
                           low_bound = low_bound,
                           up_bound = up_bound),
               "layeredBB::min_sampler: t <= s. Must have s < t")
})

test_that("Bessel bridge with minimum sampler returns errors if arguments are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  m <- min(x,y)-1
  tau <- runif(1, times[1], times[2])
  q <- runif(1, times[1], times[2])
  # should have error if s >= t
  # expect error if s > t
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[2],
                                         t = times[1],
                                         m = m,
                                         tau = tau,
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: t <= s. Must have s < t")
  # expect error if t==s
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[1],
                                         m = m,
                                         tau = tau,
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: t <= s. Must have s < t")
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[2],
                                         t = times[2],
                                         m = m,
                                         tau = tau,
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: t <= s. Must have s < t")
  # should have error if q is outside [s,t]
  # expect error if q < s
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = tau,
                                         q = times[1]-1),
               "layeredBB::min_Bessel_bridge_sampler: requested simulation time q < s")
  # expect error if q > t
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = tau,
                                         q = times[2]+1),
               "layeredBB::min_Bessel_bridge_sampler: requested simulation time q > t")
  # should have error if minimum time tau is outside [s,t]
  # expect error if tau < s
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[1]-1,
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: time of minimum tau < s")
  # expect error if tau > t
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[2]+1,
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: time of minimum tau > t")
  # if tau == s or tau == t, should have error if min != x or min != y
  # expect error if tau == s and min != x
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[1],
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: tau == s and minimum point m != x")
  # expect error if tau == t and min != y
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = m,
                                         tau = times[2],
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: tau == t and minimum point m != y")
  # should have error if minimum is not less than or equal to min(x,y)
  # expect error if min > min(x,y)
  expect_error(min_Bessel_bridge_sampler(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         m = max(x,y),
                                         tau = tau,
                                         q = q),
               "layeredBB::min_Bessel_bridge_sampler: m > min(x,y). Must have m <= min(x,y)",
               fixed = TRUE)
})

test_that("Bessel bridge with minimum path sampler returns errors if arguents are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  m <- min(x,y)-1
  tau <- runif(1, times[1], times[2])
  sim_times <- seq(times[1], times[2], 0.01)
  # should have error if s >= t
  # expect error if s > t
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                         y = y,
                                         s = times[2],
                                         t = times[1],
                                         m = m,
                                         tau = tau,
                                         times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: t <= s. Must have s < t")
  # expect errror if s == t
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[1],
                                              m = m,
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: t <= s. Must have s < t")
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[2],
                                              t = times[2],
                                              m = m,
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: t <= s. Must have s < t")
  # should have error if any of the requested times are outside [s,t]
  # expect error if min(times) < s
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = tau,
                                              times = c(times[1]-1, sim_times)),
               "layeredBB::min_Bessel_bridge_path_sampler: minimum of specified times is less than s")
  # expect error if max(times) > t
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = tau,
                                              times = c(times[2]+1, sim_times)),
               "layeredBB::min_Bessel_bridge_path_sampler: maximum of specified times is greater than t")
  # should have error if minimum time tau is outside [s,t]
  # expect error if tau < s
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[1]-1,
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: time of minimum tau < s")
  # expect error if tau > t
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[2]+1,
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: time of minimum tau > t")
  # if tau == s or tau == t, should have error if min != x or min != y 
  # expect error if tau == s and min != x
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[1],
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: tau == s and minimum point m != x")
  # expect error if tau == t and min != y
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = m,
                                              tau = times[2],
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler:: tau == t and minimum point m != y")
  # should have error if minimum is not less than min(x,y)
  # expect error if min > min(x,y)
  expect_error(min_Bessel_bridge_path_sampler(x = x,
                                              y = y,
                                              s = times[1],
                                              t = times[2],
                                              m = max(x,y),
                                              tau = tau,
                                              times = sim_times),
               "layeredBB::min_Bessel_bridge_path_sampler: m > min(x,y). Must have m <= min(x,y)",
               fixed = TRUE)
})
