context("Output checking for Brownian bridge minimum and Bessel bridge with minimum path samplers")

test_that("Brownian bridge minimum point sampler returns expected results", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  low_bound <- min(x,y) - 2
  up_bound <- min(x,y) - 1
  times <- sort(runif(2, 0, 10))
  # if low_bound < up_bound < min(x,y)
  returned <- min_sampler(x = x,
                          y = y,
                          s = times[1],
                          t = times[2],
                          low_bound = low_bound,
                          up_bound = up_bound)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_equal(names(returned), c("min", "tau"))
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_true(returned['min'] >= low_bound & returned['min'] <= up_bound)
  expect_true(returned['tau'] >= times[1] & returned['tau'] <= times[2])
  # if low_bound < up_bound == min(x,y)
  returned <- min_sampler(x = x,
                          y = y,
                          s = times[1],
                          t = times[2],
                          low_bound = low_bound,
                          up_bound = min(x,y))
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_equal(names(returned), c("min", "tau"))
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_true(returned['min'] >= low_bound & returned['min'] <= min(x,y))
  expect_true(returned['tau'] >= times[1] & returned['tau'] <= times[2])
})

test_that("Bessel bridge with minimum sampler returns expected results", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  m <- min(x,y)-1
  tau <- runif(1, times[1], times[2])
  q <- runif(1, times[1], times[2])
  # if q in (s,t) and q != tau
  returned <- min_Bessel_bridge_sampler(x = x,
                                        y = y,
                                        s = times[1],
                                        t = times[2],
                                        m = m,
                                        tau = tau,
                                        q = q)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_true(returned >= m)
  # if q == s, should return x
  returned <- min_Bessel_bridge_sampler(x = x,
                                        y = y,
                                        s = times[1],
                                        t = times[2],
                                        m = m,
                                        tau = tau,
                                        q = times[1])
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned, x)
  # if q == t, should return y
  returned <- min_Bessel_bridge_sampler(x = x,
                                        y = y,
                                        s = times[1],
                                        t = times[2],
                                        m = m,
                                        tau = tau,
                                        q = times[2])
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned, y)
  # if q == tau, should return m
  returned <- min_Bessel_bridge_sampler(x = x,
                                        y = y,
                                        s = times[1],
                                        t = times[2],
                                        m = m,
                                        tau = tau,
                                        q = tau)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned, m)
  # if tau == s and m == x (or tau == t and m == y), should still work
  if (min(x,y) == x) {
    tau <- times[1]
  } else {
    tau <- times[2]
  }
  returned <- min_Bessel_bridge_sampler(x = x,
                                        y = y,
                                        s = times[1],
                                        t = times[2],
                                        m = min(x,y),
                                        tau = tau,
                                        q = q)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_true(returned >= m)
})

test_that("Bessel bridge with minimum path sampler retruns expected results", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  m <- min(x,y)-1
  tau <- runif(1, times[1], times[2])
  sim_times <- seq(times[1]+0.01, times[2]-0.01, 0.01)
  returned <- min_Bessel_bridge_path_sampler(x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             m = m,
                                             tau = tau,
                                             sim_times)
  expect_is(returned, "list")
  expect_equal(length(returned), 3)
  expect_equal(names(returned), c("full_path", "simulated_path", "remove_m_path"))
  # first item should be a numeric matrix with 2 rows and
  # number of columns = the number of unique times passed into the function
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$full_path), c(2, length(unique(c(times, sim_times, tau)))))
  expect_false(any(is.na(returned$full_path), is.nan(returned$full_path), is.infinite(returned$full_path)))
  expect_true(all(returned$full_path[1,] >= m))
  # second item in list should be a numeric matrix with 2 rows and
  # number of columns = length(times)
  expect_is(returned$simulated_path, "matrix")
  expect_equal(dim(returned$simulated_path), c(2, length(sim_times)))
  expect_false(any(is.na(returned$simulated_path), is.nan(returned$simulated_path), is.infinite(returned$simulated_path)))
  expect_true(all(returned$simulated_path[1,] >= m))
  expect_equal(returned$simulated_path[2,], sim_times)
  # third item in the list should be a numeric matrix with 2 rows and
  # number of columns = the number of unique times passed into the function minus 1 (removed minimum)
  expect_is(returned$remove_m_path, "matrix")
  expect_equal(dim(returned$remove_m_path), c(2, length(unique(c(times, sim_times)))))
  expect_false(any(is.na(returned$remove_m_path), is.nan(returned$remove_m_path), is.infinite(returned$remove_m_path)))
  expect_true(all(returned$remove_m_path[1,] >= m))
})
