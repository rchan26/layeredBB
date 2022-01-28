context("Output checking for Brownian bridge path samplers")

test_that("Brownian bridge path sampler returns expected results", {
  # should return a list of length 2
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  sim_times <- seq(times[1]+0.01, times[2]-0.01, 0.01)
  returned <- Brownian_bridge_path_sampler(x = x,
                                           y = y,
                                           s = times[1],
                                           t = times[2],
                                           times = sim_times)
  # should return list of length 2
  expect_is(returned, "list")
  expect_equal(length(returned), 2)
  expect_equal(names(returned), c("full_path", "simulated_path"))
  # first item in list should be a numeric matrix with 2 rows and 
  # number of columns = the number of unique times passed into function
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$full_path), c(2, length(unique(c(times, sim_times)))))
  expect_false(any(is.na(returned$full_path), is.nan(returned$full_path), is.infinite(returned$full_path)))
  # second item in list should be a numeric matrix with 2 rows and
  # number of columns = length(times)
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$simulated_path), c(2, length(sim_times)))
  expect_false(any(is.na(returned$simulated_path), is.nan(returned$simulated_path), is.infinite(returned$simulated_path)))
  expect_equal(returned$simulated_path[2,], sim_times)
})

test_that("Multivariate Brownian bridge path sampler returns expected results", {
  # should return a list of length 2
  dim <- sample(x = 1:10, size = 1)
  x <- rnorm(dim, 0, 10)
  y <- rnorm(dim, 0, 10)
  times <- sort(runif(2, 0, 10))
  sim_times <- seq(times[1]+0.01, times[2]-0.01, 0.01)
  returned <- multi_brownian_bridge(dim = dim,
                                    x = x,
                                    y = y,
                                    s = times[1],
                                    t = times[2],
                                    times = sim_times)
  # should return list of length 2
  expect_is(returned, "list")
  expect_equal(length(returned), 2)
  expect_equal(names(returned), c("full_path", "simulated_path"))
  # first item in list should be a numeric matrix with dim+1 rows and
  # number of columns = the number of unique times passed into function
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$full_path), c(dim+1, length(unique(c(times, sim_times)))))
  expect_false(any(is.na(returned$full_path), is.nan(returned$full_path), is.infinite(returned$full_path)))
  # second item in list should be a numeric matrix with dim+1 rows and
  # number of columns = length(sim_times)
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$simulated_path), c(dim+1, length(sim_times)))
  expect_false(any(is.na(returned$simulated_path), is.nan(returned$simulated_path), is.infinite(returned$simulated_path)))
  expect_equal(returned$simulated_path[dim+1,], sim_times)
})
