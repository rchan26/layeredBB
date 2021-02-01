context("Input checking for Brownian bridge path samplers")

test_that("Brownian bridge path sampler returns expected results", {
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
  returned <- layered_brownian_bridge(x = x,
                                      y = y,
                                      s = times[1],
                                      t = times[2],
                                      bessel_layer = bes_layer,
                                      times = sim_times)
  expect_is(returned, "list")
  expect_equal(length(returned), 3)
  expect_equal(names(returned), c("full_path", "simulated_path", "remove_m_path"))
  # first item should be a numeric matrix with 2 rows and
  # number of columns = the number of unique times passed into the function
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$full_path), c(2, length(unique(c(times, sim_times)))+1))
  expect_false(any(is.na(returned$full_path), is.nan(returned$full_path), is.infinite(returned$full_path)))
  expect_true(all(returned$full_path[1,] >= bes_layer$L))
  expect_true(all(returned$full_path[1,] <= bes_layer$U))
  # second item in list should be a numeric matrix with 2 rows and
  # number of columns = length(times)
  expect_is(returned$simulated_path, "matrix")
  expect_equal(dim(returned$simulated_path), c(2, length(sim_times)))
  expect_false(any(is.na(returned$simulated_path), is.nan(returned$simulated_path), is.infinite(returned$simulated_path)))
  expect_true(all(returned$simulated_path[1,] >= bes_layer$L))
  expect_true(all(returned$simulated_path[1,] <= bes_layer$U))
  expect_equal(returned$simulated_path[2,], sim_times)
  # third item in the list should be a numeric matrix with 2 rows and
  # number of columns = the number of unique times passed into the function minus 1 (removed maximum)
  expect_is(returned$remove_m_path, "matrix")
  expect_equal(dim(returned$remove_m_path), c(2, length(unique(c(times, sim_times)))))
  expect_false(any(is.na(returned$remove_m_path), is.nan(returned$remove_m_path), is.infinite(returned$remove_m_path)))
  expect_true(all(returned$remove_m_path[1,] >= bes_layer$L))
  expect_true(all(returned$remove_m_path[1,] <= bes_layer$U))
})

test_that("Multivariate Brownian bridge path sampler returns expected results", {
  dim <- sample(x = 1:10, size = 1)
  x <- rnorm(dim, 0, 10)
  y <- rnorm(dim, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  bes_layer <- multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             mult = mult)
  sim_times <- seq(times[1], times[2], 0.01)
  returned <- multi_layered_brownian_bridge(dim = dim,
                                            x = x,
                                            y = y,
                                            s = times[1],
                                            t = times[2],
                                            bessel_layer = bes_layer,
                                            times = sim_times)
  expect_is(returned, "list")
  expect_equal(length(returned), 2)
  expect_equal(names(returned), c("full_path", "simulated_path"))
  # first item should be a numeric matrix with 2 rows and
  # number of columns = the number of unique times passed into the function
  expect_is(returned$full_path, "matrix")
  expect_equal(dim(returned$full_path), c(dim+1, length(unique(c(times, sim_times)))))
  expect_false(any(is.na(returned$full_path), is.nan(returned$full_path), is.infinite(returned$full_path)))
  for (i in 1:dim) {
    expect_true(all(returned$full_path[i,] >= bes_layer[[i]]$L))
    expect_true(all(returned$full_path[i,] <= bes_layer[[i]]$U))
  }
  # second item in list should be a numeric matrix with 2 rows and
  # number of columns = length(times)
  expect_is(returned$simulated_path, "matrix")
  expect_equal(dim(returned$simulated_path), c(dim+1, length(sim_times)))
  expect_false(any(is.na(returned$simulated_path), is.nan(returned$simulated_path), is.infinite(returned$simulated_path)))
  for (i in 1:dim) {
    expect_true(all(returned$simulated_path[i,] >= bes_layer[[i]]$L))
    expect_true(all(returned$simulated_path[i,] <= bes_layer[[i]]$U))
  }
  expect_true(all(returned$simulated_path[dim+1,] == sim_times))
})
