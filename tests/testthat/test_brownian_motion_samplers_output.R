context("Output checking for Brownian motion path samplers")

test_that("Brownian motion path sampler returns expected results", {
  x <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  returned <- Brownian_motion_path_sampler(x = x, times = times)
  # should return numeric matrix with 2 rows and 
  # number of columns = length(times)
  expect_is(returned, "matrix")
  expect_equal(dim(returned), c(2, length(unique(times))))
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if times has some duplicate values and is unsorted, we should sort and remove duplicates
  times <- c(0, 1, 2, 3, 1, 5, 6, 6, 8)
  returned <- Brownian_motion_path_sampler(x = x, times = times)
  # should return numeric matrix with 2 rows and 
  # number of columns = length(times)
  expect_is(returned, "matrix")
  expect_equal(dim(returned), c(2, length(unique(times))))
  expect_equal(returned['time', ], sort(unique(times)))
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
})
