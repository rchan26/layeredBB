context("Checking functions to flip Gamma and Delta probabilities using Cauchy sums and retrospective Bernoulli sampling")

test_that("function to calculate product of vector returns expected results", {
  expect_equal(product_vector(c(1,2,3,4,5)), 120)
  expect_true(is.na(product_vector(c(1,2,3,4,5, NA))))
  expect_true(is.nan(product_vector(c(1,2,3,4,5, NaN))))
  expect_true(is.infinite(product_vector(c(1,2,3,4,5, Inf))))
})

test_that("Gamma coin flipper returns errors if arguments are outside constraints", {
  u <- runif(1, 0, 1)
  k <- sample(x = 1:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  l <- min(x,y) - 1
  v <- max(x,y) + 1
  # should have error if u not in [0,1]
  # expect error if u < 0
  expect_error(gamma_coin(u = -u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[1],
                          t = times[2],
                          l = l,
                          v = v),
               "layeredBB::gamma_coin: u must be in interval [0,1]",
               fixed = TRUE)
  # expect error if u > 1
  expect_error(gamma_coin(u = u+1,
                          k = k,
                          x = x,
                          y = y,
                          s = times[1],
                          t = times[2],
                          l = l,
                          v = v),
               "layeredBB::gamma_coin: u must be in interval [0,1]",
               fixed = TRUE)
  # should have error if s >= t
  # expect error if s > t
  expect_error(gamma_coin(u = u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[2],
                          t = times[1],
                          l = l,
                          v = v),
               "layeredBB::gamma_coin: s >= t. Must have s < t")
  # expect error if s == t
  expect_error(gamma_coin(u = u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[1],
                          t = times[1],
                          l = l,
                          v = v),
               "layeredBB::gamma_coin: s >= t. Must have s < t")
  expect_error(gamma_coin(u = u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[2],
                          t = times[2],
                          l = l,
                          v = v),
               "layeredBB::gamma_coin: s >= t. Must have s < t")
})

test_that("Gamma coin flipper returns expected results", {
  u <- runif(1, 0, 1)
  k <- sample(x = 1:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  l <- min(x,y) - 1
  v <- max(x,y) + 1
  returned <- gamma_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         l = l,
                         v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if l == min(x,y)
  returned <- gamma_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         l = min(x,y),
                         v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if v == max(x,y)
  returned <- gamma_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         l = l,
                         v = max(x,y))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if l == min(x,y) and v = max(x,y)
  returned <- gamma_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         l = min(x,y),
                         v = max(x,y))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # if min(x,y) < l (already outside layer), should be FALSE
  # if x < l
  returned <- gamma_coin(u = u,
                         k = k,
                         x = l-1,
                         y = y,
                         s = times[1],
                         t = times[2],
                         l = l,
                         v = v)
  expect_equal(returned, FALSE)
  # if y < l
  returned <- gamma_coin(u = u,
                         k = k,
                         x = x,
                         y = l-1,
                         s = times[1],
                         t = times[2],
                         l = l,
                         v = v)
  expect_equal(returned, FALSE)
  # if max(x,y) > v (already outside layer), should be FALSE
  # if x > v
  returned <- gamma_coin(u = u,
                         k = k,
                         x = v+1,
                         y = y,
                         s = times[1],
                         t = times[2],
                         l = l,
                         v = v)
  expect_equal(returned, FALSE)
  # if y > v
  returned <- gamma_coin(u = u,
                         k = k,
                         x = x,
                         y = v+1,
                         s = times[1],
                         t = times[2],
                         l = l,
                         v = v)
  expect_equal(returned, FALSE)
})

test_that("Delta coin flipper returns errors if arguments are outside constraints", {
  u <- runif(1, 0, 1)
  k <- sample(x = 1:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  min <- min(x,y) - 1
  v <- max(x,y) + 1
  # should have error if u not in [0,1]
  # expect error if u < 0
  expect_error(delta_coin(u = -u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[1],
                          t = times[2],
                          min = min,
                          v = v),
               "layeredBB::delta_coin: u must be in interval [0,1]",
               fixed = TRUE)
  # expect error if u > 1
  expect_error(delta_coin(u = u+1,
                          k = k,
                          x = x,
                          y = y,
                          s = times[1],
                          t = times[2],
                          min = min,
                          v = v),
               "layeredBB::delta_coin: u must be in interval [0,1]",
               fixed = TRUE)
  # should have error if s >= t
  # expect error if s > t
  expect_error(delta_coin(u = u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[2],
                          t = times[1],
                          min = min,
                          v = v),
               "layeredBB::delta_coin: s >= t. Must have s < t")
  # expect error if s == t
  expect_error(delta_coin(u = u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[1],
                          t = times[1],
                          min = min,
                          v = v),
               "layeredBB::delta_coin: s >= t. Must have s < t")
  expect_error(delta_coin(u = u,
                          k = k,
                          x = x,
                          y = y,
                          s = times[2],
                          t = times[2],
                          min = min,
                          v = v),
               "layeredBB::delta_coin: s >= t. Must have s < t")
})

test_that("Delta coin flipper returns expected results", {
  u <- runif(1, 0, 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  min <- min(x,y) - 1
  v <- max(x,y) + 1
  # calculate k
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- delta_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min,
                         v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if min == min(x,y) too
  min <- min(x,y)
  # calculate k
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- delta_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min(x,y),
                         v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if v == max(x,y) and min < min(x,y)
  min <- min(x,y) - 1
  v <- max(x,y)
  # calculate k
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- delta_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min,
                         v = v)
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # should work if min == min(x,y) and v == max(x,y)
  min <- min(x,y)
  # calculate k
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- delta_coin(u = u,
                         k = k,
                         x = x,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min(x,y),
                         v = max(x,y))
  expect_is(returned, "logical")
  expect_equal(length(returned), 1)
  # if min(x,y) < min (already outside layer), should be FALSE
  # if x < min
  returned <- delta_coin(u = u,
                         k = k,
                         x = min-1,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min,
                         v = v)
  expect_equal(returned, FALSE)
  # if y < min
  returned <- delta_coin(u = u,
                         k = k,
                         x = min-1,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min,
                         v = v)
  expect_equal(returned, FALSE)
  # if max(x,y) > v (already outside layer), should be FALSE
  # if x > v
  returned <- delta_coin(u = u,
                         k = k,
                         x = v+1,
                         y = y,
                         s = times[1],
                         t = times[2],
                         min = min,
                         v = v)
  expect_equal(returned, FALSE)
  # if y > v
  returned <- delta_coin(u = u,
                         k = k,
                         x = x,
                         y = v+1,
                         s = times[1],
                         t = times[2],
                         min = min,
                         v = v)
  expect_equal(returned, FALSE)
})
