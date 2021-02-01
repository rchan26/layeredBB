context("Checking functions to calculate Cauchy sums")

test_that("Check easigma_bar, easigma, eaphi_bar, eaphi, eapsi, eachi return expected results", {
  j <- sample(x = 0:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  # ----- if l == min(x,y) and v == max(x,y)
  l <- min(x,y)
  v <- max(x,y)
  # easigma_bar
  returned <- easigma_bar(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # easigma
  returned <- easigma(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eaphi_bar
  returned <- eaphi_bar(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eaphi
  returned <- eaphi(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eapsi
  returned <- eapsi(j = j, xoy = max(x,y), s = times[1], t = times[2], min = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eachi
  returned <- eachi(j = j, xoy = max(x,y), s = times[1], t = times[2], min = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # ----- if l < min(x,y) and v > max(x,y)
  l <- min(x,y)-1
  v <- max(x,y)+1
  # easigma_bar
  returned <- easigma_bar(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # easigma
  returned <- easigma(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eaphi_bar
  returned <- eaphi_bar(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eaphi
  returned <- eaphi(j = j, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eapsi
  returned <- eapsi(j = j, xoy = max(x,y), s = times[1], t = times[2], min = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # eachi
  returned <- eachi(j = j, xoy = max(x,y), s = times[1], t = times[2], min = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
})

test_that("Check eagamma returns expected results", {
  n <- sample(x = 1:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  # if l == min(x,y) and v == max(x,y)
  l <- min(x,y)
  v <- max(x,y)
  returned <- eagamma(n = n, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if l < min(x,y) and v > max(x,y)
  l <- min(x,y)-1
  v <- max(x,y)+1
  returned <- eagamma(n = n, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if min(x,y) < l (already outside layer), should return 0
  # expect 0 if x < l
  returned <- eagamma(n = n, x = l-1, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, 0)
  # expect 0 if y < l
  returned <- eagamma(n = n, x = x, y = l-1, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, 0)
  # if max(x,y) > v (already outside layer), should return 0
  # expect 0 if x > v
  returned <- eagamma(n = n, x = v+1, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, 0)
  # expect 0 if y > v
  returned <- eagamma(n = n, x = x, y = v+1, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, 0)
})

test_that("Check eadelta1 returnS expected results", {
  n <- sample(x = 0:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  min <- min(x,y) - 1
  # if v == max(x,y), should work normally
  v <- max(x,y)
  returned <- eadelta1(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if v > max(x,y), should work normally
  v <- max(x,y)+1
  returned <- eadelta1(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if min == min(x,y), should return an NaN (should be using eadelta2 in this case)
  returned <- eadelta1(n = n, x = x, y = y, s = times[1], t = times[2], min = min(x,y), v = v)
  expect_true(is.nan(returned))
  # if min(x,y) < min (already outside layer), should return 0
  # expect 0 if x < min
  returned <- eadelta1(n = n, x = min-1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
  # expect 0 if y < min
  returned <- eadelta1(n = n, x = x, y = min-1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
  # if max(x,y) > v (already outside layer), should return 0
  # expect 0 if x > v
  returned <- eadelta1(n = n, x = v+1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
  # expect 0 if y > v
  returned <- eadelta1(n = n, x = x, y = v+1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
})

test_that("Check that eadelta2 returns expected results", {
  n <- sample(x = 0:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  min <- min(x,y)
  # if v == max(x,y), should work normally
  v <- max(x,y)
  returned <- eadelta2(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if v > max(x,y), should work normally
  v <- max(x,y)+1
  returned <- eadelta2(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 1)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  # if min < min(x,y), should return an NaN (should be using eadelta1 in this case)
  returned <- eadelta2(n = n, x = x, y = y, s = times[1], t = times[2], min = min(x,y)-1, v = v)
  expect_true(is.nan(returned))
  # if min(x,y) < min (already outside layer), should return 0
  # expect 0 if x < min
  returned <- eadelta2(n = n, x = min-1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
  # expect 0 if y < min
  returned <- eadelta2(n = n, x = x, y = min-1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
  # if max(x,y) > v (already outside layer), should return 0
  # expect 0 if x > v
  returned <- eadelta2(n = n, x = v+1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
  # expect 0 if y > v
  returned <- eadelta2(n = n, x = x, y = v+1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, 0)
})

test_that("Check that eadelta returns expected results", {
  n <- sample(x = 0:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  # if min(x,y) > min, eadelta should equal eadelta1
  min <- min(x,y)-1
  # if v == max(x,y)
  v <- max(x,y)
  returned_delta <- eadelta(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta1 <- eadelta1(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta1)
  # if v > max(x,y)
  v <- max(x,y)+1
  returned_delta <- eadelta(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta1 <- eadelta1(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta1)
  # if min(x,y) == min, eadelta should equal eadelta2
  min <- min(x,y)
  # if v == max(x,y)
  v <- max(x,y)
  returned_delta <- eadelta(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta2 <- eadelta2(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta2)
  # if v > max(x,y)
  v <- max(x,y)+1
  returned_delta <- eadelta(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta2 <- eadelta2(n = n, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta2)
})

test_that("Check eagamma_intervals returns expected results", {
  k <- sample(x = 0:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  # if l == min(x,y) and v == max(x,y)
  l <- min(x,y)
  v <- max(x,y)
  returned <- eagamma_intervals(k = k, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  S_2k <- eagamma(n = 2*k, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  S_2k_plus_1 <- eagamma(n = 2*k+1, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_equal(returned[1], S_2k_plus_1)
  expect_equal(returned[2], S_2k)
  # if l < min(x,y) and v > max(x,y)
  l <- min(x,y)-1
  v <- max(x,y)+1
  returned <- eagamma_intervals(k = k, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  S_2k <- eagamma(n = 2*k, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  S_2k_plus_1 <- eagamma(n = 2*k+1, x = x, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_equal(returned[1], S_2k_plus_1)
  expect_equal(returned[2], S_2k)
  # if min(x,y) < l (already outside layer), should return c(0,0)
  # expect c(0,0) if x < l
  returned <- eagamma_intervals(k = k, x = l-1, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, c(0,0))
  # expect c(0,0) if y < l
  returned <- eagamma_intervals(k = k, x = x, y = l-1, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, c(0,0))
  # if max(x,y) > v (already outside layer), should return c(0,0)
  # expect c(0,0) if x > v
  returned <- eagamma_intervals(k = k, x = v+1, y = y, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, c(0,0))
  # expect c(0,0) if y > v
  returned <- eagamma_intervals(k = k, x = x, y = v+1, s = times[1], t = times[2], l = l, v = v)
  expect_equal(returned, c(0,0))
})

test_that("Check eadelta1_intervals returns expected results", {
  k <- sample(x = 0:10, size = 1)
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  min <- min(x,y) - 1
  # if v == max(x,y), should work normally
  v <- max(x,y)
  returned <- eadelta1_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k <- eadelta1(n = 2*k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k_plus_1 <- eadelta1(n = 2*k+1, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned[1], S_2k_plus_1)
  expect_equal(returned[2], S_2k)
  # if v > max(x,y), should work normally
  v <- max(x,y)+1
  returned <- eadelta1_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k <- eadelta1(n = 2*k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k_plus_1 <- eadelta1(n = 2*k+1, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned[1], S_2k_plus_1)
  expect_equal(returned[2], S_2k)
  # if min == min(x,y), should return c(NaN, NaN) (should be using eadelta2_intervals in this case)
  returned <- eadelta1_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min(x,y), v = v)
  expect_equal(returned, c(NaN, NaN))
  # if min(x,y) < min (already outside layer), should return 0
  # expect 0 if x < min
  returned <- eadelta1_intervals(k = k, x = min-1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
  # expect 0 if y < min
  returned <- eadelta1_intervals(k = k, x = x, y = min-1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
  # if max(x,y) > v (already outside layer), should return 0
  # expect 0 if x > v
  returned <- eadelta1_intervals(k = k, x = v+1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
  # expect 0 if y > v
  returned <- eadelta1_intervals(k = k, x = x, y = v+1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
})

test_that("Check eadelta2_intervals returns expected results", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  min <- min(x,y)
  # if v == max(x,y), should work normally
  v <- max(x,y)
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- eadelta2_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k <- eadelta2(n = 2*k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k_plus_1 <- eadelta2(n = 2*k+1, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned[1], S_2k_plus_1)
  expect_equal(returned[2], S_2k)
  # if v > max(x,y), should work normally
  v <- max(x,y)+1
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned <- eadelta2_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k <- eadelta2(n = 2*k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  S_2k_plus_1 <- eadelta2(n = 2*k+1, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_is(returned, "numeric")
  expect_equal(length(returned), 2)
  expect_false(any(is.na(returned), is.nan(returned), is.infinite(returned)))
  expect_equal(returned[1], S_2k_plus_1)
  expect_equal(returned[2], S_2k)
  # if min < min(x,y), should return an NaN (should be using eadelta1 in this case)
  returned <- eadelta2_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min(x,y)-1, v = v)
  expect_equal(returned, c(NaN, NaN))
  # if min(x,y) < min (already outside layer), should return 0
  # expect 0 if x < min
  returned <- eadelta2_intervals(k = k, x = min-1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
  # expect 0 if y < min
  returned <- eadelta2_intervals(k = k, x = x, y = min-1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
  # if max(x,y) > v (already outside layer), should return 0
  # expect 0 if x > v
  returned <- eadelta2_intervals(k = k, x = v+1, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
  # expect 0 if y > v
  returned <- eadelta2_intervals(k = k, x = x, y = v+1, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned, c(0,0))
})

test_that("Check eadelta_intervals returns expected results", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  # if min(x,y) > min, eadelta should equal eadelta1_intervals
  min <- min(x,y)-1
  k <- sample(x = 0:10, size = 1)
  # if v == max(x,y)
  v <- max(x,y)
  returned_delta <- eadelta_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta1 <- eadelta1_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta1)
  # if v > max(x,y)
  v <- max(x,y)+1
  returned_delta <- eadelta_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta1 <- eadelta1_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta1)
  # if min(x,y) == min, eadelta should equal eadelta2_intervals
  min <- min(x,y)
  # if v == max(x,y)
  v <- max(x,y)
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned_delta <- eadelta_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta2 <- eadelta2_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta2)
  # if v > max(x,y)
  v <- max(x,y)+1
  D <- abs(v-min)
  k_min <- sqrt(times[2]-times[1] + D^2) / (2*D)
  k <- max(sample(x = 1:10, size = 1), ceiling(k_min))
  returned_delta <- eadelta_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  returned_delta2 <- eadelta2_intervals(k = k, x = x, y = y, s = times[1], t = times[2], min = min, v = v)
  expect_equal(returned_delta, returned_delta2)
})
