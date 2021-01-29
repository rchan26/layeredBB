context("Checking Inverse Gaussian sampler")

test_that("Inverse Gaussian sampler returns NaN values for mu = 0 or lambda = 0", {
  mu <- runif(1, 0.01, 10)
  lambda <- runif(1, 0.01, 10)
  expect_equal(inv_gauss_sampler(mu = 0, lambda = lambda), NaN)
  expect_equal(inv_gauss_sampler(mu = mu, lambda = 0), NaN)
  expect_equal(inv_gauss_sampler(mu = 0, lambda = 0), NaN)
})

test_that("Inverse Gaussian sampler returns a numeric value", {
  expect_is(inv_gauss_sampler(mu = runif(1, 0, 10), lambda = runif(1, 0, 10)), "numeric")
  # even in cases where returns NaN, should be numeric
  expect_is(inv_gauss_sampler(mu = runif(1, 0, 10), lambda = 0), "numeric")
  expect_is(inv_gauss_sampler(mu = 0, lambda = runif(1, 0, 10)), "numeric")
  expect_is(inv_gauss_sampler(mu = 0, lambda = 0), "numeric")
})
