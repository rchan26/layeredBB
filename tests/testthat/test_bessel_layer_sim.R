context("Checking Bessel layer simulation functions")

test_that("One-dimensional Bessel layer simulation returns errors if argumments are outside constraints", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  # should have error if s >= t
  # expect error if s > t
  expect_error(bessel_layer_simulation(x = x,
                                       y = y,
                                       s = times[2],
                                       t = times[1],
                                       mult = mult),
               "layeredBB::bessel_layer_simulation: s >= t. Must have s < t")
  # expect error if s == t
  expect_error(bessel_layer_simulation(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[1],
                                       mult = mult),
               "layeredBB::bessel_layer_simulation: s >= t. Must have s < t")
  expect_error(bessel_layer_simulation(x = x,
                                       y = y,
                                       s = times[2],
                                       t = times[2],
                                       mult = mult),
               "layeredBB::bessel_layer_simulation: s >= t. Must have s < t")
  # should have error if mult <= 0
  # expect error if mult < 0
  expect_error(bessel_layer_simulation(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[2],
                                       mult = -mult),
               "layeredBB::bessel_layer_simulation: mult <= 0. Must have mult > 0")
  # expect error if mult = 0
  expect_error(bessel_layer_simulation(x = x,
                                       y = y,
                                       s = times[1],
                                       t = times[2],
                                       mult = 0),
               "layeredBB::bessel_layer_simulation: mult <= 0. Must have mult > 0")
})

test_that("One-dimensional Bessel layer simualtion returns expected results", {
  x <- rnorm(1, 0, 10)
  y <- rnorm(1, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  returned <- bessel_layer_simulation(x = x,
                                      y = y,
                                      s = times[1],
                                      t = times[2],
                                      mult = mult)
  # check the returned class and length
  expect_is(returned, "list")
  expect_equal(length(returned), 4)
  expect_equal(names(returned), c("L", "l", "u", "U"))
  # check that it satisfies the constraints of what a Bessel layer should be
  expect_true(all(returned$L < min(x,y),
                  returned$l <= min(x,y),
                  returned$u >= max(x,y),
                  returned$U > max(x,y)))
  expect_true(all(returned$L < returned$l,
                  returned$l < returned$u,
                  returned$u < returned$U))
  # check no NA, NaN, Inf
  returned_ <- unlist(returned)
  expect_false(any(is.na(returned_), is.nan(returned_), is.infinite(returned_)))
})

test_that("Multi-dimensional Bessel layer simulation returns errors if arguments are outside constraints", {
  dim <- sample(x = 1:10, size = 1)
  x <- rnorm(dim, 0, 10)
  y <- rnorm(dim, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  # should have error if x and y aren't vectors of length dim
  # expect error if x is not vector of length dim
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = rnorm(dim-1),
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: length of x is not equal to dim")
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = rnorm(dim+1),
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: length of x is not equal to dim")
  # expect error if y is not vector of length dim
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = rnorm(dim-1),
                                             s = times[1],
                                             t = times[2],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: length of y is not equal to dim")
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = rnorm(dim+1),
                                             s = times[1],
                                             t = times[2],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: length of y is not equal to dim")
  # should have error if s >= t
  # expect error if s > t
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[2],
                                             t = times[1],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: s >= t. Must have s < t")
  # expect error if s == t
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[1],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: s >= t. Must have s < t")
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[2],
                                             t = times[2],
                                             mult = mult),
               "layeredBB::multi_bessel_layer_simulation: s >= t. Must have s < t")
  # should have error if mult <= 0
  # expect error if mult < 0
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             mult = -mult),
               "layeredBB::multi_bessel_layer_simulation: mult <= 0. Must have mult > 0")
  # expect error if mult == 0
  expect_error(multi_bessel_layer_simulation(dim = dim,
                                             x = x,
                                             y = y,
                                             s = times[1],
                                             t = times[2],
                                             mult = 0),
               "layeredBB::multi_bessel_layer_simulation: mult <= 0. Must have mult > 0")

})

test_that("Multi-dimensional Bessel layer simulation returns expected results", {
  dim <- sample(x = 1:10, size = 1)
  x <- rnorm(dim, 0, 10)
  y <- rnorm(dim, 0, 10)
  times <- sort(runif(2, 0, 10))
  mult <- runif(1, 0.01, 10)
  returned <- multi_bessel_layer_simulation(dim = dim,
                                            x = x,
                                            y = y,
                                            s = times[1],
                                            t = times[2],
                                            mult = mult)
  # check the returned class and length
  expect_is(returned, "list")
  expect_equal(length(returned), dim)
  expect_true(all(sapply(returned, is.list)))
  # for each dimension, check that it satisfies the constraints of what a Bessel layer should be
  for (i in 1:dim) {
    item <- returned[[i]]
    expect_is(item, "list")
    expect_equal(length(item), 4)
    expect_equal(names(item), c("L", "l", "u", "U"))
    expect_true(all(item$L < min(x[i],y[i]),
                    item$l <= min(x[i],y[i]),
                    item$u >= max(x[i],y[i]),
                    item$U > max(x[i],y[i])))
    expect_true(all(item$L < item$l,
                    item$l < item$u,
                    item$u < item$U))
  }
  # check no NA, NaN, Inf
  returned_ <- unlist(returned)
  expect_false(any(is.na(returned_), is.nan(returned_), is.infinite(returned_)))
})
