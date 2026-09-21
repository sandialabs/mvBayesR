# Tests for utility helpers in R/utils.R

test_that("ndims returns the number of dimensions", {
  expect_equal(ndims(matrix(0, 2, 3)), 2)
  expect_equal(ndims(array(0, dim = c(2, 3, 4))), 3)
  # A plain vector has no dim attribute
  expect_equal(ndims(1:5), 0)
})

test_that("normalize_column produces a unit vector", {
  x <- c(3, 4)
  nx <- normalize_column(x)
  expect_equal(sqrt(sum(nx^2)), 1)
  expect_equal(nx, c(0.6, 0.8))
})

test_that("normalize_column leaves a zero-magnitude vector unchanged", {
  x <- c(0, 0, 0)
  expect_equal(normalize_column(x), x)
})

test_that("repmat tiles a matrix like MATLAB's repmat", {
  X <- matrix(1:4, nrow = 2)
  out <- repmat(X, 2, 3)
  expect_equal(dim(out), c(4, 6))
  # The top-left tile is the original matrix
  expect_equal(out[1:2, 1:2], X)
})

test_that("cumtrapz integrates a constant correctly", {
  x <- seq(0, 1, length.out = 11)
  y <- rep(2, length(x))
  z <- cumtrapz(x, y)
  # Integral of the constant 2 from 0 to 1 is 2
  expect_equal(as.numeric(z[length(z)]), 2, tolerance = 1e-8)
  # Starts at 0
  expect_equal(as.numeric(z[1]), 0)
})

test_that("cumtrapz integrates a linear function correctly", {
  x <- seq(0, 1, length.out = 101)
  y <- x
  z <- cumtrapz(x, y)
  # Integral of x from 0 to 1 is 0.5
  expect_equal(as.numeric(z[length(z)]), 0.5, tolerance = 1e-4)
})
