# Tests for R/legendre.R

test_that("basisLegendre returns a matrix of the expected shape", {
  M <- 20
  N <- 4
  B <- basisLegendre(M, N)
  # (N + 1) polynomial orders evaluated at M grid points
  expect_equal(dim(B), c(N + 1, M))
})

test_that("basisLegendre reproduces the low-order Legendre polynomials", {
  M <- 101
  N <- 2
  B <- basisLegendre(M, N)
  x <- seq(-1, 1, length.out = M)

  # P0(x) = 1
  expect_equal(as.numeric(B[1, ]), rep(1, M), tolerance = 1e-8)
  # P1(x) = x
  expect_equal(as.numeric(B[2, ]), x, tolerance = 1e-8)
  # P2(x) = (3x^2 - 1) / 2
  expect_equal(as.numeric(B[3, ]), (3 * x^2 - 1) / 2, tolerance = 1e-8)
})
