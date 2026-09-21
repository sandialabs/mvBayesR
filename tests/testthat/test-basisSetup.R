# Tests for R/basisSetup.R (and its S3 methods)

# A small deterministic response matrix used across tests.
make_Y <- function(n = 40, p = 8, seed = 1) {
  set.seed(seed)
  t_grid <- seq(0, 1, length.out = p)
  scores <- matrix(rnorm(n * 3), n, 3)
  loadings <- rbind(sin(2 * pi * t_grid),
                    cos(2 * pi * t_grid),
                    t_grid)
  scores %*% loadings + matrix(rnorm(n * p, sd = 0.01), n, p)
}

test_that("basisSetup returns an object of the expected class and structure", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 3)

  expect_s3_class(bs, "basisSetup")
  expect_equal(bs$nMV, ncol(Y))
  expect_equal(bs$nBasis, 3)
  expect_equal(dim(bs$basis), c(3, ncol(Y)))
  expect_equal(dim(bs$coefs), c(nrow(Y), 3))
})

test_that("PCA basis vectors are orthonormal", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 3)
  gram <- bs$basis %*% t(bs$basis)
  expect_equal(gram, diag(3), tolerance = 1e-8)
})

test_that("propVarExplained selects enough components to reach the threshold", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", propVarExplained = 0.99)
  expect_gte(bs$propVarExplained, 0.99)
  # The three-component signal should be recoverable in a few bases
  expect_lte(bs$nBasis, ncol(Y))
})

test_that("nBasis larger than the number of available bases warns and is clamped", {
  Y <- make_Y(p = 4)
  expect_warning(
    bs <- basisSetup(Y, basisType = "pca", nBasis = 100),
    "larger than"
  )
  expect_lte(bs$nBasis, ncol(Y))
})

test_that("an unsupported basisType raises an error", {
  Y <- make_Y()
  expect_error(basisSetup(Y, basisType = "not-a-basis"), "un-supported")
})

test_that("getCoefs recovers the stored coefficients when Ytest is NULL", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 3)
  expect_equal(getCoefs(bs), bs$coefs)
})

test_that("getCoefs on the training data reproduces the stored coefficients", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 3)
  # Projecting the original Y back should recover the first nBasis columns
  coefsTest <- getCoefs(bs, Ytest = Y)[, 1:bs$nBasis, drop = FALSE]
  expect_equal(coefsTest, bs$coefs, tolerance = 1e-8)
})

test_that("getYtrunc reconstructs Y closely when using enough components", {
  Y <- make_Y()
  # Use as many components as columns so reconstruction is (near) exact
  bs <- basisSetup(Y, basisType = "pca", nBasis = ncol(Y), propVarExplained = 1)
  Yhat <- getYtrunc(bs)
  expect_equal(dim(Yhat), dim(Y))
  expect_equal(Yhat, Y, tolerance = 1e-6)
})

test_that("truncError equals Y minus the reconstruction", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 2)
  expect_equal(bs$truncError, Y - getYtrunc(bs), tolerance = 1e-10)
})

test_that("centering stores the column means", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 3, center = TRUE)
  expect_equal(bs$Ycenter, apply(Y, 2, mean), tolerance = 1e-10)
})

test_that("scaling stores the column standard deviations", {
  Y <- make_Y()
  bs <- basisSetup(Y, basisType = "pca", nBasis = 3, scale = TRUE)
  expect_equal(bs$Yscale, apply(Y, 2, sd), tolerance = 1e-10)
})

test_that("legendre basis requires nBasis and enforces the minimum", {
  Y <- make_Y()
  expect_error(basisSetup(Y, basisType = "legendre"), "nBasis")
  expect_error(basisSetup(Y, basisType = "legendre", nBasis = 2), "nBasis >= 3")
})

test_that("custom basis must match the number of response columns", {
  Y <- make_Y(p = 8)
  goodBasis <- diag(8)[1:3, ]
  badBasis <- diag(5)[1:3, ]

  expect_error(basisSetup(Y, basisType = "custom"), "customBasis")
  expect_error(
    basisSetup(Y, basisType = "custom", customBasis = badBasis),
    "must match"
  )
  bs <- basisSetup(Y, basisType = "custom", customBasis = goodBasis)
  expect_s3_class(bs, "basisSetup")
})
