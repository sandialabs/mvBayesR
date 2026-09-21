# End-to-end test of the mvBayes fit + predict pipeline using a lightweight,
# dependency-free mock "Bayesian" model. The mock returns a small fixed set of
# posterior draws so the whole path (basisSetup -> fit -> predict) is exercised
# without MCMC or the BASS dependency.

# --- Mock bayesModel ---------------------------------------------------------
# Fits an ordinary least-squares line and fabricates `nDraws` posterior samples
# by jittering the OLS coefficients. Returns an object with a `predict` method
# whose output is a length nSamples*ntest vector, as mvBayes expects.
mockBayes <- function(X, y, nDraws = 5, ...) {
  X <- as.matrix(X)
  Xd <- cbind(1, X)
  beta <- qr.solve(Xd, y)
  set.seed(42)
  betaDraws <- t(replicate(nDraws, beta + rnorm(length(beta), sd = 1e-6)))
  structure(list(betaDraws = betaDraws, p = ncol(X)), class = "mockBayes")
}

predict.mockBayes <- function(object, Xtest, idxSamples = "default", ...) {
  Xtest <- as.matrix(Xtest)
  Xd <- cbind(1, Xtest)
  draws <- object$betaDraws
  if (!identical(idxSamples, "default")) {
    draws <- draws[idxSamples, , drop = FALSE]
  }
  # mvBayes expects the per-component predictions as an (nSamples x ntest)
  # matrix (rows = posterior samples, cols = test points): it uses
  # length()/ntest to infer nSamples and then indexes rows as postCoefs_k[it, ].
  draws %*% t(Xd)
}

make_data <- function(n = 30, p = 6, seed = 7, noise = 0.01) {
  set.seed(seed)
  X <- matrix(runif(n * 2), n, 2)
  t_grid <- seq(0, 1, length.out = p)
  Y <- outer(X[, 1], sin(2 * pi * t_grid)) +
    outer(X[, 2], t_grid) +
    matrix(rnorm(n * p, sd = noise), n, p)
  list(X = X, Y = Y)
}

test_that("mvBayes fits and returns a well-formed object", {
  registerS3method("predict", "mockBayes", predict.mockBayes)
  d <- make_data()

  fit <- suppressMessages(
    mvBayes(mockBayes, d$X, d$Y, nBasis = 2, nCores = 1)
  )

  expect_s3_class(fit, "mvBayes")
  expect_equal(fit$basisInfo$nBasis, 2)
  expect_length(fit$bmList, 2)
  expect_true(fit$nSamples >= 1)
})

test_that("predict.mvBayes returns an array of the expected shape", {
  registerS3method("predict", "mockBayes", predict.mockBayes)
  d <- make_data()

  fit <- suppressMessages(
    mvBayes(mockBayes, d$X, d$Y, nBasis = 2, nCores = 1)
  )
  Xtest <- d$X[1:5, , drop = FALSE]
  pred <- suppressMessages(predict(fit, Xtest))

  # array of dimension c(nSamples, ntest, nMV)
  expect_equal(dim(pred)[1], fit$nSamples)
  expect_equal(dim(pred)[2], nrow(Xtest))
  expect_equal(dim(pred)[3], ncol(d$Y))
})

test_that("predict.mvBayes recovers the training response closely", {
  registerS3method("predict", "mockBayes", predict.mockBayes)
  # Noise-free signal: Y is exactly linear in X per column, so a full-rank
  # basis + per-component OLS should reconstruct it essentially exactly.
  d <- make_data(noise = 0)

  # Use enough bases for a near-exact reconstruction
  fit <- suppressMessages(
    mvBayes(mockBayes, d$X, d$Y, nBasis = ncol(d$Y),
            propVarExplained = 1, nCores = 1)
  )
  pred <- suppressMessages(predict(fit, d$X))
  # pred is c(nSamples, ntest, nMV); average over the posterior-sample axis
  postMean <- apply(pred, c(2, 3), mean)   # ntest x nMV
  expect_equal(postMean, d$Y, tolerance = 1e-4)
})

test_that("returnPostCoefs yields both Ypost and postCoefs", {
  registerS3method("predict", "mockBayes", predict.mockBayes)
  d <- make_data()

  fit <- suppressMessages(
    mvBayes(mockBayes, d$X, d$Y, nBasis = 2, nCores = 1)
  )
  out <- suppressMessages(
    predict(fit, d$X[1:4, , drop = FALSE], returnPostCoefs = TRUE)
  )
  expect_true(all(c("Ypost", "postCoefs") %in% names(out)))
  expect_equal(dim(out$postCoefs)[3], fit$basisInfo$nBasis)
})
