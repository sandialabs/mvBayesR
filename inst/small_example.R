## simulate data (Friedman function with first variable as functional)
\donttest{f = function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - .5)^2 + 10 * x[, 4] + 5 *
    x[, 5]
}
sigma = 1 # noise sd
n = 500 # number of observations
n_grid = 50 # size of functional variable grid
t_grid = seq(0, 1, length.out = n_grid) # functional grid
X = matrix(runif(n * 9), n, 9) # 9 non-functional variables, only first 4 matter
X_expanded = cbind(rep(t_grid, each = n), kronecker(rep(1, n_grid), X)) # to get y
Y = matrix(f(X_expanded), nrow = n) + rnorm(n * n_grid, 0, sigma)

## fit Multivariate BASS
fit = mvBayes(
  BASS::bass, X, Y, nBasis = 3, # mvBayes parameters
  samplesExtract = function(bm) list(
    s2 = bm$s2,
    nbasis = bm$nbasis,
    beta.prec = bm$beta.prec
  ),
  idxSamplesArg = 'mcmc.use'
)
}
