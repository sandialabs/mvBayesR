## simulate data (Friedman function with first variable as functional)
f = function(x) {
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
  nburn = 8500, # bass parameter
  samplesExtract = function(bm) list(
    s2 = bm$s2,
    nbasis = bm$nbasis,
    beta.prec = bm$beta.prec
  ),
  idxSamplesArg = 'mcmc.use'
)
plot(fit$basisInfo, title = 'Principal Component Decomposition of Y')
traceplot(
  fit,
  modelParams = c('s2', 'nbasis', 'beta.prec'),
  labels = c(
    'sigma^2',
    'Number of BASS Basis Vectors',
    'Coefficient Precision'
  ),
  title = 'Traceplots of BASS parameters'
)
plot(fit, title = 'Multivariate BASS Model Fit')


## prediction
ntest = 1000
Xtest = matrix(runif(ntest * 9), ntest, 9)
Xtest_expanded = cbind(rep(t_grid, each = ntest), kronecker(rep(1, n_grid), Xtest))
Ytest = matrix(f(Xtest_expanded), nrow = ntest) + rnorm(ntest * n_grid, 0, sigma)
pred = predict(fit, Xtest) # posterior predictive samples (each is a curve)
plot(fit,
     Xtest = Xtest,
     Ytest = Ytest,
     title = 'Multivariate BASS Out-of-sample Predictions')


## cross-validation
cv = mvBayesCV(
  BASS::bass,
  X,
  Y,
  nTest = 100,
  nRep = 2,
  # mvBayesCV parameters
  nBasis = 3,
  # mvBayes parameters
  nburn = 8500 # bass parameters
)

## sensitivity
sens = mvSobol(fit, nMC = 2^12, totalSobol = FALSE) # use BASS functionality to compute mvSobol

## Write a wrapper for a non-conforming function to conform to mvBayes requirements
linmod = function(X, y, ...)
  structure(lm(y ~ X, ...), class = 'linmod')
predict.linmod = function(object, newdata)
  predict(structure(object, class = 'lm'), as.data.frame(newdata))
fit = mvBayes(
  linmod,
  X,
  Y,
  nBasis = 5,
  samplesExtract = function(bm) list(
    residSD = sd(bm$residuals)
  )
)
yhat = predict(fit, X)
