# Load necessary libraries
library(mvBayes)
library(BayesPPR)

# Generate Data
f = function(x){
  out = 10.0 * sin(pi * tt * x[1]) +
    20.0 * (x[2] - 0.5)^2 +
    10 * x[3] +
    5.0 * x[4]
  return(out)
}

tt = seq(0, 1, length.out = 50)  # functional variable grid
n = 500  # sample size
p = 9  # number of predictors other (only 4 are used)
X = matrix(runif(n * p), n, p)  # training inputs
Y = t(apply(X, 1, f)) + matrix(rnorm(n * length(tt), sd = 0.1), n, length(tt))  # training response

ntest = 1000
Xtest = matrix(runif(ntest * p), ntest, p)
Ytest = t(apply(Xtest, 1, f)) + matrix(rnorm(ntest * length(tt), sd = 0.1), ntest, length(tt))

# Fit a multivariate BayesPPR model
mod = mvBayes(
  bppr,
  X,
  Y,
  idxSamplesArg = 'idx_use',  # 'idx_use' is bppr's equivalent of idxSamples
  # optionally extract posterior samples of key parameters
  samplesExtract = function(bppr_out) list(
    n_ridge = bppr_out$n_ridge,
    residSD = bppr_out$sd_resid,
    var_coefs = bppr_out$var_coefs
  )
)
plot(mod$basisInfo, idxMV = tt, xlabel = "tt")  # Plot PCA decomposition
traceplot(mod)
plot(mod, idxMV = tt, xlabel = "tt") # Evaluate training data fit
plot(mod, Xtest = Xtest, Ytest = Ytest, idxMV = tt, xlabel = "tt")  # Evaluate test data fit
modSensitivity = sobol(mod, nMC = 2^12)
plot(modSensitivity, idxMV = tt, xlabel = "tt")

# All posterior predictive samples
Ytest_postSamples = predict(mod, Xtest)
# Posterior predictive mean
Ytest_postMean = apply(Ytest_postSamples, 2, mean)
# single posterior predictive sample (from MCMC iteration #429)
Ytest_postSample429 = predict(mod, Xtest, idxSamples = 429)

# Use mvBayesElastic for a joint elastic functional PCA basis
mod = mvBayesElastic(
  bppr,
  X,
  Y,
  nBasis = 3,
  idxSamplesArg = 'idx_use',
  residSDExtract = function(bppr_out) sqrt(bppr_out$samples$s2)  # optionally extract posterior samples of residual standard deviation
)
plot(mod$basisInfo)  # Plot PCA decomposition
traceplot(mod)
plot(mod)  # Evaluate training data fit
plot(mod, Xtest = Xtest, Ytest = Ytest)  # Evaluate test data fit
modSensitivity = sobol(mod)
plot(modSensitivity)

# All posterior predictive samples
Ytest_postSamples = mod$predict(Xtest)
