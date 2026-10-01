# Load necessary libraries
library(mvBayes)
library(BASS)

# Generate Data
f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}

sigma<-1
nfunc = 50
tt = seq(0, 1, length.out = nfunc)  # functional variable grid
n = 500  # sample size
p = 9  # number of predictors other (only 4 are used)
X<-matrix(runif(n*p),n,p) # 9 non-functional variables, only first 4 matter
x<-cbind(rep(tt,each=n),kronecker(rep(1,nfunc),X)) # to get y
Y<-matrix(f(x),nrow=n)+rnorm(n*nfunc,0,sigma)

ntest = 1000
Xtest = matrix(runif(ntest * p), ntest, p)
x<-cbind(rep(tt,each=ntest),kronecker(rep(1,nfunc),Xtest)) # to get y
Ytest = matrix(f(x),nrow=ntest)+rnorm(ntest*nfunc,0,sigma)

# Fit a multivariate BASS model
mod = mvBayes(
  bass,
  X,
  Y,
  nBasis=3,
  samplesExtract=function(fit){
    list(
      residSD = sqrt(fit$s2),
      log.post.cold = fit$log.post.cold[(fit$nburn + 1):fit$nmcmc],
      lam = fit$lam,
      beta.prec = fit$beta.prec
    )
  },
  idxSamplesArg='mcmc.use'
)
plot(mod)
plot(mod$basisInfo, idxMV = tt, xlabel = "tt")  # Plot PCA decomposition
traceplot(mod)
plot(mod, idxMV = tt, xlabel = "tt") # Evaluate training data fit
plot(mod, Xtest = Xtest, Ytest = Ytest, idxMV = tt, xlabel = "tt")  # Evaluate test data fit
modSensitivity = mvSobol(mod, nMC = 2^9)
plot(modSensitivity, idxMV = tt, xlabel = "tt")

# All posterior predictive samples
Ytest_postSamples = predict(mod, Xtest)
# Posterior predictive mean
Ytest_postMean = apply(Ytest_postSamples, 2, mean)
# single posterior predictive sample (from MCMC iteration #429)
Ytest_postSample429 = predict(mod, Xtest, idxSamples = 429)

