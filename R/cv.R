#' @title Cross-Validation (CV) of a Multivariate Bayesian Regression Model
#'
#' @description Wrapper to fit and perform cross-validation for a multivariate Bayesian regression model, using the mvBayes function.
#' @param bayesModel A Bayesian regression model-fitting function, with first argument taking an nxp input matrix or data.frame, and second argument taking an n-vector of numeric responses.
#' @param X A matrix of predictors of dimension nxp, where n is the total number of examples (including training and test sets) and p is the number of inputs (features).
#' @param Y A response matrix of dimension nxq, where q is the number of multivariate/functional responses.
#' @param kFolds Number of test sets to partition data. If kFolds=NULL (default), test sets are instead formed by random sample, via the nReps, nTrain, and nTest arguments.
#' @param nRep Number of repetitions of CV process. Only used when kFolds=NULL (default).
#' @param nTrain Number of examples to use in the training set. Only used when kFolds=NULL (default). If nTrain=NULL, nTrain is set to n - nTest; unless nTest is also NULL, in which case nTrain is set to ceiling(n/2).
#' @param nTest Number of examples to use in the test set. Only used when kFolds=NULL (default). If nTest=NULL, nTest is set to n - nTrain.
#' @param seed Randomization seed, for replication of the train/test split. The seed is un-initialized immediately after assigning the train/test split. If NULL, no seed is set.
#' @param coverageTarget level of coverage desired (default: 0.95)
#' @param idxSamples which samples to use in CV (default: "all")
#' @param uqTruncMethod method to use for UQ truncation (c("gaussian", "empirical"))
#' @param warpData `time_warping` object from `fdasrvf`. If supplied, \link{mvBayesElastic} is used instead of \link{mvBayes}, and each fold is aligned using the corresponding columns of `warpData`.
#' @param ... Additional arguments to mvBayes, including arguments to bayesModel.
#' @details First separates the data into randomly chosen test and training sets (user-specified train/test splits and k-fold cv are forthcoming), then fits mvBayes(bayesModel, X, Y, ...) to the training set and evaluates predictive performance on the test set. Repeats this process nRep times.
#' @return An object of class "mvBayesCV", which is a list containing out-of-sample performance metrics for each replication, including rmse, rSquared, coverage, intervalWidth, intervalScore, crps, fitting and prediction times, and the function call.
#' @seealso \link{mvBayes}, \link{predict.mvBayes} for prediction
#' @export
mvCV = function(bayesModel,
                X,
                Y,
                kFolds = NULL,
                nRep = 1,
                nTrain = NULL,
                nTest = NULL,
                seed = NULL,
                coverageTarget = 0.95,
                idxSamples = "all",
                uqTruncMethod = c("gaussian", "empirical"),
                warpData = NULL,
                ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # setup
  n = nrow(X)
  
  alpha = 1 - coverageTarget
  
  uqTruncMethod = match.arg(uqTruncMethod)
  
  if (is.null(kFolds)) {
    if (is.null(nTest)) {
      if (is.null(nTrain)) {
        nTest = floor(n / 2) # half in test set
        nTrain = n - nTest
      } else if (nTrain >= n) {
        stop('Must have nTrain < nrow(X)')
      } else{
        nTest = n - nTrain
      }
    } else{
      if (nTest >= n) {
        stop('Must have nTest < nrow(X)')
      } else if (is.null(nTrain)) {
        nTrain = n - nTest
      } else if (nTrain + nTest > n) {
        stop('Must have nTrain + nTest <= n')
      }
    }
    
    # Get fold indices
    idxTest = lapply(1:nRep, function(r)
      sample(n, size = nTest)) # different test set for every rep
    nTest = rep(nTest, nRep)
    # nTrain is honored as given: it may be smaller than n - nTest, in which
    # case the training set is a subsample of the non-test observations
    nTrain = rep(nTrain, nRep)
    idxRemaining = lapply(idxTest, function(idx)
      setdiff(1:n, idx)) # remaining indices after test set is determined
    idxTrain = lapply(1:nRep, function(r)
      sample(idxRemaining[[r]], size = nTrain[r])) # training set for each rep
    rm(idxRemaining)
  } else {
    nRep = kFolds
    idxTest = split(1:n, rep_len(1:kFolds, n))
    nTest = sapply(idxTest, length)
    nTrain = n - nTest
    idxTrain = lapply(idxTest, function(idx)
      setdiff(1:n, idx)) # remaining indices after test set is determined
  }
  
  if (!is.null(seed)) {
    set.seed(NULL) # re-set as if no seed had been set
  }
  
  # Run cv
  rmse = rSquared = coverage = intervalWidth = intervalScore = crps = fitTime = predictTime = numeric(nRep)
  for (r in 1:nRep) {
    # Set up train/test split
    Xtrain = X[idxTrain[[r]], ,drop = FALSE]
    Ytrain = Y[idxTrain[[r]], ,drop = FALSE]
    
    Xtest = X[idxTest[[r]], ,drop = FALSE]
    Ytest = Y[idxTest[[r]], ,drop = FALSE]
    
    # Fit models
    useElastic = !is.null(warpData)
    startFit = Sys.time()
    if (useElastic) {
      fit = mvBayesElastic(
        bayesModel,
        Xtrain,
        Ytrain,
        warpData = warpData,
        idx = idxTrain[[r]],
        ...
      )
    } else {
      fit = mvBayes(bayesModel, Xtrain, Ytrain, ...)
    }
    fitTime[r] = as.numeric(Sys.time() - startFit, units = "secs")
    
    # Calculate rmse of posterior mean
    start_pred = Sys.time()
    preds = predict(fit, Xtest)
    if (idxSamples[1] != "all") {
      if (identical(idxSamples, "final")) {
        idxSamplesUse = dim(preds)[1]
      } else if (is.numeric(idxSamples)) {
        idxSamplesUse = idxSamples
      } else {
        stop("'idxSamples' must be 'all', 'final', or numeric.")
      }
      # drop=FALSE keeps the nSamples dimension for a single index
      preds = preds[idxSamplesUse, , , drop = FALSE]
    }
    predictTime[r] = as.numeric(Sys.time() - start_pred, units = "secs")
    
    Yhat = matrix(apply(preds, 2:3, median),
                  nrow = dim(preds)[2],
                  ncol = dim(preds)[3])
    if (useElastic) {
      basisType = fit$basisInfo$basisType
      C = fit$basisInfo$basisConstruct$C
      id = fit$basisInfo$basisConstruct$id
      srvf = fit$basisInfo$basisConstruct$srvf
      
      fnTest = warpData$fn[, idxTest[[r]], drop = FALSE]
      if (srvf) {
        m_new = sign(fnTest[id, ]) * sqrt(abs(fnTest[id, ]))
        qn = fdasrvf::f_to_srvf(fnTest, fit$basisInfo$basisConstruct$time)
        qn1 = rbind(qn, m_new)
      } else {
        qn1 = fnTest
      }
      
      gamTest = warpData$warping_functions[, idxTest[[r]], drop = FALSE]
      
      if (basisType == "jfpca") {
        time = seq(0, 1, length.out = ncol(Ytest))
        binsize <- mean(diff(time))
        vec = matrix(0, ncol(Ytest), nrow(Ytest))
        for (i in 1:nrow(Ytest)) {
          psi = sqrt(fdasrvf::gradient(gamTest[, i], binsize))
          vec[, i] <- fdasrvf::inv_exp_map(fit$basisInfo$basisConstruct$mu_psi, psi)
        }
        Ytest = t(rbind(qn1, C * vec))
      } else if (basisType == "jfpcah") {
        h = fdasrvf::gam_to_h(gamTest)
        Ytest = t(rbind(qn1, C * h))
      }
    }
    rmse[r] = sqrt(mean((Ytest - Yhat)^2))
    rSquared[r] = 1 - mean((Ytest - Yhat)^2) / mean((t(Ytest) - fit$basisInfo$Ycenter)^2)
    
    # Get truncation error for UQ
    if (uqTruncMethod == "gaussian") {
      truncErrorVar = cov(fit$basisInfo$truncError)
      truncError = array(
        MASS::mvrnorm(
          prod(dim(preds)[1:2]),
          rep(0, dim(preds)[3]),
          truncErrorVar
        ),
        dim = dim(preds)
      )
    } else if (uqTruncMethod == "empirical") {
      idxResample = sample(nTrain[r], size = prod(dim(preds)[1:2]), replace = TRUE)
      truncError = aperm(
        array(
          t(fit$basisInfo$truncError[idxResample, ]),
          dim = dim(preds)[c(3, 1, 2)]
        ),
        c(2, 3, 1)
      )
    }
    preds = preds + truncError
    rm(truncError)
    
    # Get regression error for UQ
    coefsResidError = array(dim = c(dim(preds)[1:2], fit$basisInfo$nBasis))
    for (k in 1:fit$basisInfo$nBasis) {
      if (class(fit$bmList[[k]])[1] %in% c("gbass", "tbass", "qbass", "nwbass")) {
        w <- fit$bmList[[k]]$w
        beta <- fit$bmList[[k]]$beta
        v <- fit$bmList[[k]]$v
        coefsResidMean <- sqrt(w) * beta * v
        coefsResidSD <- sqrt(w * v)
      } else {
        coefsResidMean <- 0
        coefsResidSD = fit$bmList[[k]]$samples$residSD
      }
      coefsResidError[, , k] = rnorm(
        prod(dim(preds)[1:2]),
        mean = coefsResidMean,
        sd = coefsResidSD
      )
    }
    residError = array(dim = dim(preds))
    basisScaledT = t(t(fit$basisInfo$basis) * fit$basisInfo$Yscale)
    for (idxMCMC in 1:dim(preds)[1]) {
      coefsResidErrorMC = matrix(
        coefsResidError[idxMCMC, , ],
        nrow = dim(preds)[2],
        ncol = fit$basisInfo$nBasis
      )
      residError[idxMCMC, , ] = coefsResidErrorMC %*% basisScaledT
    }
    rm(coefsResidError)
    preds = preds + residError
    rm(residError)
    
    # Calculate CRPS
    # CRPS is computed pointwise for each response dimension and then averaged
    # over response dimensions and test observations.
    nSamples = dim(preds)[1]
    nObs = dim(preds)[2]
    nResp = dim(preds)[3]
    
    weights = 2 * (1:nSamples) - nSamples - 1
    yMat = matrix(0, nrow = nSamples, ncol = nObs)
    crps_sum = 0
    
    for (j in 1:nResp) {
      pred_j = preds[, , j, drop = FALSE][, , 1]   # nSamples x nObs
      y_j = Ytest[, j]                             # length nObs
      
      # Faster than constructing matrix(y_j, ..., byrow = TRUE) each iteration
      yMat[] = y_j
      term1 = colMeans(abs(pred_j - yMat))
      
      # Compute sorted columns with an explicit loop to avoid apply() overhead
      pred_j_sort = matrix(NA_real_, nrow = nSamples, ncol = nObs)
      for (i in 1:nObs) {
        pred_j_sort[, i] = sort.int(pred_j[, i], method = "auto")
      }
      
      # Closed-form for 0.5 * E|X - X'|
      term2 = as.numeric(crossprod(weights, pred_j_sort)) / (nSamples^2)
      
      crps_sum = crps_sum + sum(term1 - term2)
    }
    
    crps[r] = crps_sum / (nObs * nResp)
    
    # Calculate distance from posterior mean
    distBound = numeric(dim(preds)[2])
    for (idx in 1:dim(preds)[2]) {
      predsIdx = matrix(preds[, idx, ], nrow = nSamples, ncol = dim(preds)[3])
      distSamples = sqrt(rowMeans((predsIdx - matrix(
        Yhat[idx, ], nSamples, dim(preds)[3], byrow = TRUE
      ))^2))
      distBound[idx] = quantile(distSamples, coverageTarget)
    }
    distTest = sqrt(apply((Ytest - Yhat)^2, 1, mean))
    
    # Calculate UQ metrics
    distRatio = distTest / distBound
    coverage[r] = mean(distRatio <= 1)
    intervalWidth[r] = exp(mean(log(distBound)))
    intervalScore[r] = intervalWidth[r] * exp(mean(log(distRatio) * (distRatio > 1)) /
                                                alpha)
  }
  
  out = list(
    rmse = rmse,
    rSquared = rSquared,
    coverageTarget = coverageTarget,
    coverage = coverage,
    intervalWidth = intervalWidth,
    intervalScore = intervalScore,
    crps = crps,
    fitTime = fitTime,
    predictTime = predictTime,
    call = match.call()
  )
  
  return(structure(out, class = 'mvBayesCV'))
}
