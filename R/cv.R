#' @title Cross-Validation (CV) of a Multivariate Bayesian Regression Model
#'
#' @description Wrapper to fit and perform cross-validation for a multivariate Bayesian regression model, using the mvBayes function.
#' @param bayesModel A Bayesian regression model-fitting function, with first argument taking an nxp input matrix or data.frame, and second argument taking an n-vector of numeric responses.
#' @param X A matrix of predictors of dimension nxp, where n is the total number of examples (including training and test sets) and p is the number of inputs (features).
#' @param Y A response matrix of dimension nxq, where q is the number of multivariate/functional responses.
#' @param nTrain Number of examples to use in the training set. If NULL, nTrain = n - nTest; unless nTest is also NULL, in which case nTrain = ceiling(n/2).
#' @param nTest Number of examples to use in the test set. If NULL, nTest = n - nTrain.
#' @param nRep Number of repetitions of CV process.
#' @param seed Randomization seed, for replication of the train/test split. The seed is un-initialized immediately after assigning the train/test split. If NULL, no seed is set.
#' @param coverageTarget level of coverage desired (default: 0.95)
#' @param idxSamples which samples to use in CV (default: "all")
#' @param uqTruncMethod method to use for UQ truncation (c("gaussian", "empirical"))
#' @param ... Additional arguments to mvBayes, including arguments to bayesModel.
#' @details First separates the data into randomly chosen test and training sets (user-specified train/test splits and k-fold cv are forthcoming), then fits mvBayes(bayesModel, X, Y, ...) to the training set and evaluates predictive performance on the test set. Repeats this process nRep times.
#' @return An object of class "mvBayesCV", which is a list containing the out-of-sample rmse for each replication, the fitting and prediction times, and the function call. Other prediction metrics, including coverage of prediction intervals, are forthcoming.
#' @keywords multivariate Bayesian regression modeling, functional data analysis
#' @seealso \link{mvBayes}, \link{predict.mvBayes} for prediction
#' @export
mvCV = function(bayesModel,
                     X,
                     Y,
                     nTrain = NULL,
                     nTest = NULL,
                     nRep = 1,
                     seed = NULL,
                     coverageTarget = 0.95,
                     idxSamples = "all",
                     uqTruncMethod = c("gaussian", "empirical"),
                     ...) {
  # setup
  n = nrow(X)
  p = ncol(X)
  nMV = ncol(Y)

  alpha = 1 - coverageTarget

  uqTruncMethod = uqTruncMethod[1]

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
  set.seed(seed)
  idxTest = sapply(1:nRep, function(r)
    sample(n, size = nTest)) # different test set for every rep
  idxRemaining = apply(idxTest, 2, function(idx)
    setdiff(1:n, idx)) # remaining indices after test set is determined
  idxTrain = apply(idxRemaining, 2, function(idx)
    sample(idx, size = nTrain)) # training set for each rep
  rm(idxRemaining)
  set.seed(NULL) # re-set as if no seed had been set


  # Run cv
  rmse = rSquared = coverage = intervalWidth = intervalScore = fit_time = pred_time = numeric(nRep)
  for (r in 1:nRep) {
    # Set up train/test split
    Xtrain = X[idxTrain[, r], ]
    Ytrain = Y[idxTrain[, r], ]

    Xtest = X[idxTest[, r], ]
    Ytest = Y[idxTest[, r], ]

    # Fit models
    if (methods::hasArg("warp_data")) {
      start_fit = Sys.time()
      fit = mvBayesElastic(bayesModel, Xtrain, Ytrain, idx = idxTrain[, r], ...)
      fit_time[r] = as.numeric(Sys.time() - start_fit, units = "secs")
    } else {
      start_fit = Sys.time()
      fit = mvBayes(bayesModel, Xtrain, Ytrain, ...)
      fit_time[r] = as.numeric(Sys.time() - start_fit, units = "secs")
    }

    # Calculate rmse of posterior mean
    start_pred = Sys.time()
    preds = predict(fit, Xtest)
    if (idxSamples[1] != "all") {
      preds = preds[idxSamples, , ]
    }
    pred_time[r] = as.numeric(Sys.time() - start_pred, units = "secs")

    Yhat = apply(preds, 2:3, median)
    if (methods::hasArg("warp_data")) {
      tmp = as.list(match.call())[-1]
      call.envir = parent.frame(1)
      if (fit$basisInfo$basisType == "jfpca")
      {
        C = fit$basisInfo$basisConstruct$C
        id = fit$basisInfo$basisConstruct$id
        srvf = fit$basisInfo$basisConstruct$srvf
        if (srvf) {
          m_new = sign(eval(tmp$warp_data, envir = call.envir)$fn[id, idxTest[, r]]) * sqrt(abs(eval(tmp$warp_data, envir =
                                                                                                       call.envir)$fn[id, idxTest[, r]]))
          qn = fdasrvf::f_to_srvf(
            eval(tmp$warp_data, envir = call.envir)$fn[, idxTest[, r]],
            fit$basisInfo$basisConstruct$time
          )
          qn1 = rbind(qn, m_new)
        } else {
          qn1 = eval(tmp$warp_data, envir = call.envir)$fn[, idxTest[, r]]
        }

        time = seq(0, 1, length.out = ncol(Ytest))
        binsize <- mean(diff(time))
        psi = matrix(0, ncol(Ytest), nrow(Ytest))
        vec = matrix(0, ncol(Ytest), nrow(Ytest))
        for (i in 1:nrow(Ytest)) {
          psi[, i] = sqrt(fdasrvf::gradient(
            eval(tmp$warp_data, envir = call.envir)$warping_functions[, idxTest[i, r]],
            binsize
          ))
          vec[, i] <- fdasrvf::inv_exp_map(fit$basisInfo$basisConstruct$mu_psi, psi[, i])
        }

        Ytest = t(rbind(qn1, C * vec))

      } else if (fit$basisInfo$basisType == "jfpcah") {
        C = fit$basisInfo$basisConstruct$C
        id = fit$basisInfo$basisConstruct$id
        srvf = fit$basisInfo$basisConstruct$srvf
        if (srvf) {
          m_new = sign(eval(tmp$warp_data, envir = call.envir)$fn[id, idxTest[, r]]) * sqrt(abs(eval(tmp$warp_data, envir =
                                                                                                       call.envir)$fn[id, idxTest[, r]]))
          qn = fdasrvf::f_to_srvf(
            eval(tmp$warp_data, envir = call.envir)$fn[, idxTest[, r]],
            fit$basisInfo$basisConstruct$time
          )
          qn1 = rbind(qn, m_new)
        } else {
          qn1 = eval(tmp$warp_data, envir = call.envir)$fn[, idxTest[, r]]
        }

        h = fdasrvf::gam_to_h(eval(tmp$warp_data, envir = call.envir)$warping_functions[, idxTest[, r]])

        Ytest = t(rbind(qn1, C * h))

      }
    }
    rmse[r] = sqrt(mean((Ytest - Yhat)^2))
    rSquared[r] = 1 - mean((Ytest - Yhat)^2) / mean((t(Ytest) - fit$basisInfo$Ycenter)^2)

    # Get truncation error for UQ
    if (uqTruncMethod == "gaussian") {
      truncErrorVar = cov(fit$basisInfo$truncError)
      truncError = array(MASS::mvrnorm(prod(dim(preds)[1:2]), rep(0, dim(preds)[3]), truncErrorVar), dim = dim(preds))
    } else if (uqTruncMethod == "empirical") {
      idxResample = sample(nTrain, size = prod(dim(preds)[1:2]), replace = TRUE)
      truncError = aperm(array(t(fit$basisInfo$truncError[idxResample, ]), dim =
                                 dim(preds)[c(3, 1, 2)]), c(2, 3, 1))
    }

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
      coefsResidError[, , k] = rnorm(prod(dim(preds)[1:2]), mean = coefsResidMean, sd =
                                       coefsResidSD)
    }
    residError = array(dim = dim(preds))
    for (idxMCMC in 1:dim(preds)[1]) {
      residError[idxMCMC, , ] = coefsResidError[idxMCMC, , ] %*% fit$basisInfo$basis
    }

    # Get posterior predictive samples
    postPredSamples = preds + truncError + residError

    # Calculate distance from posterior mean
    distBound = numeric(dim(preds)[2])
    for (idx in 1:dim(preds)[2]) {
      distSamples = sqrt(apply((t(
        postPredSamples[, idx, ]
      ) - Yhat[idx, ])^2, 2, mean))
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
    fit_time = fit_time,
    predict_time = pred_time,
    call = match.call()
  )

  return(structure(out, class = 'mvBayesCV'))
}
