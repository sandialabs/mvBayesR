#' @title Multivariate Bayesian Regression
#'
#' @description Wrapper to fit a multivariate response, using an arbitrary univariate Bayesian regression model to independently fit basis components (e.g., principal components) of the response.
#' @param bayesModel A Bayesian regression model-fitting function, with first argument taking an nxp input matrix or data.frame, and second argument taking an n-vector of numeric responses.
#' @param X A matrix of predictors of dimension nxp, where n is the number of training examples and p is the number of inputs (features).
#' @param Y A response matrix of dimension nxq, where q is the number of multivariate/functional responses.
#' @param basisType The type of basis functions to use. Options are `pca`, `pns`, `splinet`, `bspline`, or `legendre`.
#' @param customBasis Optional user-provided basis of dimension `c(nBasis, nMV)`, only used if `basisType==custom`.
#' @param nBasis An integer specifying the number of basis functions to use. The default is NULL, in which case propVarExplained is used to choose nBasis.
#' @param propVarExplained Proportion (between 0 and 1) of variation to explain when choosing the number of principal components. Only used if nBasis is NULL (which is the default).
#' @param center A logical argument specifying whether or not to center the responses before computing the basis decomposition.
#' @param scale A logical argument specifying whether or not to scale by the standard deviation of each component of the response before computing the basis decomposition.
#' @param nCores An integer less than or equal to nBasis, specifying the number of threads to use when fitting independent Bayesian models.
#' @param samplesExtract function taking the output of bayesModel (`bm`) and extracting posterior samples of all parameters of interest. If `NULL`, mvBayes tries to access`bm$samples`; if unsuccessful, an object called `samples` is created with attribute `residSD`.
#' @param residSDExtract function taking the output of bayesModel (`bm`) and extracting posterior samples of the residual standard deviation (`residSD`). If `NULL``, mvBayes tries to access `bm$samples$residSD`; if unsuccessful, `residSD` is the  standard deviation of the residuals.
#' @param idxSamplesArg str Name of an optional argument of `predict` controlling which posterior samples are used for posterior prediction.
#' @param ... Additional arguments to bayesModel.
#' @details First uses the basisSetup function to decompose the response into nBasis components, then independently fits bayesModel to each of those components.
#' @return An object of class "mvBayes", which is a list containing "X", an object called "basisInfo" of class "basisSetup" containing information about the basis decomposition, "bayesModel", and "bmList", which contains a list of length nBasis containing fitted model objects for each basis component.
#' @keywords multivariate Bayesian regression modeling, functional data analysis
#' @seealso \link{basisSetup} for computing the basis decomposition, \link{predict.mvBayes} for prediction, \link{plot.mvBayes} for plotting the model fit, \link{traceplot} for monitoring posterior convergence, and \link{mvSobol} for sensitivity analysis.
#' @export
#' @import parallel
#' @example inst/examples.R
#'
mvBayes = function(bayesModel,
                   X,
                   Y,
                   basisType = 'pca',
                   customBasis = NULL,
                   nBasis = NULL,
                   propVarExplained = 0.99,
                   center = TRUE,
                   scale = FALSE,
                   nCores = 1,
                   samplesExtract = NULL,
                   residSDExtract = NULL,
                   idxSamplesArg = "idxSamples",
                   ...) {
  out = structure(list(), class = 'mvBayes')

  out$X = X
  out$Y = Y
  out$nMV = ncol(Y)
  out$bayesModel = bayesModel

  basisInfo = basisSetup(
    Y = Y,
    basisType = basisType,
    customBasis = customBasis,
    nBasis = nBasis,
    propVarExplained = propVarExplained,
    center = center,
    scale = scale
  )
  out$basisInfo = basisInfo

  out$samplesExtract = samplesExtract
  out$residSDExtract = residSDExtract
  out$idxSamplesArg = idxSamplesArg

  out = fit(out, nCores, ...)

  return(out)
}

#' @export
fit.mvBayes <- function(object, nCores = 1, ...) {
  nCores <- .nCoresAdjust(object, nCores)

  fitBayesModel <- function(k) {
    tryCatch({
      out <- object$bayesModel(object$X, object$basisInfo$coefs[, k], ...)
      return(out)
    }, error = function(e) {
      message(sprintf("Error fitting model %d: %s", k, e$message))
      return(NULL)
    })
  }

  message(
    sprintf(
      "Starting mvBayes with %d components, using %d cores.",
      object$basisInfo$nBasis,
      nCores
    )
  )

  if (nCores == 1) {
    object$bmList <- lapply(1:object$basisInfo$nBasis, fitBayesModel)
  } else {
    object$bmList <- parallel::mclapply(
      1:object$basisInfo$nBasis,
      fitBayesModel,
      mc.cores = nCores,
      mc.preschedule = FALSE
    )
  }

  object <- .getSamples(object)
  object <- .getResidSD(object, nCores)

  object$nSamples <- length(object$bmList[[1]]$samples$residSD)

  return(object)
}

.nCoresAdjust.mvBayes <- function(object, nCores) {
  nCores <- min(nCores, object$basisInfo$nBasis)
  nCoresAvailable <- parallel::detectCores()

  if (nCores > nCoresAvailable) {
    message(sprintf(
      "Only %d cores are available. Using all available cores.",
      nCoresAvailable
    ))
    nCores <- nCoresAvailable
  }

  return(nCores)
}

.getSamples.mvBayes <- function(object) {
  for (k in 1:length(object$bmList)) {
    bm <- object$bmList[[k]]
    if (is.null(object$samplesExtract)) {
      if (!("samples" %in% names(bm))) {
        if (k == 1) {
          message("Generating 'samples' attribute, since it was absent in 'bmList[[1]]'")
        }
        bm$samples <- list()
        next
      }
    } else {
      bm$samples <- object$samplesExtract(bm)
    }
    object$bmList[[k]] <- bm
  }
  return(object)
}

.getResidSD.mvBayes <- function(object, nCores) {
  if (is.null(object$residSDExtract)) {
    if (!("residSD" %in% names(object$bmList[[1]]$samples))) {
      message("Approximating 'residSD', since 'residSDExtract' is None")
      postCoefs <- predict(object,
                           object$X,
                           nCores = nCores,
                           returnPostCoefs = TRUE)$postCoefs
      for (k in 1:length(object$bmList)) {
        resid <- object$basisInfo$coefs[, k] - t(postCoefs[, , k])
        object$bmList[[k]]$samples$residSD <- apply(resid, 2, sd)
      }
    }
  } else {
    for (k in 1:length(object$bmList)) {
      object$bmList[[k]]$samples$residSD <- object$residSDExtract(object$bmList[[k]])
    }
  }
  return(object)
}



#' @title Posterior Predictive Samples of a Multivariate Response
#'
#' @description Given an object of class "mvBayes" from the mvBayes() function, predict.mvBayes() obtains posterior predictive samples at specific input locations. Residual variance from the basis truncation can (optionally) be inserted as well.
#' @param object An object of class "mvBayes" containing the multivariate Bayesian model fit of a response matrix Y.
#' @param Xtest A matrix or data.frame of inputs at which predictions are desired.
#' @param idxSamples A str describing which samples to use
#' @param addResidError A logical to add back in residual error
#' @param addTruncError A logical indicating whether or not to insert random truncation error into the predictions, where the truncation comes from the reduced dimensionality in the representation of Y.
#' @param returnPostCoefs A logical indicating whether or not to output predictions of Ynew (i.e., the observation-specific basis weighting) in addition to predictions of Y.
#' @param nCores An integer less than or equal to object$basisInfo$nBasis, specifying the number of threads to use when obtaining predictions from the independent Bayesian models.
#' @param idxSamplesArg str Name of an optional argument of `predict` controlling which posterior samples are used for posterior prediction.
#' @param ... Additional arguments to predict.bayesModel, where "bayesModel" is the Bayesian model used in fitting, specified in object$bayesModel.
#' @return If getpostCoefs==FALSE, predict.mvBayes() outpus an array of dimension c(nSamples, n, q), where nSamples is the number of posterior samples obtained in fitting bayesModel, and n and q are respectively the number of rows and columns in Y. Elements of this array are posterior predictive samples of Y. Otherwise, a list of length two is output, with elements Ypost giving samples from the posterior predictive distribution of Y, and postCoefs givine samples from the posterior predictive distribution of Ynew, i.e., from the observation-specific basis weights.
#' @keywords multivariate bayesian regression modeling, functional data analysis
#' @seealso See \link{mvBayes}
#' @export
#'
predict.mvBayes = function(object,
                           Xtest,
                           idxSamples = "default",
                           addResidError = FALSE,
                           addTruncError = FALSE,
                           returnPostCoefs = FALSE,
                           nCores = 1,
                           idxSamplesArg = NULL,
                           ...) {
  n = nrow(Xtest)
  q = object$basisInfo$nMV
  if (is.null(idxSamplesArg)) {
    idxSamplesArg = object$idxSamplesArg
  }

  predictArgs = list(...)

  if (idxSamples[1] == "default") {
    # Do nothing
  } else if (!(idxSamplesArg %in% names(formals(
    utils::getS3method("predict", class(object$bmList[[1]])[1])
  )))) {
    warning(
      paste0(
        "'",
        idxSamplesArg,
        "' is not an argument of the bayesModel predict function...setting idxSamples='default'"
      )
    )
    idxSamples = "default"
  } else {
    if (idxSamples[1] == "final") {
      idxSamples = object$nSamples
    } else if (!is.numeric(idxSamples)) {
      stop("'idxSamples' must be 'default', 'final', or numeric.")
    }
    predictArgs[[idxSamplesArg]] = idxSamples
  }

  nCores = .nCoresAdjust(object, nCores)

  ntest = nrow(Xtest)
  nBasis = object$basisInfo$nBasis
  nMV = object$basisInfo$nMV
  basis = object$basisInfo$basis

  predictBayesModel = function(k) {
    return(do.call(predict, c(
      list(object$bmList[[k]], Xtest), predictArgs
    )))
  }

  if (nCores == 1) {
    postCoefs = lapply(1:nBasis, predictBayesModel)
  } else {
    postCoefs = parallel::mclapply(1:nBasis,
                                   predictBayesModel,
                                   mc.cores = nCores,
                                   mc.preschedule = FALSE)
  }
  postCoefs = matrix(unlist(postCoefs), ncol = nBasis)

  if (addResidError) {
    residError <- sapply(1:nBasis, function(k) {
      rnorm(nrow(postCoefs), sd = object$bmList[[k]]$samples$residSD)
    })
    postCoefs = postCoefs + residError
  }

  nSamples = nrow(postCoefs) / ntest
  postCoefs <- array(postCoefs, dim = c(nSamples, ntest, nBasis))

  if (object$basisInfo$basisType == "pns") {
    PNS = object$basisInfo$basisConstruct
    N = dim(postCoefs)[1] * dim(postCoefs)[2]
    nBasis = object$basisInfo$nBasis
    inmat = matrix(0, length(PNS$PNS$radii), N)
    inmat[1:nBasis, ] = t(array(postCoefs, dim = c(N, nBasis)))
    YtruncStandard = fdasrvf:::fastPNSe2s(inmat, PNS)
    if (nSamples == 1){
      YpostT = array(YtruncStandard, dim = c(n, q)) * object$basisInfo$radius
      YpostT = t(YpostT)
    } else {
      YpostT = array(YtruncStandard, dim = c(nSamples, n, q)) * object$basisInfo$radius
      YpostT = aperm(YpostT, c(3, 2, 1))
    }

    rm(YtruncStandard)
  } else {
    if (nSamples == 1){
      if (n == 1){
        YstandardPostT =  t(basis) %*% t(t(postCoefs[1, , ]))
      } else {
        YstandardPostT =  t(basis) %*% t(postCoefs[1, , ])
      }
    } else {
      YstandardPostT = array(0, dim = c(q, n, nSamples))
      for (it in 1:nSamples){
        if (n == 1){
          YstandardPostT[, , it] =  t(basis) %*% t(t(postCoefs[it, , ]))
        } else {
          YstandardPostT[, , it] =  t(basis) %*% t(postCoefs[it, , ])
        }

      }
    }
    YpostT = YstandardPostT * object$basisInfo$Yscale + object$basisInfo$Ycenter
  }

  if (nSamples == 1){
    YpostT = t(YpostT)
  } else {
    YpostT = aperm(YpostT, 3:1)
  }

  if (addTruncError) {
    idx = sample(nrow(object$basisInfo$truncError),
                 size = nSamples * ntest,
                 replace = TRUE)
    if (nSamples == 1){
      truncError = array(object$basisInfo$truncError[idx, ], dim = c(ntest, nMV))
    } else {
      truncError = array(object$basisInfo$truncError[idx, ], dim = c(nSamples, ntest, nMV))
    }

    YpostT = YpostT + truncError
  }

  if (returnPostCoefs) {
    return(list(Ypost = YpostT, postCoefs = postCoefs))
  } else {
    return(YpostT)
  }
}


#' @title Traceplots from the Bayesian Model Fit of a Multivariate Response
#'
#' @description Given an object of class "mvBayes" from the mvBayes() function, traceplot() plots traceplots of user-specified variables from the Bayesian model fit, colored by basis.
#' @param object An object of class "mvBayes" containing the multivariate Bayesian model fit of a response matrix Y.
#' @param modelParams A character vector of length <= 8 containing the names of the parameters for which traceplots will be made. If NULL, selects "plottable" attributes of object$bmList[[1]].
#' @param labels A character vector of the same length as modelParams, containing the y-axis labels that will be printed for each traceplot. The default is to simply use "modelParams".
#' @param file An optional location at which the traceplots will be saved. If NULL, no file is saved.
#' @param title An optional title to be printed at the top of the traceplots.
#' @param ... additional plot arguments
#' @keywords multivariate bayesian regression modeling, functional data analysis
#' @seealso See \link{mvBayes}
#' @export
#' @import graphics
#'
traceplot = function(object,
                     modelParams = NULL,
                     labels = NULL,
                     title = NULL,
                     file = NULL,
                     ...) {
  # Extract bmList and nBasis from the object
  nBasis <- object$basisInfo$nBasis

  # Default parameter selection
  if (is.null(modelParams)) {
    isModelParam <- function(obj) {
      if (is.null(obj) || is.function(obj) || is.character(obj)) {
        return(FALSE)
      }
      tryCatch({
        objUnlist <- unlist(obj)
        return(length(objUnlist) == object$nSamples)
      }, error = function(e) {
        return(FALSE)
      })
    }
    modelParams <- names(object$bmList[[1]]$samples)
    modelParams <- modelParams[sapply(modelParams, function(p)
      isModelParam(object$bmList[[1]]$samples[[p]]))]
  }

  if (is.null(labels)) {
    labels <- modelParams
  }

  nParams <- length(modelParams)
  if (nParams > 8) {
    warning("Currently, must have length(modelParams) <= 8. Plotting the first 8.")
    modelParams <- modelParams[1:8]
    labels <- labels[1:8]
    nParams <- 8
  }

  nrow <- ceiling(nParams / 2)
  ncol <- ifelse(nParams == 1, 1, 2)

  cmap = c(
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#aec7e8",
    "#ffbb78",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#c49c94",
    "#f7b6d2",
    "#c7c7c7",
    "#dbdb8d",
    "#9edae5"
  )

  par(mfrow = c(nrow, ncol), mar = c(5, 5, 1, 1), oma=c(0, 0, 2, 0))

  for (j in 1:length(modelParams)) {
    param_list <- lapply(object$bmList, function(bm)
      bm$samples[[modelParams[j]]])
    ylim <- range(unlist(param_list))

    plot(
      param_list[[1]],
      ylim = ylim,
      xlab = 'MCMC Iteration',
      ylab = labels[j],
      col = cmap[1],
      type = 'l'
    )
    if (nBasis > 1) {
      for (k in 2:nBasis) {
        lines(param_list[[k]], col = cmap[k])
      }
    }
  }

  if (!is.null(title)) {
    mtext(title, outer = TRUE, font = 2) # Add title
  }

  if (!is.null(file)) {
    grDevices::dev.copy(grDevices::png, file, ...)
    grDevices::dev.off()
  }
}


#' @title Plot the Bayesian Model Fit of a Multivariate Response
#'
#' @description Given an object of class "mvBayes" from the mvBayes() function, plot.mvBayes() plots a few aspects of the Bayesian model fit, colored by basis.
#' @param object An object of class "mvBayes" containing the multivariate Bayesian model fit of a response matrix Y.
#' @param Xtest A matrix or data.frame of inputs at which the model fit will be evaluated. If NULL, the training inputs X will be used.
#' @param Ytest A matrix of responses for which the model fit will be evaluated. Should correspond to Xtest. If NULL, the training responses Y will be used.
#' @param idxSamples str which samples to use
#' @param nPlot A positive integer specifying the number of samples to plot. Default is min(n, 1000), where n is nrow(Xtest) (or nrow(object$X) if Xtest is NULL).
#' @param idxMV A vector describing the sample points
#' @param xscale 'linear' or 'log' of x-axis
#' @param xlabel str for x-axis label
#' @param title An optional title to be printed at the top of the diagnostic plots.
#' @param file An optional location at which the plot will be saved. If NULL, no file is saved.
#' @param ... additional plot arguments
#' @return Top left: Residuals after applying mvBayes(), plotted on top of Ytest (centered if centered=TRUE in the mvBayes() call). Top right: Residuals decomposed and colored by basis. Bottom left: R^2 values for each component. Bottom right: %Residual variance attribution for each component.
#' @keywords multivariate bayesian regression modeling, functional data analysis
#' @seealso See \link{mvBayes}, \link{predict.mvBayes}
#' @export
#' @import graphics
#'
plot.mvBayes <- function(object,
                         Xtest = NULL,
                         Ytest = NULL,
                         idxSamples = "final",
                         nPlot = NULL,
                         idxMV = NULL,
                         xscale = "linear",
                         xlabel = "Multivariate Index",
                         title = NULL,
                         file = NULL,
                         ...) {
  # Check if training data should be used
  useXtrain <- is.null(Xtest)
  useYtrain <- is.null(Ytest)

  if (useXtrain) {
    Xtest <- object$X
  } else if (useYtrain) {
    message("Model output at user-specified Xtest is being compared to training responses Y.")
  }

  if (useYtrain) {
    Ytest <- .getY(object$basisInfo)
    coefs <- object$basisInfo$coefs
    truncError <- object$basisInfo$truncError
  } else {
    if (useXtrain) {
      message("Model output at training inputs X is being compared to user-specified Ytest.")
    }

    coefs <- getCoefs(object$basisInfo, Ytest)
    Ytest <- preprocessY(object$basisInfo, Ytest)
    Ytrunc <- getYtrunc(object$basisInfo, coefs = coefs)
    truncError <- Ytest - Ytrunc
  }

  if (is.null(nPlot)) {
    nPlot <- min(nrow(Xtest), 1000)
  } else if (nPlot > nrow(Xtest)) {
    warning("nPlot should be at most n=len(Xtest) (or n=len(X) if Xtest is NULL). Using nPlot=n.")
    nPlot <- nrow(Xtest)
  }
  idxPlot <- sample(nrow(Xtest), nPlot, replace = FALSE)

  if (is.null(idxMV)) {
    idxMV <- 1:object$basisInfo$nMV
  }

  if (object$basisInfo$basisType == "pns") {
    Ycentered = t(t(Ytest) - colMeans(Ytest))
  } else {
    Ycentered = t(t(Ytest) - object$basisInfo$Ycenter)
  }

  # Get predictions and residuals
  YpostObj <- predict(object,
                      Xtest,
                      returnPostCoefs = TRUE,
                      idxSamples = idxSamples)
  Ypost <- YpostObj$Ypost
  postCoefs <- YpostObj$postCoefs

  if (ndims(Ypost)==3){
    Ypred <- apply(Ypost, 2:3, mean)
  } else {
    Ypred <- Ypost
  }
  coefsPred <- apply(postCoefs, 2:3, mean)

  R <- Ytest - Ypred
  RbasisCoefs <- coefs - coefsPred

  cmap = c(
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#aec7e8",
    "#ffbb78",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#c49c94",
    "#f7b6d2",
    "#c7c7c7",
    "#dbdb8d",
    "#9edae5"
  )

  # Create plot
  par(mfrow = c(2, 2), mar = c(5, 5, 1, 1), oma = c(0, 0, 2, 0))

  # Residuals plot
  mseOverall <- mean(R^2)
  if (length(object$basisInfo$Ycenter) == 1 &
      length(object$basisInfo$Yscale) == 1) {
    legendLab = "Y"
  } else{
    legendLab = "Y (Centered)"
  }
  matplot(
    idxMV,
    t(Ycentered[idxPlot, , drop=FALSE]),
    type = "l",
    lty = 1,
    ylim = range(Ycentered, R),
    col = grDevices::rgb(0.7, 0.7, 1.0, alpha = 0.5),
    xlab = xlabel,
    ylab = "Residuals",
    log = ifelse(xscale == "log", "x", ""),
    cex.lab = 0.9
  )
  mtext(
    sprintf("Overall MSE = %.4g", mseOverall),
    cex=0.85
  )
  matlines(idxMV,
           t(R),
           lty = 1,
           col = grDevices::rgb(0, 0, 0, alpha = 0.5))
  legend(
    "bottomleft",
    legend = c("Residual", legendLab),
    col = c(grDevices::rgb(0, 0, 0), grDevices::rgb(0.7, 0.7, 1.0)),
    lwd=2,
    cex = 0.65,
    y.intersp = 0.75,
    text.width=0.25*diff(range(idxMV)),
    bg = grDevices::rgb(1, 1, 1, 0.2)
  )

  # Residual decomposition plot
  RbasisScaled <- list()
  mseBasis <- numeric(object$basisInfo$nBasis)
  varBasis <- numeric(object$basisInfo$nBasis)
  if (object$basisInfo$basisType == "pns") {
    mseTrunc = sum(object$basisInfo$varExplained[object$basisInfo$nBasis:object$basisInfo$nMV])
    for (k in 1:object$basisInfo$nBasis) {
      PNS = object$basisInfo$basisConstruct
      N = length(RbasisCoefs[, k])
      inmat = matrix(0, length(PNS$PNS$radii), N)
      inmat[k, ] = coefsPred[, k]
      basisScaled = as.matrix(fdasrvf:::fastPNSe2s(inmat, PNS)) * object$basisInfo$radius
      RbasisScaled[[k]] = Ytest - basisScaled
      mseBasis[k] = mean(RbasisCoefs[, k]^2)
      varBasis[k] = mean(coefs[, k]^2)
    }
  } else {
    mseTrunc = sum(object$basisInfo$varExplained[object$basisInfo$nBasis:object$basisInfo$nMV])*(nrow(Ytest)-1)/nrow(Ytest)
    for (k in 1:object$basisInfo$nBasis) {
      RbasisScaled[[k]] = t(t(outer(
        RbasisCoefs[idxPlot, k], object$basisInfo$basis[k, ]
      )) * object$basisInfo$Yscale) # all residuals
      mseBasis[k] = mean(RbasisScaled[[k]]^2)

      basisScaled = t(t(outer(coefs[, k], object$basisInfo$basis[k, ])) * object$basisInfo$Yscale)
      varBasis[k] = mean(basisScaled^2)
    }
  }

  mseOrder <- order(mseBasis, decreasing = TRUE)

  rgbCmap = grDevices::col2rgb('grey')
  col = grDevices::rgb(rgbCmap[1] / 255, rgbCmap[2] / 255, rgbCmap[3] / 255, alpha = 0.5)
  matplot(
    idxMV,
    t(truncError[idxPlot, , drop=FALSE]),
    type = "l",
    lty = 1,
    col = col,
    ylim = range(RbasisScaled, truncError),
    xlab = xlabel,
    ylab = "Basis Projection",
    log = ifelse(xscale == "log", "x", ""),
    cex.lab = 0.9
  )
  for (k in 1:object$basisInfo$nBasis) {
    rgbCmap = grDevices::col2rgb(cmap[(k - 1) %% 20 + 1])
    col = grDevices::rgb(rgbCmap[1] / 255, rgbCmap[2] / 255, rgbCmap[3] / 255, alpha = 0.5)
    matlines(
      idxMV,
      t(RbasisScaled[[k]]),
      type = 'l',
      col = col,
      lty = 1,
      xlab = xlabel,
      log = ifelse(xscale == 'log', 'x', '')
    )
  }

  # R^2 plot
  r2Basis <- 1 - mseBasis / varBasis
  varOverall <- mean(Ycentered^2)
  r2Overall <- 1 - mseOverall / varOverall

  plot(
    1:object$basisInfo$nBasis,
    r2Basis,
    pch = 19,
    xlab = "Component",
    ylab = expression(R^2),
    col = cmap[(1:object$basisInfo$nBasis - 1) %% 20 + 1],
    cex.lab = 0.9
  )
  mtext(
    latex2exp::TeX(sprintf("Overall $R^2$ = %.3g", r2Overall)),
    outer=FALSE,
    cex=0.85
  )
  axis(
    side = 1,
    at = 1:(object$basisInfo$nBasis + 1),
    labels = c(1:object$basisInfo$nBasis, 'T')
  )
  abline(h = r2Overall,
         col = "grey",
         lty = 2)

  # Residual variance plot
  mseTruncProp <- (object$basisInfo$varExplained[(object$basisInfo$nBasis + 1):object$basisInfo$nMV] /
                     sum(object$basisInfo$varExplained[(object$basisInfo$nBasis + 1):object$basisInfo$nMV]))
  mseTruncCS <- cumsum(mseTruncProp * mseTrunc / mseOverall)
  mseExplained = mseBasis / mseOverall

  plot(
    1:object$basisInfo$nBasis,
    100 * mseExplained,
    pch = 19,
    xlab = "Component",
    ylab = "%Residual Variance",
    col = cmap[(1:object$basisInfo$nBasis - 1) %% 20 + 1],
    cex.lab = 0.9
  )
  points(
    rep(object$basisInfo$nBasis + 1, length(mseTruncCS)),
    100 * mseTruncCS,
    col = "grey",
    pch = 19
  )

  if (!is.null(title)) {
    mtext(title, outer = TRUE, font = 2) # Add title
  }

  # Save or display plot
  if (!is.null(file)) {
    grDevices::dev.copy(grDevices::png, file, ...)
    grDevices::dev.off()
  }
}
