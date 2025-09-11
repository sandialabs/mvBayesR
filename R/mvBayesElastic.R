#' @title Multivariate Elastic FDA Bayesian Regression
#'
#' @description Wrapper to fit a elastic functional response, using an arbitrary univariate Bayesian regression model to independently fit basis components (e.g., principal components) of the response.
#' @param bayesModel A Bayesian regression model-fitting function, with first argument taking an nxp input matrix or data.frame, and second argument taking an n-vector of numeric responses.
#' @param X A matrix of predictors of dimension nxp, where n is the number of training examples and p is the number of inputs (features).
#' @param Y A response matrix of dimension nxq, where q is the number of multivariate/functional responses.
#' @param warpData `time_warping` object from `fdasrvf` if `basisType=jfpca` or `jfpcah`
#' @param basisType The type of basis functions to use. Options are `jfpca`, `jfpcah`.
#' @param nBasis An integer specifying the number of basis functions to use. The default is NULL, in which case propVarExplained is used to choose nBasis.
#' @param propVarExplained Proportion (between 0 and 1) of variation to explain when choosing the number of principal components. Only used if nBasis is NULL (which is the default).
#' @param nCores An integer less than or equal to nBasis, specifying the number of threads to use when fitting independent Bayesian models.
#' @param samplesExtract function taking the output of bayesModel (`bm`) and extracting posterior samples of all parameters of interest. If `NULL`, mvBayes tries to access`bm$samples`; if unsuccessful, an object called `samples` is created with attribute `residSD`.
#' @param residSDExtract function taking the output of bayesModel (`bm`) and extracting posterior samples of the residual standard deviation (`residSD`). If `NULL``, mvBayes tries to access `bm$samples$residSD`; if unsuccessful, `residSD` is the  standard deviation of the residuals.
#' @param idxSamplesArg str Name of an optional argument of `predict` controlling which posterior samples are used for posterior prediction.
#' @param srvf use SRVF if `basisType=jfpca` or `jfpcah` default = TRUE
#' @param idx vector of indices to subset `warpData`
#' @param ... Additional arguments to bayesModel.
#' @details First uses the basisSetup function to decompose the response into nBasis components, then independently fits bayesModel to each of those components.
#' @return An object of class "mvBayes", which is a list containing "X", an object called "basisInfo" of class "basisSetup" containing information about the basis decomposition, "bayesModel", and "bmList", which contains a list of length nBasis containing fitted model objects for each basis component.
#' @keywords multivariate bayesian regression modeling, functional data analysis
#' @seealso \link{basisSetup} for computing the basis decomposition, \link{predict.mvBayes} for prediction, \link{plot.mvBayes} for plotting the model fit, \link{traceplot} for monitoring posterior convergence, and \link{mvSobol} for sensitivity analysis.
#' @export
#' @import parallel
#' @example inst/examples.R
#'
mvBayesElastic = function(bayesModel,
                          X,
                          Y,
                          warpData = NULL,
                          basisType = 'jfpcah',
                          nBasis = NULL,
                          propVarExplained = 0.99,
                          nCores = 1,
                          samplesExtract = NULL,
                          residSDExtract = NULL,
                          idxSamplesArg = "idxSamples",
                          srvf = FALSE,
                          idx = NULL,
                          ...) {
  out = structure(list(), class = 'mvBayesElastic')

  if (!is.null(warpData) && !is.null(idx)) {
    warpData$fn = warpData$fn[, idx]
    warpData$f0 = warpData$f0[, idx]
    warpData$q0 = warpData$q0[, idx]
    warpData$qn = warpData$qn[, idx]
    warpData$warping_functions = warpData$warping_functions[, idx]
  }

  out$X = X
  out$Y = Y
  out$nMV = ncol(Y)
  out$bayesModel = bayesModel

  basisInfo = basisSetupElastic(
    Y = Y,
    warpData = warpData,
    basisType = basisType,
    nBasis = nBasis,
    propVarExplained = propVarExplained,
    srvf = srvf
  )
  out$basisInfo = basisInfo

  out$samplesExtract = samplesExtract
  out$residSDExtract = residSDExtract
  out$idxSamplesArg = idxSamplesArg

  out = fit(out, nCores, ...)

  return(out)
}

.nCoresAdjust.mvBayesElastic <- function(object, nCores) {
  return(.nCoresAdjust.mvBayes(object, nCores))
}

.getSamples.mvBayesElastic <- function(object) {
  return(.getSamples.mvBayes(object))
}

.getResidSD.mvBayesElastic <- function(object, nCores) {
  return(.getResidSD.mvBayes(object, nCores))
}

#' @export
fit.mvBayesElastic = function(object, nCores = 1, ...) {
  return(fit.mvBayes(object, nCores, ...))
}

#' @export
predict.mvBayesElastic = function(object,
                                  Xtest,
                                  idxSamples = "all",
                                  addResidError = FALSE,
                                  addTruncError = FALSE,
                                  returnPostCoefs = FALSE,
                                  nCores = 1,
                                  ...) {
  return(
    predict.mvBayes(
      object,
      Xtest,
      idxSamples,
      addResidError,
      addTruncError,
      returnPostCoefs,
      nCores,
      ...
    )
  )
}

#' @export
plot.mvBayesElastic = function(x,
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
  return(plot.mvBayes(x,
                      Xtest = NULL,
                      Ytest = NULL,
                      idxSamples = "final",
                      nPlot = NULL,
                      idxMV = NULL,
                      xscale = "linear",
                      xlabel = "Multivariate Index",
                      title = NULL,
                      file = NULL,
                      ...))
}
