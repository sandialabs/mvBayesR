# for basisSetup and its descendants
.getY = function(object) {
  UseMethod(".getY", object)
}
getYtrunc = function(object,
                     Ytest = NULL,
                     coefs = NULL,
                     ...) {
  # Generic getYtrunc function
  UseMethod("getYtrunc", object)
}
getCoefs = function(object, Ytest = NULL) {
  UseMethod("getCoefs", object)
}
preprocessY = function(object, ...) {
  UseMethod("preprocessY", object)
}


# for mvBayes and its descendants
fit = function(object, ...) {
  UseMethod("fit", object)
}
.getSamples = function(object) {
  UseMethod(".getSamples", object)
}
.getResidSD = function(object) {
  UseMethod(".getResidSD", object)
}
.nCoresAdjust = function(object, nCores) {
  UseMethod(".nCoresAdjust", object)
}


#' @title Traceplots from a Bayesian Model Fit
#'
#' @description Generic function for plotting traceplots of the parameters of a
#'   fitted Bayesian model.
#' @param object A fitted model object, e.g. of class "mvBayes".
#' @param ... Additional arguments passed to methods.
#' @return no return value
#' @seealso \link{traceplot.mvBayes}
#' @export
#'
traceplot = function(object, ...) {
  UseMethod("traceplot")
}
