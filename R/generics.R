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
.getResidSD = function(object, nCores) {
  UseMethod(".getResidSD", object)
}
.nCoresAdjust = function(object, nCores) {
  UseMethod(".nCoresAdjust", object)
}


traceplot = function(object, ...) {
  UseMethod("traceplot")
}

