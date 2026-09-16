#' @title Elastic FDA Basis Expansion Calculation
#'
#' @description Calculates a basis expansion of a functional data. Used by the mvBayes function to calculate a basis expansion of a response matrix.
#' @param Y A matrix of dimension nxq, where n is the number of observations of q multivariate/functional variables
#' @param warpData `time_warping` object from `fdasrvf` if `basisType=jfpca` or `jfpcah`
#' @param basisType The type of basis functions to use. Options are `jfpca` or `jfpcah`.
#' @param nBasis An integer specifying the number of basis functions to use. The default is NULL, in which case propVarExplained is used to choose nBasis.
#' @param propVarExplained Proportion (between 0 and 1) of variation to explain when choosing the number of principal components. Only used if nBasis is NULL (which is the default).
#' @param srvf use SRVF if `basisType=jfpca` or `jfpcah` default = TRUE
#' @return An object of class "basisSetup" containing information about the basis decomposition: "Y" is the original response matrix, "Ycenter" and "Yscale" are vectors of length q specifying the centering and scaling values used for each dimension of Y, "basisType" is the user-specified type of basis, the q-vector "varExplained" specifies the amount of variance explained by each basis in the expansion, matrix "basis" is the qxq basis matrix (e.g., the principal component vectors in the case of pca), "coefs" is the nxq matrix of observation-specific basis weights (e.g., the "scores" in the case of pca), and "truncError" is the nxq matrix of residuals, after accounting for the first nBasis bases.
#' @seealso Used by the main function \link{mvBayes}
#' @import stats
#' @export
#'
basisSetupElastic = function(Y,
                             warpData = NULL,
                             basisType = "jfpcah",
                             nBasis = NULL,
                             propVarExplained = 0.99,
                             srvf = FALSE) {
  if (!requireNamespace("fdasrvf", quietly = TRUE)) {
    stop("The 'fdasrvf' package is not installed. Cannot use basisSetupElastic.")
  }

  out = structure(list(), class = 'basisSetupElastic')

  if (is.null(warpData)) {
    warning("warpData not provided. Computing default warpData from Y")
    idxMV = seq(0, 1, length = ncol(Y))
    warpData = fdasrvf::time_warping(t(Y), idxMV)
  }

  out$warpData = warpData
  out$basisType = basisType
  out$Y = Y
  out$Yscale = 1

  if (basisType == 'jfpca') {
    basisConstruct = fdasrvf::jointFPCA(
      warpData,
      no = length(warpData$time),
      srvf = srvf,
      showplot = FALSE
    )
    out$varExplained = basisConstruct$eigs

    basis = basisConstruct$U
    out$varExplained = basisConstruct$eigs
    out$Ycenter = basisConstruct$mu_g
    out$YjointElastic = t(basisConstruct$g)
  } else if (basisType == 'jfpcah') {
    if (!is.null(nBasis)) {
      propVarExplained = 0.999 # would set to 1.0, but this breaks fdasrsf
    }
    basisConstruct = fdasrvf::jointFPCAh(warpData,
                                         var_exp = propVarExplained,
                                         srvf = srvf,
                                         showplot = FALSE)
    out$varExplained = basisConstruct$eigs

    no_q = ncol(basisConstruct$U)
    Psi_q = basisConstruct$U %*% basisConstruct$Uz[1:no_q, ]
    Psi_h = basisConstruct$U1 %*% basisConstruct$Uz[(no_q + 1):nrow(basisConstruct$Uz), ]
    basis = rbind(Psi_q, Psi_h)
    rm(Psi_q, Psi_h)
    out$varExplained = basisConstruct$eigs
    out$Ycenter = c(basisConstruct$mqn, basisConstruct$mh)
    out$YjointElastic = t(rbind(basisConstruct$qn1, basisConstruct$C * basisConstruct$h))
  } else{
    stop("un-supported basisType")
  }

  propVarCumSum = cumsum(out$varExplained) / sum(out$varExplained)
  nBasisAvailable = min(ncol(basis), ncol(basisConstruct$coef))
  if (is.null(nBasis)) {
    if (basisType == 'jfpcah') {
      # jointFPCAh already truncated to propVarExplained
      nBasis = nBasisAvailable
    } else {
      nBasis = which(propVarCumSum >= propVarExplained)[1]
      if (is.na(nBasis)) {
        nBasis = nBasisAvailable
      }
      nBasis = min(nBasis, nBasisAvailable)
    }
  } else if (nBasis > nBasisAvailable) {
    warning(
      paste0(
        "User-specified 'nBasis' larger than the number of ",
        basisType,
        " bases. Setting nBasis=",
        nBasisAvailable,
        "."
      )
    )
    nBasis = nBasisAvailable
  }

  out$nBasis = nBasis
  out$propVarExplained = propVarCumSum[nBasis]

  out$basisConstruct = basisConstruct
  coefs = basisConstruct$coef
  out$nMV = ncol(out$YjointElastic)

  # basis is stored with one row per component, matching basisSetup(), so that
  # getYtrunc() can compute coefs %*% basis
  out$basis = t(basis[, 1:nBasis, drop = FALSE])
  out$coefs = coefs[, 1:nBasis, drop = FALSE]
  YjointElasticTrunc = getYtrunc(out)
  out$truncError = out$YjointElastic - YjointElasticTrunc

  return(out)
}

.getY.basisSetupElastic = function(object) {
  return(object$YjointElastic)
}

#' @export
getYtrunc.basisSetupElastic = function(object,
                                       Ytest = NULL,
                                       coefs = NULL,
                                       ...) {
  return(getYtrunc.basisSetup(object, Ytest, coefs))
}

#' @export
getCoefs.basisSetupElastic = function(object, Ytest = NULL) {
  if (is.null(Ytest)) {
    return(object$coefs)
  } else{
    stop(
      "Projecting new responses onto an elastic basis is not supported: 'fdasrvf' provides no projection method for jfpca/jfpcah bases. Re-run the elastic alignment on the new data instead.",
      call. = FALSE
    )
  }
}

#' @export
preprocessY.basisSetupElastic = function(object,
                                         Ytest = NULL,
                                         ...) {
  if (is.null(Ytest)) {
    return(.getY(object))
  } else{
    stop(
      "Preprocessing new responses with an elastic basis is not supported: 'fdasrvf' provides no projection method for jfpca/jfpcah bases. Re-run the elastic alignment on the new data instead.",
      call. = FALSE
    )
  }
}

#' @export
plot.basisSetupElastic = function(x, ...) {
  plot.basisSetup(x, ...)
}
