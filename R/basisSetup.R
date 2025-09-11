#' @title Basis Expansion Calculation
#'
#' @description Calculates a basis expansion of a specified type.
#' @param Y A matrix of dimension `c(n, nMV)`, where `n` is the number of observations of `nMV` response variables
#' @param basisType The type of basis functions to use. Options are `pca`, `pns`, `jfpca`, `jfpcah`, `splinet`, `bspline`, or `legendre`.
#' @param customBasis Optional user-provided basis of dimension `c(nBasis, nMV)`, only used if `basisType==custom`.
#' @param nBasis integer number of basis components to use (optional). By default, `propVarExplained` is used to choose nBasis.
#' @param propVarExplained Proportion (between 0 and 1) of variation to explain when choosing the number of basis components. Only used if `is.null(nBasis)`.
#' @param center logical: whether or to center `Y` before basis computations.
#' @param scale logical: whether to scale `Y` before basis computations.
#' @return object of class `basisSetup` with plot method.
#' @keywords basis expansion, principal component analysis, dimension reduction
#' @seealso Used by the main function \link{mvBayes}
#' @import stats
#' @export
#'
basisSetup = function(Y,
                      basisType = "pca",
                      customBasis = NULL,
                      nBasis = NULL,
                      propVarExplained = 0.99,
                      center = TRUE,
                      scale = FALSE) {
  out = structure(list(), class = "basisSetup")
  out$Y = Y
  out$nMV = ncol(Y)
  out$basisType = basisType
  out$Ycenter = 0
  out$Yscale = 1
  if (center) {
    out$Ycenter = apply(Y, 2, mean)
  }
  if (scale) {
    out$Yscale = apply(Y, 2, sd)
    out$Yscale[out$Yscale == 0] = 1
  }
  Ystandard = t((t(Y) - out$Ycenter) / out$Yscale)

  if (basisType == "pca") {
    pca = eigen(cov(Ystandard))
    basis = t(pca$vectors)
    coefs = Ystandard %*% t(basis)
    out$varExplained = pca$values
  } else if (basisType == 'pns') {
    is_shapes_available <- requireNamespace("fdasrvf", quietly = TRUE)
    if (!is_shapes_available) {
      stop('install shapes for this basis option')
    }

    # Apply normalization to each column
    pnsdat <- apply(t(Y), 2, normalize_column)
    radius = mean(sqrt(apply(t(Y)^2, 2, sum)))

    # get rough estimate of n.pc
    pca = eigen(cov(t(pnsdat)))
    varExplained.psi = pca$values

    cs.psi = cumsum(varExplained.psi) / sum(varExplained.psi)
    n.pc = which(cs.psi >= 0.99)[1]
    cli::cli_alert_info("Setting n.pc to {n.pc}...")

    basisConstruct = fdasrvf:::fastpns(pnsdat,
                                       n.pc = n.pc,
                                       sphere.type = "small",
                                       output = FALSE)
    coefs = t(basisConstruct$resmat)

    basis = array(0, dim = c(out$nMV, out$nMV))

    out$basisConstruct = basisConstruct
    out$varExplained = basisConstruct$percent / 100
    out$radius = radius

  } else if (basisType %in% c('legendre', 'bspline', 'splinet', 'custom')) {
    if (basisType %in% c('legendre', 'bspline')) {
      if (is.null(nBasis)) {
        stop(paste0(
          "need to specificy nBasis for basisType='",
          basisType,
          "'"
        ))
      }
      if (nBasis < 3) {
        stop("Must have nBasis >= 3 for basisType='", basisType, "'")
      }
      if (basisType == 'legendre') {
        customBasis = basisLegendre(out$nMV, nBasis - 1)
      } else if (basisType == 'bspline') {
        customBasis = t(splines::bs(seq(0, 1, length.out = out$nMV), degree = nBasis))
      }
    } else if (basisType == 'splinet') {
      if (!is.null(nBasis)) {
        stop(paste0(
          "Cannot specificy nBasis for basisType='",
          basisType,
          "'"
        ))
      }
      customBasis = basisSplinet(Ystandard)
      nBasis = nrow(customBasis)
    } else if (basisType == 'custom') {
      if (is.null(customBasis)) {
        stop("Must provide customBasis if basisType=='custom'")
      }
      if (ncol(customBasis) != out$nMV) {
        stop("ncol(customBasis) must match ncol(Y)")
      }
    }
    basisConstruct = customBasisConstruct(customBasis, Ystandard)
    basis = basisConstruct$basis
    coefs = basisConstruct$coefs
    out$varExplained = basisConstruct$varExplained
  } else{
    stop("un-supported basisType")
  }

  propVarCumSum = cumsum(out$varExplained) / sum(out$varExplained)
  if (is.null(nBasis)) {
    nBasis = min(which(propVarCumSum > propVarExplained)[1], out$nMV)
  }
  out$nBasis = nBasis
  out$propVarExplained = propVarCumSum[nBasis]

  out$basis = basis[1:nBasis, , drop = FALSE]
  out$coefs = coefs[, 1:nBasis, drop = FALSE]
  Ytrunc = getYtrunc(out)
  out$truncError = out$Y - Ytrunc

  return(out)
}

.getY.basisSetup = function(object) {
  return(object$Y)
}

#' @export
getYtrunc.basisSetup = function(object,
                                Ytest = NULL,
                                coefs = NULL,
                                nBasis = NULL,
                                ...) {
  if (is.null(coefs)) {
    coefs = getCoefs(object, Ytest)
  }
  if (is.null(nBasis) || (nBasis > object$nBasis)) {
    nBasis = object$nBasis
  }
  if (object$basisType == 'pns') {
    inmat = array(0, dim = c(length(object$basisConstruct$PNS$radii), nrow(coefs)))
    inmat[1:nBasis, ] = t(coefs[, 1:nBasis])
    YtruncStandard = fdasrvf:::fastPNSe2s(inmat, object$basisConstruct) * object$radius
  } else {
    YtruncStandard = coefs %*% object$basis
  }

  if (length(object$Yscale) > 1) {
    YscaleMat = t(replicate(nrow(YtruncStandard), object$Yscale))
  } else {
    YscaleMat = object$Yscale
  }
  if (length(object$Ycenter) > 1) {
    YcentMat = t(replicate(nrow(YtruncStandard), object$Ycenter))
  } else {
    YcentMat = object$Ycenter
  }
  Ytrunc = t(t(as.array(YtruncStandard))) * YscaleMat + YcentMat

  return(Ytrunc)
}

#' @export
getCoefs.basisSetup = function(object, Ytest = NULL) {
  if (is.null(Ytest)) {
    return(object$coefs)
  } else{
    YtestStandard = t((t(Ytest) - object$Ycenter) / object$Yscale)
    return(YtestStandard %*% t(object$basis))
  }
}

#' @export
preprocessY.basisSetup = function(object, Ytest = NULL, ...) {
  if (is.null(Ytest)) {
    return(.getY(object))
  } else{
    return(Ytest)
  }
}

#' @title Visualizing the Basis Decomposition
#'
#' @description Given an object of class "basisSetup" from the basisSetup() function, plot.basisSetup() visualizes the decomposition, colored by basis.
#' @param object An object of class "basisSetup" containing information about a basis decomposition of a response matrix Y.
#' @param nBasis An integer specifying the number of basis functions to plot. If both nBasis and propVarExplained are NULL, object$nBasis is used, but anything <= object$nbasis is allowed.
#' @param propVarExplained Proportion (between 0 and 1) of variation to represent when choosing the number of bases to plot. It is only used if nBasis is NULL. If both are NULL, object$nBasis bases are plotted.
#' @param nPlot A positive integer specifying the number of samples to plot. Default is min(n, 1000), where n is nrow(Xtest) (or nrow(object$X) if Xtest is NULL).
#' @param file An optional location at which the plot will be saved. If NULL, no file is saved.
#' @param title An optional title to be printed at the top of the plot.
#' @param ... additional plot arguments
#' @return Top left: response Y. Top right: Y decomposed and colored by basis. Bottom left: Truncation error due to dimension reduction in the representation of Y. Bottom right: Percent variance explained by each basis.
#' @keywords multivariate Bayesian regression modeling, functional data analysis
#' @seealso See \link{mvBayes}
#' @export
#' @import graphics
#'
plot.basisSetup = function(object,
                           nBasis = NULL,
                           propVarExplained = NULL,
                           nPlot = NULL,
                           idxMV = NULL,
                           xscale = "linear",
                           xlabel = "Multivariate Index",
                           file = NULL,
                           title = NULL,
                           ...) {
  cs = cumsum(object$varExplained) / sum(object$varExplained)

  if (is.null(nBasis) & is.null(propVarExplained)) {
    nBasis = object$nBasis
  }
  if (is.null(nBasis)) {
    nBasis = which(cs >= propVarExplained)[1]
  }
  if (nBasis > length(object$varExplained)) {
    nBasis = object$nBasis
  }

  n = nrow(.getY(object))

  if (is.null(nPlot)) {
    nPlot = min(n, 1000)
  } else if (nPlot > n) {
    warning(
      "nPlot should be at most n, where n=nrow(Xtest) (or n=nrow(object$X) if Xtest is NULL). Using nPlot=n."
    )
    nPlot = n
  }

  if (nPlot < n) {
    idxPlot = sample(n, nPlot)
  } else{
    idxPlot = 1:n
  }

  if (is.null(idxMV)) {
    idxMV = 1:object$nMV
  }

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

  par(mfrow = c(2, 2), mar = c(5, 5, 1, 1), oma=c(0, 0, 2, 0))

  rgbCmap = grDevices::col2rgb('darkblue')
  col = grDevices::rgb(rgbCmap[1] / 255, rgbCmap[2] / 255, rgbCmap[3] / 255, alpha = 0.5)
  matplot(
    idxMV,
    t(.getY(object)[idxPlot, , drop=FALSE]),
    type = 'l',
    col = col,
    lty = 1,
    ylab = "Y",
    xlab = xlabel,
    log = ifelse(xscale == "log", "x", ""),
    cex.lab = 0.9
  )

  basisScaled = vector('list', nBasis)
  for (k in 1:nBasis) {
    basisScaled[[k]] = t(outer(object$coefs[idxPlot, k], object$basis[k, ])) * object$Yscale
  }
  ylim = range(basisScaled, object$truncError)
  for (k in 1:nBasis) {
    rgbCmap = grDevices::col2rgb(cmap[(k - 1) %% 20 + 1])
    col = grDevices::rgb(rgbCmap[1] / 255, rgbCmap[2] / 255, rgbCmap[3] / 255, alpha = 0.5)
    matplot(
      idxMV,
      basisScaled[[k]],
      type = 'l',
      add = (k > 1),
      col = col,
      lty = 1,
      ylim = ylim,
      ylab = "Basis Projection",
      xlab = xlabel,
      log = ifelse(xscale == "log", "x", ""),
      cex.lab = 0.9
    )
  }

  rgbCmap = grDevices::col2rgb('grey')
  col = grDevices::rgb(rgbCmap[1] / 255, rgbCmap[2] / 255, rgbCmap[3] / 255, alpha = 0.5)
  matplot(
    idxMV,
    t(object$truncError[idxPlot, , drop=FALSE]),
    type = 'l',
    col = col,
    lty = 1,
    ylim = ylim,
    ylab = "Truncation Error",
    xlab = xlabel,
    cex.lab = 0.9,
    log = ifelse(xscale == "log", "x", "")
  )

  # Percent variance explained plot
  nMV = ncol(object$Y)
  varTotal = sum(object$varExplained)
  propVarTruncCS = cumsum(object$varExplained[(nBasis + 1):nMV] / varTotal)
  propVarBasis = object$varExplained[1:nBasis] / varTotal
  plot(
    100 * propVarBasis,
    xlim = c(1, nBasis + 1),
    ylim = c(0, 100 * max(propVarBasis, propVarTruncCS, na.rm = TRUE)),
    xaxt = 'n',
    col = cmap[(1:nBasis - 1) %% 20 + 1],
    pch = 19,
    xlab = 'Component',
    ylab = '%Variance',
    main = paste0(sprintf(
      "%.3g",
      100 * object$propVarExplained
    ), "% Variance Explained"),
    cex.main = 1,
    cex.lab = 0.9
  )
  points(
    rep(nBasis + 1, nMV - nBasis),
    100 * propVarTruncCS,
    col = 'grey',
    pch = 19
  )
  axis(
    side = 1,
    at = 1:(nBasis + 1),
    labels = c(1:nBasis, 'T')
  )

  if (!is.null(title)) {
    mtext(title, outer = TRUE, font = 2) # Add title
  }

  if (!is.null(file)) {
    grDevices::dev.copy(grDevices::png, file, ...)
    grDevices::dev.off()
  }
}
