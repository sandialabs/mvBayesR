#' @title Sobol' Indices for Multivariate Response
#'
#' @description Given an object of class "mvBayes" from the mvBayes() function, mvSobol() calculates the Sobol' indices for each input variable, assuming the response is multivariate. Currently only works if object$bayesModel is the BASS::bass() function, or if object$bayesModel is compatible with the anova() function.
#' @param object An object of class "mvBayes" containing the multivariate Bayesian model fit for a specific basis expansion of a response matrix Y.
#' @param idxSamples Use final MCMC sample
#' @param nMC Number of Monte Carlo iterations to compute, if None will attempt `anova` if available. Parameter is log2( n ), i.e., it makes 2^n points
#' @param totalSobol whether to compute the total Sobol' index (summing all interactions)
#' @param showPlot whether to compute generate a plot along with the returned Sobol' object
#' @param ... Additional arguments to BASS::sobolBasis (other emulator-specific Sobol' calculations to be implemented in the future). Unnecessary for anova-compatible models.
#' @return An object of class "bassSob" if object$bayesModel is the BASS::bass() function, or "mvSobol" otherwise. Contains information about the Sobol' decomposition: See ?BASS::sobolBasis for more info.
#' @seealso See \link{mvBayes}
#' @export
#' @import BASS
#'
mvSobol = function (object,
                  totalSobol = TRUE,
                  idxSamples = "final",
                  nMC = NULL,
                  showPlot = FALSE,
                  ...) {
  if (!(methods::is(object, 'mvBayes') ||
        methods::is(object, 'mvBayesElastic'))) {
    stop("'object' must be of class 'mvBayes' or 'mvBayesElastic'")
  }

  predictArgs = list(...)

  if (object$basisInfo$basisType == "pns" && (is.null(nMC))) {
    nMC = 2^12
  }

  if (!is.null(nMC)) {
    if (length(nMC) != 1 || !is.numeric(nMC) || nMC < 2) {
      stop("'nMC' must be a single number >= 2.")
    }
    # rsobol() generates 2^m points, so nMC must be a power of two
    m = round(log2(nMC))
    if (!isTRUE(all.equal(2^m, as.numeric(nMC)))) {
      nMC = 2^m
      warning(sprintf(
        "'nMC' must be a power of 2. Using nMC=%d.",
        nMC
      ))
    }
  }

  if (idxSamples[1] == "default") {
    # Do nothing
  } else {
    if (idxSamples[1] == "final") {
      idxSamples = object$nSamples
    } else if (!is.numeric(idxSamples)) {
      stop("'idxSamples' must be 'default', 'final', or numeric.")
    }
    predictArgs[[object$idxSamplesArg]] = idxSamples
  }
  p = ncol(object$X)
  nMV = object$basisInfo$nMV

  if ((class(object$bmList[[1]])[1] == 'bass') && (is.null(nMC))) {
    class(object) = 'bassBasis'
    names(object)[names(object) == 'bmList'] = 'mod.list'
    names(object)[names(object) == 'basisInfo'] = 'dat'
    object$dat$basis = t(object$dat$basis)

    # BASS::sobolBasis takes a single posterior sample at a time (and names the
    # argument 'mcmc.use' regardless of object$idxSamplesArg), so loop over the
    # requested samples
    if (idxSamples[1] == "default") {
      idxUse = list(NULL)
    } else {
      idxUse = as.list(idxSamples)
    }
    nUse = length(idxUse)

    firstOrder = array(0, dim=c(nUse, p, nMV))
    varTotal = matrix(0, nUse, nMV)
    for (idx in 1:nUse){
      bassArgs = list(object, int.order = 1, ...)
      if (!is.null(idxUse[[idx]])) {
        bassArgs$mcmc.use = idxUse[[idx]]
      }
      out.bass = do.call(BASS::sobolBasis, bassArgs)
      firstOrder[idx,,] = out.bass$S.var[1,1:p,]
      varTotal[idx, ] = out.bass$S.var[1,1,] / out.bass$S[1,1,]
    }

    varTotal = colMeans(varTotal)

    totalOrder = NULL
    if (totalSobol) {
      totalSobol = FALSE
      warning(" 'BASS' has not implemented totalOrder")
    }

    out = list(
      firstOrderSobol = firstOrder,
      totalOrderSobol = totalOrder,
      varTotal = varTotal,
      nMV = nMV,
      p = p
    )
  } else {
    if (is.null(nMC)) {
      nMC = 2^12
    }
    # Generate random samples of parameters according to Saltelli
    # (2010) method.
    baseSequence = rsobol(m = log2(nMC), s = p * 2)
    nMC = nrow(baseSequence)
    A = baseSequence[, 1:p, drop = FALSE]
    B = baseSequence[, (p + 1):(p * 2), drop = FALSE]
    rm("baseSequence")
    X = rbind(A, B)
    for (i in 1:p) {
      AB = A
      AB[, i] = B[, i]
      X = rbind(X, AB)
    }
    AB = X[(2 * nMC + 1):nrow(X), , drop = FALSE]

    saltelliSequence = rbind(A, B, AB)
    rm("A", "B", "AB", "X")

    xmin = apply(object$X, 2, min)
    xrange = apply(object$X, 2, max) - xmin
    # matrix() recycles column-wise, which is correct for a per-input vector
    # (and unlike t(replicate(...)) it stays a matrix when p == 1)
    nSaltelli = nrow(saltelliSequence)
    xmintmp = matrix(xmin, nSaltelli, p, byrow = TRUE)
    xrange = matrix(xrange, nSaltelli, p, byrow = TRUE)
    saltelliSequence = saltelliSequence * xrange + xmintmp

    # evaluate model at those param values
    saltelliMC = do.call(
      predict,
      c(
        list(object, saltelliSequence),
        predictArgs
      )
    )
    rm("saltelliSequence")

    # Normalize to c(nSamplesMC, N, nMV) so each posterior sample gets its own
    # set of Sobol' indices
    if (ndims(saltelliMC) == 3) {
      nSamplesMC = dim(saltelliMC)[1]
    } else {
      nSamplesMC = 1
      saltelliMC = array(saltelliMC, dim = c(1, dim(saltelliMC)))
    }

    # If predict.bayesModel ignored idxSamplesArg it returned every draw, so
    # honor the request here rather than computing indices for all of them
    if (is.numeric(idxSamples) && (nSamplesMC != length(idxSamples)) &&
        all(idxSamples <= nSamplesMC)) {
      saltelliMC = saltelliMC[idxSamples, , , drop = FALSE]
      nSamplesMC = length(idxSamples)
    }

    basisType = object$basisInfo$basisType

    firstOrder = NULL
    totalOrder = NULL
    varTotal = NULL

    for (idxMC in 1:nSamplesMC) {
      mcMat = matrix(saltelliMC[idxMC, , ],
                     nrow = dim(saltelliMC)[2],
                     ncol = dim(saltelliMC)[3])
      meanS = colMeans(mcMat)
      mcMat = mcMat - matrix(meanS, nrow(mcMat), ncol(mcMat), byrow = TRUE)

      if ((basisType == 'jfpca') || (basisType == 'jfpcah')) {
        mcMat = .sobolElasticWarp(mcMat, object, basisType, nMV)
      }
      nMVout = ncol(mcMat)

      if (is.null(firstOrder)) {
        firstOrder = array(0, dim = c(nSamplesMC, p, nMVout))
        if (totalSobol) {
          totalOrder = array(0, dim = c(nSamplesMC, p, nMVout))
        }
        varTotal = matrix(0, nSamplesMC, nMVout)
      }

      modA = mcMat[1:nMC, , drop = FALSE]
      modB = mcMat[(nMC + 1):(nMC * 2), , drop = FALSE]

      varTotal[idxMC, ] = apply(mcMat, 2, var)

      for (j in 1:p) {
        modAB = mcMat[((2 + (j - 1)) * nMC + 1):((3 + (j - 1)) * nMC), , drop = FALSE]
        firstOrder[idxMC, j, ] = colMeans(modB * (modAB - modA))
        if (totalSobol) {
          totalOrder[idxMC, j, ] = 0.5 * colMeans((modA - modAB)^2)
        }
      }
      rm(mcMat, modA, modB)
    }
    nMV = dim(firstOrder)[3]
    rm("saltelliMC")

    # Truncate at zero (machine precision can give negative values)
    firstOrder[firstOrder < 0] <- 0

    varTotal = colMeans(varTotal)

    # Ensure total variance >= sum of first order terms
    firstOrderSum = apply(apply(firstOrder, 2:3, mean), 2, sum)
    varTotal = apply(rbind(varTotal, firstOrderSum), 2, max)

    out = list(
      firstOrderSobol = firstOrder,
      totalOrderSobol = totalOrder,
      varTotal = varTotal,
      nMV = nMV,
      p = p
    )
  }

  out = structure(out, class = 'sobol')

  if (showPlot) {
    plot(out, totalSobol = totalSobol)
  }

  return(out)
}

# Map jfpca/jfpcah posterior draws back to the aligned function space, so that
# Sobol' indices are computed on the functions rather than on the joint
# (function, warping) representation.
.sobolElasticWarp = function(mcMat, object, basisType, nMV) {
  C = object$basisInfo$basisConstruct$C
  N = nrow(mcMat)

  if (nMV %% 2 == 0) {
    M = floor(nMV / 2)
    idxWarp = (M + 1):nMV
  } else {
    M = floor((nMV - 1) / 2)
    idxWarp = (M + 2):nMV
  }
  time = seq(0, 1, length.out = M)
  postSamples = array(0, dim = c(N, M))

  if (basisType == 'jfpca') {
    gamtmp = fdasrvf::v_to_gam(t(mcMat[, idxWarp, drop = FALSE] / C))
  } else {
    gamtmp = fdasrvf::h_to_gam(t(mcMat[, idxWarp, drop = FALSE] / C))
  }

  if (nMV %% 2 == 0) {
    for (jj in 1:N) {
      postSamples[jj, ] = fdasrvf::warp_f_gamma(mcMat[jj, 1:M], time, gamtmp[, jj])
    }
  } else {
    mididx = object$basisInfo$basisConstruct$id
    for (jj in 1:N) {
      ftmp = mcMat[jj, 1:(M + 1)]
      ftmp = cumtrapzmid(time,
                         ftmp[1:M] * abs(ftmp[1:M]),
                         sign(ftmp[M + 1]) * (ftmp[M + 1]^2),
                         mididx)
      postSamples[jj, ] = fdasrvf::warp_f_gamma(ftmp, time, gamtmp[, jj])
    }
  }

  return(postSamples)
}

#' @title Plot Sobol Decomposition
#'
#' @description Given an object of class "mvSobol" from the mvSobol() function
#' @param x An object of class "mvSobol" containing the Sobol Indices
#' @param totalSobol A boolean to plot the total sobol (default = `TRUE`)
#' @param labels A character vector of length <= 8 containing the names of the parameters
#' @param idxMV A vector defining the time points
#' @param xscale string whether to plot on a "linear" scale or "log"
#' @param xlabel string for the xlabel
#' @param yOverlay A boolean if to overlay the experimental function (default = `NULL`)
#' @param yOverlayLabel A string defining the overlay label
#' @param waterfall bool whether to plot sobol as a waterfall (functional pie-chart) default (`FALSE`)
#' @param file An optional location at which the plots will be saved. If NULL, no file is saved.
#' @param title An optional title to be printed at the top of the traceplots.
#' @param ... additional plot arguments
#' @return no return value
#' @seealso See \link{mvBayes}
#' @export
#' @import graphics
#'
plot.sobol = function(x,
                      totalSobol = TRUE,
                      labels = NULL,
                      idxMV = NULL,
                      xscale = "linear",
                      xlabel = "Multivariate Index",
                      yOverlay = NULL,
                      yOverlayLabel = "Overlay",
                      waterfall = FALSE,
                      title = NULL,
                      file = NULL,
                      ...) {
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

  p = x$p

  # Decided up front, since it determines the number of panels below
  if (totalSobol && is.null(x$totalOrderSobol)) {
    warning("'x' does not have totalOrderSobol")
    totalSobol = FALSE
  }

  if (is.null(idxMV)) {
    idxMV = 1:x$nMV
  }

  if (is.null(labels)) {
    labels = paste0("X", 1:p)
  }
  labels = c(labels, "Higher-Order")

  lty = c(rep(1:4, ceiling(p/4))[1:p], 1)
  if (waterfall) {
    ltyLegend = rep(1, p+1)
  } else {
    ltyLegend = lty
  }

  cmap = c(cmap[0:(p - 1) %% length(cmap) + 1], "grey")
  cmapLegend = cmap

  lwdLegend <- rep(2, p+1)

  if (!is.null(yOverlay)) {
    labels = c(labels, yOverlayLabel)
    ltyLegend = c(ltyLegend, 2)
    cmapLegend = c(cmapLegend, "black")
    lwdLegend = c(lwdLegend, 1)
  }

  # matrix() keeps p x nMV even when p == 1
  firstOrder = matrix(apply(x$firstOrderSobol, c(2, 3), mean),
                      nrow = p,
                      ncol = x$nMV)

  firstOrderRel = t(t(firstOrder) / x$varTotal)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(
    mfrow = c(1, 2 + totalSobol),
    mar = c(5, 5, 1, 1),
    oma=c(0, 0, 2, 0)
  )

  if (waterfall) {
    meanX = rbind(
      firstOrderRel,
      1.0 - apply(firstOrderRel, 2, sum)
    )
    # t(apply(..., cumsum)) collapses to a vector when there is only one row
    sens = matrix(apply(meanX, 2, cumsum), nrow = p + 1, ncol = x$nMV)
    sens = t(sens)

    plot(
      idxMV,
      rep(0, x$nMV),
      type = "l",
      col = cmap[1],
      ylim = c(0, 1),
      xlim = c(min(idxMV), max(idxMV)),
      log = ifelse(xscale == 'log', 'x', ''),
      xlab = xlabel,
      ylab = "Relative First-Order Sobol' Index"
    )
    polygon(c(idxMV, rev(idxMV)), c(rep(0, x$nMV), rev(sens[, 1])), col = cmap[1])
    for (j in seq_len(p - 1) + 1) {
      polygon(c(idxMV, rev(idxMV)), c(sens[, j - 1], rev(sens[, j])), col = cmap[j])
    }
    lines(idxMV, sens[, p + 1], type = "l", col = "grey")
    polygon(c(idxMV, rev(idxMV)), c(sens[, p], rev(sens[, p + 1])), col = "grey")

  } else {
    plot(
      idxMV,
      firstOrderRel[1, ],
      type = "l",
      lwd = 2,
      lty = lty[1],
      col = cmap[1],
      xlab = xlabel,
      ylab = "Relative First-Order Sobol' Index",
      log = ifelse(xscale == 'log', 'x', ''),
      ylim = c(0, 1),
      xlim = c(min(idxMV), max(idxMV))
    )
    for (j in seq_len(p - 1) + 1) {
      lines(
        idxMV,
        firstOrderRel[j, ],
        lwd = 2,
        lty = lty[j],
        col = cmap[j]
      )
    }
    lines(idxMV,
          1.0 - apply(firstOrderRel, 2, sum),
          lwd = 2,
          col = "grey")
  }

  if (!is.null(yOverlay)) {
    par(new=TRUE)
    plot(
      idxMV, yOverlay,
      type = 'l', lty = 2,
      col = grDevices::rgb(0, 0, 0, 0.7),
      axes = FALSE, xlab = "", ylab = ""
    )
    axis(side = 4, at = pretty(range(yOverlay)))
  }

  if (waterfall) {
    sensVar = t(rbind(matrix(apply(firstOrder, 2, cumsum), nrow = p, ncol = x$nMV),
                      x$varTotal))

    plot(
      idxMV,
      rep(0, x$nMV),
      type = "l",
      col = cmap[1],
      ylim = c(0, max(sensVar)*1.1),
      xlim = c(min(idxMV), max(idxMV)),
      log = ifelse(xscale == 'log', 'x', ''),
      xlab = xlabel,
      ylab = "First-Order Sobol' Index"
    )
    polygon(c(idxMV, rev(idxMV)), c(rep(0, x$nMV), rev(sensVar[, 1])), col = cmap[1])
    for (j in seq_len(p - 1) + 1) {
      polygon(c(idxMV, rev(idxMV)), c(sensVar[, j - 1], rev(sensVar[, j])), col = cmap[j])
    }
    lines(idxMV, sensVar[, p + 1], type = "l", col = "grey")
    polygon(c(idxMV, rev(idxMV)), c(sensVar[, p], rev(sensVar[, p + 1])), col = "grey")
  } else {
    plot(
      idxMV,
      firstOrder[1, ],
      type = "l",
      lwd = 2,
      lty = lty[1],
      col = cmap[1],
      xlab = xlabel,
      ylab = "First-Order Sobol' Index",
      log = ifelse(xscale == 'log', 'x', ''),
      ylim = c(0, max(firstOrder)*1.1),
      xlim = c(min(idxMV), max(idxMV))
    )
    for (j in seq_len(p - 1) + 1) {
      lines(
        idxMV,
        firstOrder[j, ],
        lwd = 2,
        lty = lty[j],
        col = cmap[j]
      )
    }
    lines(
      idxMV,
      x$varTotal - apply(firstOrder, 2, sum),
      lwd = 2,
      col = "grey"
    )
  }

  if (!is.null(yOverlay)) {
    par(new=TRUE)
    plot(
      idxMV, yOverlay,
      type = 'l', lty = 2,
      col = grDevices::rgb(0, 0, 0, 0.7),
      axes = FALSE, xlab = "", ylab = ""
    )
    axis(side = 4, at = pretty(range(yOverlay)))
  }

  legend(
    "topleft",
    legend = labels,
    col = cmapLegend,
    lwd = lwdLegend,
    lty = ltyLegend,
    xpd = NA,
    cex = 0.65,
    y.intersp = 0.75,
    text.width=0.25*diff(range(idxMV)),
    bg = grDevices::rgb(1, 1, 1, 0.2)
  )

  # total sobol plot
  if (totalSobol) {

    totalOrder = matrix(apply(x$totalOrderSobol, c(2, 3), mean),
                        nrow = p,
                        ncol = x$nMV)

    plot(
      idxMV,
      totalOrder[1, ],
      type = "l",
      lwd = 2,
      lty = lty[1],
      col = cmap[1],
      xlab = xlabel,
      ylab = "Total-Order Sobol' Index",
      log = ifelse(xscale == 'log', 'x', ''),
      ylim = c(0, max(totalOrder) * 1.05),
      xlim = c(min(idxMV), max(idxMV))
    )
    for (j in seq_len(p - 1) + 1) {
      lines(
        idxMV,
        totalOrder[j, ],
        lwd = 2,
        lty = lty[j],
        col = cmap[j]
      )
    }
  }

  if (!is.null(yOverlay)) {
    par(new=TRUE)
    plot(
      idxMV, yOverlay,
      type = 'l', lty = 2,
      col = grDevices::rgb(0, 0, 0, 0.7),
      axes = FALSE, xlab = "", ylab = ""
    )
    axis(side = 4, at = pretty(range(yOverlay)))
  }

  if (!is.null(title)) {
    mtext(title, outer = TRUE, font = 2) # Add title
  }

  # Save or display plot
  if (!is.null(file)) {
    grDevices::dev.copy(grDevices::png, file, ...)
    grDevices::dev.off()
  }
}
