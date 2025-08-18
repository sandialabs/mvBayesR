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
#' @keywords Sobol' decomposition of variance, BASS
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
  if (!methods::is(object, 'mvBayes')) {
    stop("'object' must be of class 'mvBayes'")
  }

  predictArgs = list(...)

  if (object$basisInfo$basisType == "pns" && (is.null(nMC))) {
    nMC = 2^12
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
  nSamples = length(idxSamples)   # stub out to use more indices
  p = ncol(object$X)
  nMV = object$basisInfo$nMV

  if ((class(object$bmList[[1]])[1] == 'bass') && (is.null(nMC))) {
    class(object) = 'bassBasis'
    names(object)[names(object) == 'bmList'] = 'mod.list'
    names(object)[names(object) == 'basisInfo'] = 'dat'

    out.bass = BASS::sobolBasis(object, int.order=1, ...)

    firstOrder = array(0, dim = c(nSamples, p, nMV))
    varTotal = matrix(0, nSamples, nMV)

    firstOrder = array(0, dim = c(nSamples, p, nMV))
    varTotal = matrix(0, nSamples, nMV)

    totalOrder = NULL
    if (totalSobol) {
      warning(" 'BASS' has not implemented totalOrder")
    }

    totalOrder = NULL
    if (totalSobol) {
      warning(" 'BASS' has not implemented totalOrder")
    }

    out = list(
      firstOrderSobol = firstOrder,
      totalOrderSobol = totalOrder,
      varTotal = varTotal,
      # sobolProp = firstOrder / aperm(replicate(p, varTotal), c(1, 3, 2)),
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
    xmintmp = t(replicate(nrow(saltelliSequence), xmin))
    xrange = apply(object$X, 2, max) - xmin
    xrange = t(replicate(nrow(saltelliSequence), xrange))
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

    if (ndims(saltelliMC) == 3){
      meanS = colMeans(saltelliMC[1, , ])
      meanS = t(replicate(nrow(saltelliMC[1, , ]), meanS))
      saltelliMC = saltelliMC[1, , ] - meanS
    } else {
      meanS = colMeans(saltelliMC)
      meanS = t(replicate(nrow(saltelliMC), meanS))
      saltelliMC = saltelliMC - meanS
    }

    basisType = object$basisInfo$basisType

    if ((basisType == 'jfpca') || (basisType == 'jfpcah')) {
      C = object$basisInfo$basisConstruct$C
      if (nMV %% 2 == 0) {
        M = floor(nMV / 2)
        N = nrow(saltelliMC)
        time = seq(0, 1, length.out = M)
        post_samples = array(0, dim = c(N, M))
        if (basisType == 'jfpca') {
          gamtmp = fdasrvf::v_to_gam(t((saltelliMC[, (M + 1):nMV] / C)))
        } else if (basisType == 'jfpcah') {
          gamtmp = fdasrvf::h_to_gam(t((saltelliMC[, (M + 1):nMV] / C)))
        }
        for (jj in 1:N) {
          ftmp = saltelliMC[jj, 1:M]
          post_samples[jj, ] = fdasrvf::warp_f_gamma(ftmp, time, gamtmp[, jj])
        }

        saltelliMC = post_samples
        nMV = M
        rm("post_samples")

      } else {
        M = floor((nMV - 1) / 2)
        N = nrow(saltelliMC)
        time = seq(0, 1, length.out = M)
        post_samples = array(0, dim = c(N, M))
        mididx = object$basisInfo$basisConstruct$id
        if (basisType == 'jfpca') {
          gamtmp = fdasrvf::v_to_gam(t((saltelliMC[, (M + 2):nMV] / C)))
        } else if (basisType == 'jfpcah') {
          gamtmp = fdasrvf::h_to_gam(t((saltelliMC[, (M + 2):nMV] / C)))
        }
        for (jj in 1:N) {
          ftmp = saltelliMC[jj, 1:(M + 1)]
          ftmp = cumtrapzmid(time,
                             ftmp[1:M] * abs(ftmp[1:M]),
                             sign(ftmp[M + 1]) * (ftmp[M + 1]^2),
                             mididx)
          post_samples[jj, ] = fdasrvf::warp_f_gamma(ftmp, time, gamtmp[, jj])
        }

        saltelliMC = post_samples
        nMV = M
        rm("post_samples")

      }
    }

    modA = saltelliMC[1:nMC, ]
    modB = saltelliMC[(nMC + 1):(nMC * 2), ]
    modAB = list()
    for (i in 1:p) {
      modAB[[i]] = saltelliMC[((2 + (i - 1)) * nMC + 1):((3 + (i - 1)) * nMC), ]
    }

    varTotal = apply(saltelliMC, 2, var)
    rm("saltelliMC")

    firstOrder = array(0, dim = c(nSamples, p, nMV))
    totalOrder = NULL
    if (totalSobol) {
      totalOrder = array(0, dim = c(nSamples, p, nMV))
    }
    for (j in 1:p) {
      firstOrder[, j, ] = apply(modB * (modAB[[j]] - modA), 2, mean)
      if (totalSobol) {
        totalOrder[, j, ] = 0.5 * apply((modA - modAB[[j]])^2, 2, mean)
      }
    }

    tmp = t(replicate(p, c(varTotal)))
    out = list(
      firstOrderSobol = firstOrder,
      totalOrderSobol = totalOrder,
      varTotal = varTotal,
      # sobolProp = firstOrder / aperm(replicate(nSamples, tmp), c(3, 1, 2)),
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

#' @title Plot Sobol Decomposition
#'
#' @description Given an object of class "mvSobol" from the mvSobol() function
#' @param object An object of class "mvSobol" containing the Sobol Indices
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
#' @keywords multivariate Bayesian regression modeling, functional data analysis
#' @seealso See \link{mvBayes}
#' @export
#' @import graphics
#'
plot.sobol = function(object,
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

  p = object$p

  if (is.null(idxMV)) {
    idxMV = 1:object$nMV
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

  firstOrder = apply(object$firstOrderSobol, c(2, 3), mean)
  firstOrderRel = t(t(firstOrder) / object$varTotal)

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
    sens = t(apply(meanX, 2, cumsum))

    plot(
      idxMV,
      rep(0, object$nMV),
      type = "l",
      col = cmap[1],
      ylim = c(0, 1),
      xlim = c(min(idxMV), max(idxMV)),
      log = ifelse(xscale == 'log', 'x', ''),
      xlab = xlabel,
      ylab = "Relative First-Order Sobol' Index"
    )
    polygon(c(idxMV, rev(idxMV)), c(rep(0, object$nMV), rev(sens[, 1])), col = cmap[1])
    for (j in 2:p) {
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
    for (j in 2:p) {
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
    sensVar = t(rbind(apply(firstOrder, 2, cumsum), object$varTotal))

    plot(
      idxMV,
      rep(0, object$nMV),
      type = "l",
      col = cmap[1],
      ylim = c(0, max(sensVar) + 3),
      xlim = c(min(idxMV), max(idxMV)),
      log = ifelse(xscale == 'log', 'x', ''),
      xlab = xlabel,
      ylab = "First-Order Sobol' Index"
    )
    polygon(c(idxMV, rev(idxMV)), c(rep(0, object$nMV), rev(sensVar[, 1])), col = cmap[1])
    for (j in 2:p) {
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
      ylim = c(0, max(firstOrder) * 1.05),
      xlim = c(min(idxMV), max(idxMV))
    )
    for (j in 2:p) {
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
      object$varTotal - apply(firstOrder, 2, sum),
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
    if (is.null(object$totalOrderSobol)) {
      warning(
        "'object' does not have totalOrderSobol"
      )
    }
    totalOrder = apply(object$totalOrderSobol, c(2, 3), mean)

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
    for (j in 2:p) {
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
