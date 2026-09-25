evaluateExtdepthCoverage <- function(preds, Ytest, alpha = 0.05) {
  # preds: nSamples x nTest x nMV array
  # Ytest: nTest x nMV matrix
  #
  # extdepth::edepth_set() expects fmat to be nMV x nSamples,
  # with each column representing one function.
  
  if (!is.array(preds) || length(dim(preds)) != 3L) {
    stop("preds must be a 3D array with dimensions nSamples x nTest x nMV.")
  }
  
  if (!is.matrix(Ytest)) {
    stop("Ytest must be a matrix with dimensions nTest x nMV.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in [0, 1).")
  }
  
  nSamples <- dim(preds)[1]
  nTest <- dim(preds)[2]
  nMV <- dim(preds)[3]
  
  stopifnot(nrow(Ytest) == nTest)
  stopifnot(ncol(Ytest) == nMV)
  
  covered <- logical(nTest)
  intervalWidths <- numeric(nTest)
  
  lowerMat <- matrix(NA_real_, nrow = nTest, ncol = nMV)
  upperMat <- matrix(NA_real_, nrow = nTest, ncol = nMV)
  
  useMatrixStats <- requireNamespace("matrixStats", quietly = TRUE)
  
  for (i in seq_len(nTest)) {
    # predsIdx: nSamples x nMV
    predsIdx <- matrix(
      preds[, i, , drop = FALSE],
      nrow = nSamples,
      ncol = nMV
    )
    
    # extdepth::edepth_set() expects nMV x nSamples.
    fmat <- t(predsIdx)
    
    eDepths <- extdepth::edepth_set(fmat)
    
    keep <- eDepths > alpha
    
    if (!any(keep)) {
      stop(
        paste0(
          "No posterior predictive samples retained in central region ",
          "for test function ", i, ". Consider using a smaller alpha."
        )
      )
    }
    
    # Since predsIdx is nSamples x nMV and keep indexes posterior samples,
    # this is equivalent to applying central_region() to fmat, but avoids
    # an extra subsetting/transposition path.
    retainedPreds <- predsIdx[keep, , drop = FALSE]
    
    if (useMatrixStats) {
      lower <- matrixStats::colMins(retainedPreds)
      upper <- matrixStats::colMaxs(retainedPreds)
    } else {
      lower <- apply(retainedPreds, 2L, min)
      upper <- apply(retainedPreds, 2L, max)
    }
    
    lowerMat[i, ] <- lower
    upperMat[i, ] <- upper
    
    insideBand <- Ytest[i, ] >= lower & Ytest[i, ] <= upper
    
    covered[i] <- all(insideBand)
    intervalWidths[i] <- mean(upper - lower)
  }
  
  list(
    coverageTarget = 1 - alpha,
    simultaneousCoverage = mean(covered),
    intervalWidths = intervalWidths,
    lower = lowerMat,
    upper = upperMat
  )
}
