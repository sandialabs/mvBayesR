customBasisConstruct <- function(customBasis, Ystandard) {
  out = structure(list(), class = "customBasisConstruct")
  nBasis <- nrow(customBasis)
  nMV <- ncol(customBasis)

  qr_decomp <- qr(t(customBasis))
  Q <- qr.Q(qr_decomp)

  out$basis <- t(Q)
  out$coefs <- Ystandard %*% t(out$basis)
  out$varExplained <- apply(out$coefs, 2, var)

  if (nBasis < nMV) {
    totalVar <- sum(apply(Ystandard, 2, var))
    out$varExplained <- c(out$varExplained, totalVar - sum(out$varExplained))
  }

  return(out)
}
