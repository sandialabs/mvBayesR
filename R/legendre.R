basisLegendre <- function(M, N) {
  coeff = matrix(0, N + 1, N + 1)
  coeff[c(1, N + 3)] = 1
  for (ii in 3:(N + 1)) {
    coeff[ii, ] = (2 - 1 / (ii - 1)) * coeff[ii - 1, c(N + 1, 1:(N + 1 - 1))] -
      (1 - 1 / (ii - 1)) * coeff[ii - 2, ]
  }
  x = seq(-1, 1, length.out = M)
  basismat = coeff %*% apply(cbind(rep(1, M), replicate(N, x)), 1, cumprod)

  return(basismat)
}
