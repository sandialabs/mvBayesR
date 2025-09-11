ndims <- function(x) {
  return(length(dim(x)))
}

normalize_column <- function(x) {
  magnitude <- sqrt(sum(x^2))
  if (magnitude == 0) {
    return(x)  # Handle zero-magnitude columns
  }
  x / magnitude
}

repmat <- function(X, m, n) {
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  if (is.null(mx)) {
    mx = 1
    nx = length(X)
    mat = matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
  } else {
    nx = dim(X)[2]
    mat = matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
  }

  return(mat)
}

cumtrapz <- function(x, y, dims = 1) {
  if ((dims - 1) > 0) {
    perm = c(dims:max(ndims(y), dims), 1:(dims - 1))
  } else {
    perm = c(dims:max(ndims(y), dims))
  }

  if (ndims(y) == 0) {
    n = 1
    m = length(y)
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
    n = ncol(y)
  }

  if (n == 1) {
    dt = diff(x) / 2.0
    z = c(0, cumsum(dt * (y[1:(m - 1)] + y[2:m])))
    dim(z) = c(m, 1)
  } else {
    tmp = diff(x)
    dim(tmp) = c(m - 1, 1)
    dt = repmat(tmp / 2.0, 1, n)
    z = rbind(rep(0, n), apply(dt * (y[1:(m - 1), ] + y[2:m, ]), 2, cumsum))
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    z = aperm(z, perm2)
  }

  return(z)
}

cumtrapzmid <- function(x, y, c, mid) {
  a = length(x)

  # case < mid
  fn = rep(0, a)
  tmpx = x[seq(mid - 1, 1, -1)]
  tmpy = y[seq(mid - 1, 1, -1)]
  tmp = c + cumtrapz(tmpx, tmpy)
  fn[1:(mid - 1)] = rev(tmp)

  # case >= mid
  fn[mid:a] = c + cumtrapz(x[mid:a], y[mid:a])

  return(fn)
}
