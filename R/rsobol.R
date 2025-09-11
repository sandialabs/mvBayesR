# This code computes scrambled Sobol' points.
# It has been used for academic research purposes.
# It comes with no guarantees.  It might be
# much slower than other implementations.
# The primary purpose is to provide source
# for `nested uniform scrambling' (often called
# `Owen scrambling') of Sobol' points.
#
# All functionality is through a main function rsobol
#   the default is nested uniform scrambling
#   set rand = FALSE to get unrandomized Sobol'
#   set type = "mato" gets Matousek's linear scramble
#
# Basic use:
#    place files rsobol.R, fiftysobol.col, sobol_Cs.col in a directory
#    From within R in that directory run:
#      source("rsobol.R")
#      x <- rsobol(m=15,s=10,seed=1)
#   Then x will have n=2^15 scrambled Sobol' points in s=10 dimensions.
#   Repeating with the same seed gives back the same points
#   unless the underyling R pseudo-random number structure has changed.
#   There will be n rows and s columns.  For s < 50 <= 20201 use
#   x <- rsobol(fn="sobol_Cs.col",m,s,seed)
#
# Some helper functions are available with names starting with .rsobol
# That way they are available to users who want them but they don't clutter
# the namespace.  Prepending .rsobol in a few places makes some of the
# code a bit less readable as part of that tradeoff.

# This code uses generating matrices from Dirk Nuyens' magic point shop
# described in this article:
#   F.Y. Kuo & D. Nuyens. Application of quasi-Monte Carlo methods
#   to elliptic PDEs with random diffusion coefficients - a survey of
#   analysis and implementation, Foundations of Computational Mathematics,
#   16(6):1631-1696, 2016.
#
# The underlying direction numbers come from
#   Joe, Stephen, and Frances Y. Kuo.
#   "Constructing Sobol sequences with better two-dimensional projections."
#   SIAM Journal on Scientific Computing 30, no. 5 (2008): 2635-2654.
# Those papers should be cited in publications based on this code.
#
# fiftysobol.col provides generating matrices for up to 50 dimensions.
# sobol_Cs.col handles up to 21021 dimensions.
# The code reads one of them from disk at each use.
#
# The first 2^m points of rsobol(fn,m+1,s) are not the same as
# the output from rsobol(fn,m,s).
#

rsobol = function(fn = "fiftysobol.col",
                  m = 10,
                  s = 5,
                  rand = TRUE,
                  type = "nestu",
                  M = 32,
                  seed = 20171215) {
  # Scramble Sobol' points.  Uses nested uniform points by default.
  # fn   = file of matrix values from Dirk Nuyens' magic point shop
  # m    = log2( n ), i.e., it makes 2^n points
  # s    = dimension of points
  # rand = TRUE for randomized Sobol' (the default), FALSE for original Sobol' (use with care)
  # type = nestu (for nested uniform scramble) or mato (for Matousek's linear scramble)
  # M    = number of bits to produce; must match number of columns in fn
  # seed = one random seed for the whole matrix
  #
  # returns matrix in [0,1]^{n * s }of scrambled Sobol' points
  #

  if (m > M)
    stop(paste(
      sep = "",
      "Cannot deliver 2**",
      m,
      " points. Note: parameter m is log2( n )."
    ))
  # This is to catch users who ask for m=1000 or m=1024 or more thinking that m is the sample size
  # Also, plain Sobol' points will repeat with a period of 2^M
  # Randomizing beyond that point will be just like replicates of 2^M points (at best)
  # and it would be better to just change the seed and get replicates.

  set.seed(seed)

  if (!rand)
    # User wants unrandomized Sobol'
    return(.rsobol.sobolpts(fn, m = m, s = s, M = M))

  if (type == "mato")
    # User wants Matousek's random linear scrambles
    return(.rsobol.rsobolptsmato(
      fn,
      m = m,
      s = s,
      M = M,
      seed = seed
    ))

  if (type != "nestu")
    stop(paste("Unknown scramble type: ", type, ". Should be nestu or mato"))

  if (m == 0)
    return(matrix(runif(s), nrow = 1))  # Avoids a corner case in workhorse function

  n = 2^m
  ans = matrix(0, n, s)
  newbits = .rsobol.rsobolbits(fn, m, s, M, seed)
  for (i in 1:n)
    for (j in 1:s)
      ans[i, j] = .rsobol.bits2unif(newbits[i, j, ])
  ans
}


# The next functions are helper functions for rsobol.
# Their names start with .rsobol.
# Then come test functions.

.rsobol.int2bits = function(x, M = 32) {
  # Convert an integer x into M bits
  # For large x (or small M) you get bits of x modulo 2^M
  # as.vector() is there to force answer to be a vector.
  # This does just one integer (not a vector of them).

  ans = rep(0, M)
  for (j in 1:M) {
    ans[j] = as.vector(x %% 2)
    x = (x - ans[j]) / 2
  }
  ans
}

.rsobol.sobolbits = function(fn = "fiftysobol.col",
                             m = 8,
                             s = 20,
                             M = 32) {
  # Get array of Sobol' bits.
  # Calls .rsobol.int2bits and .rsobol.sobomats
  # M must be the number of columns of data in the file named fn

  n   = 2^m
  ans = array(0, c(n, s, M))             # n obs x s dimensions x M bits
  a   = .rsobol.sobomats(fn, m, s, M)

  bits = rep(0, M)
  for (i in 1:n) {
    bitsi = .rsobol.int2bits(i - 1, M)   # bits of integer i-1 used for observation i
    for (j in 1:s) {
      bitsj = a[j, , ] %*% bitsi        # a[j,,] is the Sobol' matrix for dimension j
      ans[i, j, ] = bitsj %% 2
    }
  }

  ans
}

.rsobol.sobomats = function(fn = "fiftysobol.col",
                            m = 8,
                            s = 20,
                            M = 32) {
  # main workhorse function for sobolbits
  # Return an array of Sobol' matrices.
  # Not tuned for efficiency in time or space.
  # For instance, the matrices have many zeros.
  #
  data <- system.file("extdata", fn, package = "mvBayes")
  col = utils::read.table(data)
  if (s > nrow(col))
    stop(paste(
      sep = "",
      "Not enough colunns in file fn = ",
      fn,
      ". There should be at least s = ",
      s,
      "."
    ))
  a = array(0, c(s, M, M))
  for (j in 1:s) {
    for (k in 1:M)
      a[j, , k] = .rsobol.int2bits(col[j, k], M)
  }
  a
}

.rsobol.bits2unif = function(bits) {
  # Turn sequence of bits into a point in [0,1)
  # First bits are highest order
  ans = 0
  for (j in length(bits):1) {
    ans = (ans + bits[j]) / 2
  }
  ans
}

.rsobol.getpermset2 = function(J) {
  # Get 2^(j-1) random binary permutations for j=1 ... J
  # J will ordinarily be m when there are n=2^m points
  #
  # A nuisance is that m=0 gives J=0, requiring a list of length 0
  # that the for loop doesn't do as desired.
  # The caller will handle that corner case a different way.
  #
  # Caller has set the seed
  ans = list()
  for (j in 1:J) {
    ans[[j]] = as.integer(runif(2^(j - 1)) > 1 / 2)
  }
  ans
}

.rsobol.bits2int = function(b, M = 32) {
  # Convert bits b into integers.
  # Inverse of int2bits
  # This is vectorized: each row of the matrix b is a vector of bits
  #
  if (is.vector(b))
    b = matrix(b, nrow = 1)

  n = nrow(b)
  p = ncol(b)
  ans = rep(0, n)

  for (j in p:1)
    ans = ans * 2 + b[, j]
  ans
}


# Workhorse function for rsobol with nested uniform scramble
.rsobol.rsobolbits = function(fn = "fiftysobol.col",
                              m = 10,
                              s = 5,
                              M = 32,
                              seed = 20171215) {
  # Scramble Sobol' bits; nested uniform.
  # uses .rsobol.sobolbits and .rsobol.getpermset2
  #
  set.seed(seed)

  if (m < 1)
    stop("We need m >= 1") # m=0 causes awkward corner case below.  Caller handles that case specially.

  thebits = .rsobol.sobolbits(fn, m, s, M)         # The unrandomized Sobol' bits
  newbits = thebits                            # Initialize the randomized bits

  n = 2^m
  for (j in 1:s) {
    theperms = .rsobol.getpermset2(m)          # Permutations to apply to bits 1:m
    for (k in 1:m) {
      # Here is where we want m > 0 so the loop works ok
      if (k > 1) {
        bitmat = thebits[, j, ]                  # slice on dim j to get a matrix
        bitmat = bitmat[, 1:(k - 1), drop = FALSE]   # drop=FALSE will keep this as a matrix, even for k=2
        indices = .rsobol.bits2int(bitmat)     # index of which perms to use at bit k for each i
      } else{
        indices = rep(0, n)                     # same permutation for all observations i
      }
      newbits[, j, k] = (thebits[, j, k] + theperms[[k]][1 + indices]) %% 2   # permutation by adding a bit modulo 2
    }
  }
  if (M > m)
    # Paste in random entries for bits after m'th one
    newbits[, , (m + 1):M] = runif(n * s * (M - m)) > 0.5
  newbits
}

.rsobol.getmatousek2 = function(J) {
  # Genereate the Matousek linear scramble in base 2 for one of the s components
  # We need a J x J bit matrix M and a length J bit vector C
  #
  M = diag(J) # Identity
  C = runif(J) > 0.5
  for (i in 2:J)
    for (j in 1:(i - 1))
      M[i, j] = runif(1) > 0.5
  list(M = M, C = C^2)             # squaring turns boolean into binary
}


.rsobol.rsobolbitsmato = function(fn = "fiftysobol.col",
                                  m = 10,
                                  s = 5,
                                  M = 32,
                                  seed = 20171215) {
  # workhorse function for randsobolbitsmato
  # Scramble Sobol' bits using Matousek's linear scramble
  # Uses .rsobol.sobolbits and .rsobol.getmatousek2

  if (m < 1)
    stop("Need m >= 1") # m=0 causes awkward corner case below; caller avoids doing this
  set.seed(seed)

  thebits = .rsobol.sobolbits(fn, m, s, M)
  newbits = thebits

  n = 2^m
  for (j in 1:s) {
    themato = .rsobol.getmatousek2(m)
    for (k in 1:m) {
      # not correct for m=1
      if (k > 1) {
        bitmat = thebits[, j, ]                  # slice on dim j to matrix
        bitmat = bitmat[, 1:(k - 1), drop = FALSE]   # keep as matrix, even for k=1
        newbits[, j, k] = (thebits[, j, k] + bitmat %*% themato$M[k, 1:(k -
                                                                          1)]) %% 2
      }
      newbits[, j, k] = (newbits[, j, k] + themato$C[k]) %% 2
    }
  }
  if (M > m)
    newbits[, , (m + 1):M] = runif(n * s * (M - m)) > 0.5
  newbits
}

.rsobol.rsobolptsmato = function(fn = "fiftysobol.col",
                                 m = 10,
                                 s = 5,
                                 M = 32,
                                 seed = 20171215) {
  # Scramble Sobol' points; Matousek's nested linear scrambling
  set.seed(seed)
  if (m == 0)
    return(matrix(runif(s), nrow = 1))  # Avoids a corner case in workhorse function

  n = 2^m
  ans = matrix(0, n, s)
  newbits = .rsobol.rsobolbitsmato(fn, m, s, M, seed)
  for (i in 1:n)
    for (j in 1:s)
      ans[i, j] = .rsobol.bits2unif(newbits[i, j, ])
  ans
}

.rsobol.sobolpts = function(fn = "fiftysobol.col",
                            m = 8,
                            s = 2,
                            M = 32) {
  # Plain Sobol' points
  # Caller gets them as an option to rsobol
  #
  # Warning: these points are not randomized!

  data <- system.file("data", fn, package = "mvBayes")
  col = utils::read.table(data)
  if (s > nrow(col))
    stop(paste(
      sep = "",
      "Not enough colunns in file fn = ",
      fn,
      ". There should be at least s = ",
      s,
      "."
    ))
  a = array(0, c(s, M, M))
  for (j in 1:s) {
    for (k in 1:M)
      a[j, , k] = .rsobol.int2bits(col[j, k], M)
  }

  n = 2^m
  ans = matrix(0, n, s)
  bits = rep(0, M)
  for (i in 1:n) {
    bitsi = .rsobol.int2bits(i - 1)
    for (j in 1:s) {
      bitsj = a[j, , ] %*% bitsi
      for (k in 1:M)
        bitsj[k] = bitsj[k] %% 2
      ans[i, j] = .rsobol.bits2unif(bitsj)
    }
  }
  ans
}


# Below are functions to test whether things are going right
# They do not furnish a proof of correctness, but they can catch many errors.

.rsobol.testrate = function(mset = 5:18,
                            R = 50,
                            type = "nestu",
                            seed = 20210120) {
  # For a smooth function in modest dimensions we should see
  # nearly 1/n^3 variance.
  #
  # This function puts an example plot on screen.
  #
  # This function plots some empirical variances along with
  # reference curves as sample size varies.  The selected
  # integrand is smooth but not antithetic or polynomial
  # or spiky or in any way sensitive to base b=2.

  # Test integrand on [0,1]^2; smooth but not specially
  # tuned for Sobol points
  g = function(x) {
    d = length(x)
    exp(sum(x)) - (exp(1) - 1)^d
  }

  cummean = function(v) {
    cumsum(v) / (1:length(v))
  } # Cumulative means

  vals = matrix(0, R, length(mset))
  n = 2^mset
  for (r in 1:R) {
    x = rsobol(
      m = max(mset),
      s = 2,
      type = type,
      seed = seed + r
    )
    y = apply(x, 1, g)
    vals[r, ] = cummean(y)[n]
  }
  varhat = apply(vals^2, 2, mean)
  plot(
    n,
    varhat,
    xlab = "n",
    ylab = "Estimated variance",
    log = "xy",
    main = "Reference curves through first point\nparallel to 1/n^3 and log(n)/n^3"
  )

  lines(n,
        (1 / n^3) * varhat[1] * n[1]^3,
        col = "blue",
        lwd = 2)
  lines(n,
        (log(n) / n^3) * varhat[1] * n[1]^3 / log(n[1]),
        col = "red",
        lwd = 2)
}


.rsobol.canvas = function(...) {
  # Empty uncluttered plot initialization, aspect ratio is 1
  plot(
    c(0, 1),
    c(0, 1),
    type = "n",
    asp = 1,
    xlab = "",
    ylab = ""
  )
}

.rsobol.glines = function(hn, vn, ...) {
  # Apply uniform grid lines in [0,1]^2 to check elementary intervals
  for (i in 0:hn)
    lines(c(0, 1), rep(i / hn, 2), ...)
  for (j in 0:vn)
    lines(rep(j / vn, 2), c(0, 1), ...)
}

.rsobol.testgrids = function(fn = "rsobolgridtest.pdf", m = 8) {
  # This plots some pairs of variables with grid lines
  # Grid cells do not always have one point each because
  # the Sobol' nets may have t>0
  #
  print(paste("Sending test plots to", fn))

  grDevices::pdf(fn, 6, 7)
  par(xaxt = "n", yaxt = "n", bty = "n")
  par(mar = c(0, 0, 0, 0) + .1)
  par(oma = c(0, 0, .5, 0))

  dims = c(1, 2, 5, 10, 20, 30, 40, 50)
  vals = rsobol(m = m, s = 50, seed = 1)[, dims]

  s = ncol(vals)
  for (i in 1:(s - 1))
    for (j in (i + 1):s)
      for (k in 1:(m - 1)) {
        .rsobol.canvas()
        points(vals[, i], vals[, j], pch = 20)
        mtext(
          outer = TRUE,
          side = 3,
          line = -2,
          paste(sep = "", "Input ", dims[j], " vs ", dims[i])
        )
        n1 = 2^k
        n2 = 2^(m - k)
        .rsobol.glines(n1, n2)
      }

  grDevices::dev.off()
}

.rsobol.verifylhs = function(mset = 3:12,
                             s = 1000,
                             type = "nestu",
                             verbose = TRUE) {
  # Verify that each variable is stratified

  lenu = function(v) {
    length(unique(v))
  } # number of unique elements in a vector
  allok = TRUE
  for (m in mset) {
    if (verbose)
      print(paste("Doing m =", m))
    n = 2^m
    vals = rsobol(
      fn = "sobol_Cs.col",
      m = m,
      s = s,
      type = type,
      seed = 1
    )
    vals = floor(n * vals) # should be 0,1, ... , 2^m-1 in every column
    lenus = apply(vals, 2, lenu)
    if (any(lenus != n)) {
      allok = FALSE
      print(paste("LHS problems at m =", m))
      print(which(lenus != n))
    }
    if (verbose && allok)
      print(paste("All ok up to m =", m))
  }
  if (allok)
    print("LHS properties look ok")

}
