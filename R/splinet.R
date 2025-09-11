data_prepare <- function(f_data, t_data) {
  colnames(f_data) <- NULL
  colnames(t_data) <- NULL
  f_data <- t(f_data)
  ready_data <- cbind(t_data, f_data)
  ready_data <- as.matrix(ready_data)
  return(ready_data)
}

# Prepare knots
Knots_prepare <- function(selected_knots, Time) {
  knots_normalized <- selected_knots / max(selected_knots)
  knots_normalized = knots_normalized * (max(Time) - min(Time)) + min(Time)
  # taking the first three decimal number
  knots_normalized = as.numeric(format(round(knots_normalized, 4), nsmall = 1))
  return(knots_normalized)
}


basisSplinet <- function(Ystandard, nBasis) {
  nMV = ncol(Ystandard)
  is_Splinets_available <- requireNamespace("Splinets", quietly = TRUE)
  if (!is_Splinets_available) {
    stop("For basisType=='splinets', install Splinets")
  }
  is_DDK_available <- requireNamespace("DDK", quietly = TRUE)
  if (!is_DDK_available) {
    stop(
      "For basisType=='splinets', install DDK (devtools::install_github('ranibasna/ddk'))"
    )
  }
  # knot selection
  initial_knots <- c(0, nMV)
  # selecting the knots with add_knots function
  KS <- DDK::add_knots(f = Ystandard,
                       knots = initial_knots,
                       L = 20,
                       M = 1)
  DDKnots <- Knots_prepare(selected_knots = KS[[1]],
                           Time = seq(0, 1, length.out = nMV))
  f_ready_data <- data_prepare(f_data = Ystandard, t_data = seq(0, 1, length.out = nMV))
  ProjObj <- Splinets::project(f_ready_data, DDKnots)
  basis = t(Splinets::evspline(ProjObj$basis, x = seq(0, 1, length.out = nMV))[, -1])

  return(basis)
}
