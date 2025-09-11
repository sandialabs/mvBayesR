library(fdasrvf)
f <- function(x) {
  dnorm(seq(0, 1, length.out = 99),
        sin(2 * pi * x[1] ^ 2) / 4 - x[1] / 10 + 0.5,
        0.05) * x[2]
}

n = 100
M = 99
p = 3
sim_variables = matrix(runif(n * p), n)
x_test = matrix(runif(1000 * p), 1000)
e = rnorm(n * 99)
f_sim = matrix(0, n, M)
for (i in 1:n) {
  f_sim[i, ] = f(sim_variables[i, ])
}
y_test = matrix(0, 1000, M)
for (i in 1:1000) {
  y_test[i, ] = f(x_test[i, ])
}

x_true = c(0.1028, 0.5930)
tt = seq(0, 1, length.out = M)
f_exp = f(x_true)

tt = seq(0, 1, length.out = M)
warp_list = multiple_align_functions(
  t(f_sim),
  tt,
  f_exp,
  lambda = .01,
  parallel = TRUE,
  showplot = FALSE,
  verbose = FALSE
)
gam_sim = warp_list$warping_functions
psi_sim = gam_to_psi(warp_list$warping_functions)
ftilde_sim = warp_list$fn

psi_obs = gam_to_psi(seq(0, 1, length.out=M))
gam_obs = seq(0, 1, length.out=M)
ftilde_obs = f_exp

emu_vv = mvBayesElastic(
  BASS::bass,
  sim_variables,
  t(ftilde_sim),
  warpData = warp_list
)

tmp = predict(emu_vv, sim_variables)

plot(emu_vv)
