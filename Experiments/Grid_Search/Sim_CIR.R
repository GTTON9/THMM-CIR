install.packages("yuima")
library(yuima)

simCIR(time.points, n, h, alpha, beta, gamma, equi.dist = FALSE)




ln_d_Heston <- function(V_t, V_t_plus_k, k = 1, kappa, theta, sigma) {
  
  # time step
  dt <- k
  
  # CIR parameters for non-central chi-square
  c  <- (sigma^2 * (1 - exp(-kappa * dt))) / (4 * kappa)
  df <- 4 * kappa * theta / sigma^2
  nc <- (4 * kappa * exp(-kappa * dt) / (sigma^2 * (1 - exp(-kappa * dt)))) * V_t
  
  # Safety: force non-negative
  V_t_plus_k[V_t_plus_k < 0] <- 0
  
  # Equivalent noncentral chi-square variable
  x <- V_t_plus_k / c
  
  # log density using dchisq()
  ld <- dchisq(x, df = df, ncp = nc, log = TRUE) - log(c)
  
  return(ld)
}













