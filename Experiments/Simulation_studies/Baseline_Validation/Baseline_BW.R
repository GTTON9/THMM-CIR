source("Gen_fit.R")

N <- 250
v0 <- 0.03
S0 <- 100
nstates = 2
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(10, 5)
theta <- c(0.03, 0.6)
sigma <- c(0.05, 0.05)
rho <- c(-0.1, -0.1)

BW_result <- Gen_fit(mu, kappa, theta, sigma, rho, V_ = TRUE , plot_path = TRUE)
BW_result$states_estimate








