source("Fit_HMM.R")
source("simulate_heston.R")
source("simulate_heston.R")
set.seed(999)

N <- 250
Reg_chain <- simulate_Reg(series_length = N)
plot(Reg_chain)



Reg_param_theta_sig <- matrix(
  c(# mu, kappa, theta, sigma, rho
    0.05, 2.0, 0.03, 0.15, -0.1,  0.05, 2.0, 0.2, 0.6, -0.1
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)



sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param_theta_sig, T, N, M=2, method = "M")
S_simulated <- sim_series$S_paths[1,]
plot(S_simulated, type = "l")

S_model <- fit_HMM(S_simulated, Reg_chain)
S_model$estimate
get_HMM_param(S_model)

V_simulated <- sim_series$V_paths[1,]
V_model <- fit_HMM(V_simulated, Reg_chain, sdds = "gamma")
get_HMM_param(V_model)
summary(V_model)

plot(V_simulated, type = "l")








