
source("Fit_HMM.R")
source("simulate_heston.R")

N <- 250
set.seed(999)
Reg_chain <- simulate_Reg(series_length = N)
plot(Reg_chain)



Reg_param_mu <- matrix( # only mu is regime dependent
  c( # mu, kappa, theta, sigma, rho
    0.05, 2.0, 0.03, 0.15, -0.1,  5, 2.0, 0.03, 0.15, -0.1
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)

sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param_mu, T, N, M=2, method = "M")
S_simulated <- sim_series$S_paths[1,]
S_model <- fit_HMM(S_simulated, Reg_chain)

get_HMM_param(S_model)


V_simulated <- sim_series$V_paths[1,]
V_model <- fit_HMM(V_simulated, Reg_chain)

get_HMM_param(V_model)









