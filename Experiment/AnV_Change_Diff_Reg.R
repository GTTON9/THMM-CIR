source("Fit_HMM.R")
source("simulate_heston.R")
source("simulate_Reg.R")

N <- 250
set.seed(999)
V_chain <- simulate_Reg(series_length = N)

A_chain <- simulate_Reg(series_length = N, seed = 299)


Reg_param_AnV <- matrix( # all parameters are different 
  c(
    0.05, 2.0, 0.03, 0.15, -0.1, 5, 2.0, 0.2, 0.6, -0.1
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
) 


set.seed(999)
sim_series  <- simulate_heston_diff(S0, v0, V_chain, A_chain, Reg_param_AnV, T, N, M=2, method = "M")
S_simulated <- sim_series$S_paths[1,]
S_model <- fit_HMM(S_simulated, A_chain)
S_model$estimate

V_simulated <- sim_series$V_paths[1,]
V_model <- fit_HMM(V_simulated, V_chain)
V_model$estimate




