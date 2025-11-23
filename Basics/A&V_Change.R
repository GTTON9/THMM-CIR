source("Fit_HMM.R")
source("simulate_heston.R")
library(fHMM)
set.seed(999)
Reg_chain <- simulate_Reg(series_length = 250)
plot(Reg_chain)



Reg_param_AnV <- matrix(
  c(
    0.05, 2.0, 0.03, 0.15, -0.1, 3, 2.0, 0.2, 0.6, -0.1
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)



sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param_AnV, T, N, M=2, method = "M")
S_simulated <- sim_series$S_paths[1,]
S_mu_model <- fit_HMM(S_simulated, Reg_chain)



V_simulated <- sim_series$V_paths[1,]
V_mu_model <- fit_HMM(V_simulated, Reg_chain)

  