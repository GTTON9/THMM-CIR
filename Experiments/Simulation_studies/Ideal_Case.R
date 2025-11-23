source("Model_Simulation.R")

# perfectly detectable case
# Large parameter differences, low driving noise variance, and a good initialization.
set.seed(999)
N <- 250
v0 <- 1
S0 <- 100

Reg_chain <- simulate_Reg(series_length = N)




Reg_param_theta <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta,  sigma,  rho
    0.5,   10,    0.3,   0.2,   -0.1, # calm
    0.5,   5,     0.6,    0.2,   -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)


# 62-144 regime change
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param_theta, T, N, M=1, method = "M")
S_simulated <- sim_series$S_paths[1,]
plot(S_simulated,  type = "l")
# S_mu_model <- fit_HMM(S_simulated, Reg_chain)
# summary(S_mu_model)
# S_mu_model


V_simulated <- sim_series$V_paths[1,]
# V_mu_model <- fit_HMM(V_simulated, Reg_chain)
plot(V_simulated, type = "l")
# plot(S_simulated, type = "l")
# plot(diff(V_simulated), type = "l")
# plot(diff(log(S_simulated)), type = "l")



set.seed(999)
series_length <- length(V_simulated)
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)



my_data_df <- data.frame(
  Date = date_sequence,
  price = V_simulated
)


source("Heston_controls.R")
source("Heston_data.R")
source("Heston_fit_model.R")
source("Heston_parameters.R")
source("Heston_get_init.R")
source("Heston_Conversion.R")
source("LL_HMM_R.R")
source("Heston_reorder_states.R")
source("fHMM_model.R")
source("Heston_decode.R")

# followed the example
series_control <- Heston_set_controls( 
  states      = 2,     # 2 state
  sdds        = "Heston",         
  date_column = "Date",
  file        = my_data_df, 
  data_column = "price",      
  logreturns  = FALSE,         
  from        = date_sequence[1],             
  to          = date_sequence[length(date_sequence)]
)


data_hmm <- prepare_data(series_control)
model_hmm <- Heston_fit_model(data_hmm) 

final_model <- decode_states_heston(model_hmm) 
states_estimate <- final_model$decoding
states_estimate

plot(states_estimate, type = "s", col = "blue")
lines(Reg_chain+1)
# final estimation
parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)




