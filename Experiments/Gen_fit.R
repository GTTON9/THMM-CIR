source("Model_Simulation.R")
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
source("Heston_likelihood.R")

Gen_fit <- function(mu, kappa, theta, sigma, rho, gen_length = 250, V_0 = 0.03, S0 = 100,  V_ = TRUE , plot_path = TRUE){
  # V_: whether we know the volatility process
  
  set.seed(999)
  Reg_chain <- simulate_Reg(series_length = N)
  Reg_param <- cbind(mu, kappa, theta, sigma, rho)
  sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E_C")
  S_simulated <- sim_series$S_paths[1,]
  
  V_simulated <- sim_series$V_paths[1,]
  
  
  par(mfrow = c(1, 2))
  if(plot_path){
    plot(
      S_simulated,
      type = "l",
      main = "S path",
      col = "blue")
    
    plot(
      V_simulated,
      type = "l",
      main = "V path",
      col = "red"
    )
  }
  par(mfrow = c(1, 1))
  
  if(V_){
    
  }
  
  
  series_length <- length(V_simulated)
  start_date <- as.Date("2024-01-01")
  date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
  
  
  
  my_data_df <- data.frame(
    Date = date_sequence,
    price = V_simulated
  )
  
  
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
  
  plot(states_estimate, col = "blue")
  lines(Reg_chain+1)
  param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
  
  
  return(list(states_estimate = states_estimate, 
              param = param))
}



mu <- c(0.5, 0.5)
kappa <- c(10, 5)
theta <- c(0.03, 0.6)
sigma <- c(0.05, 0.05)
rho <- c(-0.1, -0.1)


Gen_fit(mu, kappa, theta, sigma, rho, V_ = TRUE , plot_path = TRUE)




