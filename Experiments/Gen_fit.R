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
source("Viterbi_Visual.R")

Gen_fit <- function(Gamma, mu, kappa, theta, sigma, rho, gen_length = 252, V_0 = 0.03, S0 = 100,  V_ = TRUE , plot_path = TRUE, seed = 999){
  # V_: whether we know the volatility process
  
  set.seed(999)
  N <- gen_length
  Reg_chain <- simulate_Reg(series_length = N, Reg_tran = Gamma)
  Reg_param <- cbind(mu, kappa, theta, sigma, rho)
  sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E", seed = seed)
  S_simulated <- sim_series$S_paths
  print(S_simulated)
  V_simulated <- sim_series$V_paths
  
  
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
    Var = V_simulated
  )
  
  
  series_control <- Heston_set_controls( 
    states      = 2,     # 2 state
    sdds        = "Heston",         
    date_column = "Date",
    file        = my_data_df, 
    data_column = "Var",      
    logreturns  = FALSE,         
    from        = date_sequence[1],             
    to          = date_sequence[length(date_sequence)],
    runs = 10
  )
  
  
  data_hmm <- prepare_data(series_control)
  model_hmm <- Heston_fit_model(data_hmm) 
  
  final_model <- decode_states_heston(model_hmm) 
  states_estimate <- final_model$decoding
  states_estimate

  plot(states_estimate, col = "blue")
  lines(Reg_chain+1)
  param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
  
  
  p <- plot_viterbi(V_simulated, nstates, param$Gamma, 
               param$kappa, param$theta, 
               param$sigma, Reg_chain)
  plot(p)
  
  
  return(list(states_estimate = states_estimate, 
              param = param,
              fisher_inverse = final_model$inverse_fisher))
}
















Gen_fit <- function(Gamma, mu, kappa, theta, sigma, rho, input_ = "V_", gen_length = 252, V_0 = 0.03, S0 = 100,  V_ = TRUE , plot_path = TRUE, seed = 999){
  # V_: whether we know the volatility process
  
  set.seed(999)
  
  
  if(input_ == "V_"){
    N <- gen_length
    Reg_chain <- simulate_Reg(series_length = N, Reg_tran = Gamma)
    Reg_param <- cbind(mu, kappa, theta, sigma, rho)
    sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E", seed = seed)
    S_simulated <- sim_series$S_paths
    print(S_simulated)
    V_simulated <- sim_series$V_paths
    
    
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
    
    
    series_length <- length(S_simulated * 100)
    start_date <- as.Date("2024-01-01")
    date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
    
    
    
    my_data_df <- data.frame(
      Date = date_sequence,
      Var = V_simulated
    )
    
    
    series_control <- Heston_set_controls( 
      states      = 2,     # 2 state
      sdds        = "Heston",         
      date_column = "Date",
      file        = my_data_df, 
      data_column = "Var",      
      logreturns  = FALSE,         
      from        = date_sequence[1],             
      to          = date_sequence[length(date_sequence)],
      runs = 10
    )
    
    
    data_hmm <- prepare_data(series_control)
    model_hmm <- Heston_fit_model(data_hmm) 
    
    final_model <- decode_states_heston(model_hmm) 
    states_estimate <- final_model$decoding
    states_estimate
    
    plot(states_estimate, col = "blue")
    lines(Reg_chain+1)
    param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
    
    
    p <- plot_viterbi(V_simulated, nstates, param$Gamma, 
                      param$kappa, param$theta, 
                      param$sigma, Reg_chain)
    plot(p)
    
  }else if(input_ == "S_"){
    
    N <- gen_length * 100
    Reg_chain_year <- simulate_Reg(series_length = N/100, Reg_tran = Gamma)
    Reg_chain <- rep(Reg_chain_year, each = 100)
    
    Reg_param <- cbind(mu, kappa, theta, sigma, rho)
    sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E", seed = seed)
    S_simulated <- sim_series$S_paths
    print(S_simulated)
    V_simulated <- sim_series$V_paths
    
    
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
    
    
    
    
    S <- S_simulated
    n_days <- 250
    n_intraday <- 100
    
    RV_V <- numeric(n_days)
    
    for (t in 1:n_days) {
      idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
      S_day <- S[idx]
      r_day <- diff(log(S_day))      # length 3
      RV_V[t] <- sum(r_day^2) * 250 
    }
    
    
    
    plot(V_simulated[seq(1, length(V_simulated), by = 100)], type = 'l', ylim = c(0,1), col = "black")
    lines(RV_V, col = "blue")
    lines(lowess(RV_V, f = 0.1), col = "red")
    legend("topright", 
           legend = c("True Variance (V_t)", 
                      "Realized Variance Estimate", 
                      "Smoothed RV (lowess, f=0.1)"),
           col = c("black", "blue", "red"),
           lwd = c(2, 1.5, 2),
           lty = c(1, 1, 2),
           bty = "n",
           cex = 0.5)
    
    
    
    RV_V <- lowess(RV_V, f = 0.1)$y
    
    
    
    
    series_length <- length(RV_V)
    start_date <- as.Date("2024-01-01")
    date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
    
    
    my_data_df <- data.frame(
      Date = date_sequence,
      Var = RV_V
    )
    colnames(my_data_df) <- c("Date", "Var")
    
    
    source("Gen_fit.R")
    
    # followed the example
    series_control <- Heston_set_controls( 
      states      = 2,     # 2 state
      sdds        = "Heston",         
      date_column = "Date",
      file        = my_data_df, 
      data_column = "Var",      
      logreturns  = FALSE,         
      from        = date_sequence[1],             
      to          = date_sequence[length(date_sequence)],
      runs = 10
    )
    
    
    data_hmm <- prepare_data(series_control)
    model_hmm <- Heston_fit_model(data_hmm) 
    
    final_model <- decode_states_heston(model_hmm) 
    states_estimate <- final_model$decoding
    states_estimate
    
    plot(states_estimate, col = "blue")
    lines(Reg_chain_year+1)
    param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
    
    
    p <- plot_viterbi(V_simulated, nstates, param$Gamma, 
                      param$kappa, param$theta, 
                      param$sigma, Reg_chain_year)
    plot(p)
  }
  
  
  
  
  
  return(list(states_estimate = states_estimate, 
              param = param,
              fisher_inverse = final_model$inverse_fisher))
}




