
simulate_heston <- function(S0, v0, Reg_series, Reg_param, T, N, M, method = "E", interp = TRUE, seed = 999, min_var = 1e-6) {
  if (!is.null(seed)) set.seed(seed)
  
  # Check Feller (just warn)
  if (2 * Reg_param[1,2] * Reg_param[1,3] < Reg_param[1,4]^2 || 
      2 * Reg_param[2,2] * Reg_param[2,3] < Reg_param[2,4]^2) {
    message("Warning: Feller condition (2*kappa*theta > sigma^2) NOT satisfied in at least one regime.")
  }
  if(interp){
    N <-  N * 100
    Reg_series <- rep(Reg_series, each = 100)
  }
  
  dt <- (T / N) 

  sqrt_dt <- sqrt(dt)
  method <- toupper(method)
  
  S_paths <- matrix(0, nrow = M, ncol = N + 1)
  V_paths <- matrix(0, nrow = M, ncol = N + 1)
  S_paths[, 1] <- S0
  V_paths[, 1] <- v0
  for (i in 1:N) {
    reg_index <- Reg_series[i] + 1
    mu_curr <- Reg_param[reg_index, 1]
    kappa_curr <- Reg_param[reg_index, 2]
    theta_curr <- Reg_param[reg_index, 3]
    sigma_curr <- Reg_param[reg_index, 4] 
    rho_curr <- Reg_param[reg_index, 5]
    
    rho_prime <- sqrt(max(0, 1 - rho_curr^2))
    
    # draw normals
    Z1 <- rnorm(M)
    Z2 <- rnorm(M)
    dW1 <- Z1 * sqrt_dt
    dW2 <- (rho_curr * Z1 + rho_prime * Z2) * sqrt_dt
    
    V_prev_raw <- V_paths[, i]
    # avoid negative arguments to sqrt
    V_prev_pos <- pmax(V_prev_raw, 0)
    sqrt_V_prev <- sqrt(V_prev_pos)
    if (method == "E_C") {
      
      dt_curr <- dt
      kappa <- kappa_curr
      theta <- theta_curr
      sigma <- sigma_curr
      
      # df = 4*kappa*theta / sigma^2
      df_val <- 4 * kappa * theta / sigma^2
      
      # C = 2*kappa / ((1-exp(-kappa*dt))*sigma^2)
      exp_k_dt <- exp(-kappa * dt_curr)
      C_val <- (2 * kappa) / ((1 - exp_k_dt) * (sigma^2))
      
      # lambda = 2*u = 2 * C*V_t*exp(-kappa*dt)
      lambda_val <- 2 * C_val * V_prev_pos * exp_k_dt
      
      Y_sample <- rchisq(M, df = df_val, ncp = lambda_val)
      
      # 3. V_next
      V_next <- Y_sample / (2 * C_val)
      
    } else if (method == "E") {
      V_next <- V_prev_raw + kappa_curr * (theta_curr - V_prev_raw) * dt + sigma_curr * sqrt_V_prev * dW2
    } else if (method == "M") {
      V_euler_term <- V_prev_raw + kappa_curr * (theta_curr - V_prev_raw) * dt + sigma_curr * sqrt_V_prev * dW2
      V_milstein_term <- 0.25 * sigma_curr^2 * (dW2^2 - dt)
      V_next <- V_euler_term + V_milstein_term
    } else {
      stop("Unknown method. Use 'E', 'M', or 'E_C'.")
    }
    
    
    # enforce a nonnegative floor to avoid NaNs later
    V_next <- pmax(V_next, min_var)
    V_paths[, i + 1] <- V_next
    
    # asset dynamics - use variance at previous step (explicit)
    S_prev <- S_paths[, i]
    # S_paths[, i + 1] <- S_prev + mu_curr * S_prev * dt + S_prev * sqrt_V_prev * dW1
    
    # S_paths[, i + 1] <- S_prev * (1 + mu_curr * dt + sqrt_V_prev * dW1)
    S_paths[, i + 1] <- S_prev * exp((mu_curr - 0.5 * V_prev_pos) * dt + sqrt_V_prev * dW1)
    
    }
  if(interp){
    sample_indices <- seq(1, N + 1, by = 100)
    S_paths<- S_paths[, sample_indices]
    V_paths <- V_paths[, sample_indices]
  }else{
    S_paths<- as.vector(S_paths)
    V_paths <- as.vector(V_paths)
  }

  return(list(S_paths = S_paths, V_paths = V_paths, method_used = method))
}






simulate_Reg <- function(
    series_length = 250,
    initial_state = 0,
    Reg_tran = matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2),
    plot = T,
    seed = 999
) {
  set.seed(seed)
  Reg_chain <- numeric(series_length)
  Reg_chain[1] <- initial_state
  
  for (t in 2:series_length) {
    current_state <- Reg_chain[t-1]
    matrix_row_index <- current_state + 1 
    transition_probabilities <- Reg_tran[matrix_row_index, ]
    
    next_state <- sample(
      x = c(0, 1), 
      size = 1, 
      prob = transition_probabilities
    )
    Reg_chain[t] <- next_state
  }
  
  if (plot) {
    plot(Reg_chain, type = 's', 
         ylim = c(-0.1, 1.1),
         main = "Simulated Regime Chain (0=Calm, 1=Turbulent)",
         xlab = "Time Step", 
         ylab = "Regime",
         yaxt = 'n' 
    )
    axis(2, at = c(0, 1), labels = c("Calm", "Turbulent"))
  }
  
  invisible(Reg_chain) 
}



fit_HMM <- function(series, real_reg = NULL, sdds = "t"){
  set.seed(999)
  series_length <- length(series)
  start_date <- as.Date("2024-01-01")
  date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
  
  
  
  my_data_df <- data.frame(
    Date = date_sequence,
    price = series
  )

  # followed the example
  series_control <- set_controls( 
    states      = 2,     # 2 state
    sdds        = sdds,         
    date_column = "Date",
    file        = my_data_df, 
    data_column = "price",      
    logreturns  = TRUE,         
    from        = date_sequence[1],             
    to          = date_sequence[length(date_sequence)]
  )
  
  
  data_hmm <- prepare_data(series_control)

  model_hmm <- fit_model(data_hmm) # fit HMM
  model_hmm <- decode_states(model_hmm)
  
  
  plot(model_hmm$decoding - 1,          
            col = "blue",     
            xlab = "Time Step",  
            ylab = "Regime",           
            main = "Estimated Regime Chain",
            ylim = range(c(model_hmm$decoding - 1, real_reg))
  )
  
  if(! is.null(real_reg)){
    
    lines(real_reg,
          col = "red",     
          type = "l",
          lwd = 3
          
    )
    legend("topleft",              
           legend = c("Estimated Regime", "True Regime"),
           col = c("blue", "red"),
           lty = c(1, 2),
           cex = 0.5
    )
  }
  return(model_hmm)
}








  
