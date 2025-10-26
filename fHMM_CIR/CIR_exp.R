# SDE discretization to simulate Heston model
simulate_heston <- function(S0, v0, Reg_series, Reg_param, T, N, M, method = "E") {
  set.seed(999)
  if (2 * Reg_param[1,2] * Reg_param[1,3] < Reg_param[1,4]^2 || 
      2 * Reg_param[2,2] * Reg_param[2,3] < Reg_param[2,4]^2) {
    message("Warning: Feller condition (2*kappa*theta > sigma^2) NOT satisfied in at least one regime.")
  }
  
  dt <- T / N
  sqrt_dt <- sqrt(dt)
  
  method <- toupper(method) 
  
  S_paths <- matrix(0, nrow = M, ncol = N + 1)
  V_paths <- matrix(0, nrow = M, ncol = N + 1)
  
  S_paths[, 1] <- S0
  V_paths[, 1] <- v0
  
  for (i in 1:N) {
    reg_index <- Reg_series[i] + 1 # 0->1, 1->2 (different parameter by regime)
    mu_curr <- Reg_param[reg_index, 1]
    kappa_curr <- Reg_param[reg_index, 2]
    theta_curr <- Reg_param[reg_index, 3]
    sigma_curr <- Reg_param[reg_index, 4] 
    rho_curr <- Reg_param[reg_index, 5]
    
    rho_prime <- sqrt(1 - rho_curr^2)
    # stochastic term
    Z1 <- rnorm(M) 
    Z2 <- rnorm(M) 
    
    V_prev <- V_paths[, i]
    S_prev <- S_paths[, i]
    
    dW1 <- Z1 * sqrt_dt
    dW2 <- (rho_curr * Z1 + rho_prime * Z2) * sqrt_dt
    
    sqrt_V <- sqrt(V_prev)
    
    # check negative/zero variance
    
    if (method == "E") { # Euler-Maruyama Scheme
      
      V_next <- (
        V_prev 
        + kappa_curr * (theta_curr - V_prev) * dt     
        + sigma_curr * sqrt_V * dW2                  
      )
      
    } else if (method == "M") {
      V_euler_term <- (
        V_prev 
        + kappa_curr * (theta_curr - V_prev) * dt
        + sigma_curr * sqrt_V * dW2
      )
      
      V_milstein_term <- 0.25 * sigma_curr^2 * (dW2^2 - dt)
      
      V_next <- V_euler_term + V_milstein_term # 修正: 使用 sigma_curr 代替 xi
    }
    
    V_next <- pmax(V_next, 0)
    
    V_paths[, i + 1] <- V_next
    S_paths[, i + 1] <- (
      S_prev 
      + mu_curr * S_prev * dt                       
      + S_prev * sqrt_V * dW1                  
    )
  }
  
  return(list(S_paths = S_paths, V_paths = V_paths, method_used = method))
}






# SDE discretization to simulate Heston model with different regime for i and j
simulate_heston_diff <- function(S0, v0, Reg_series, A_series, Reg_param, T, N, M, method = "E") {
  set.seed(999)
  if (2 * Reg_param[1,2] * Reg_param[1,3] < Reg_param[1,4]^2 || 
      2 * Reg_param[2,2] * Reg_param[2,3] < Reg_param[2,4]^2) {
    message("Warning: Feller condition (2*kappa*theta > sigma^2) NOT satisfied in at least one regime.")
  }
  
  dt <- T / N
  sqrt_dt <- sqrt(dt)
  
  method <- toupper(method) 
  
  S_paths <- matrix(0, nrow = M, ncol = N + 1)
  V_paths <- matrix(0, nrow = M, ncol = N + 1)
  
  S_paths[, 1] <- S0
  V_paths[, 1] <- v0
  
  for (i in 1:N) {
    reg_index <- Reg_series[i] + 1 # 0->1, 1->2 (different parameter by regime)
    A_index <- A_series[i] + 1
    mu_curr <- Reg_param[A_index, 1]
    kappa_curr <- Reg_param[reg_index, 2]
    theta_curr <- Reg_param[reg_index, 3]
    sigma_curr <- Reg_param[reg_index, 4] 
    rho_curr <- Reg_param[reg_index, 5]
    
    rho_prime <- sqrt(1 - rho_curr^2)
    # stochastic term
    Z1 <- rnorm(M) 
    Z2 <- rnorm(M) 
    
    V_prev <- V_paths[, i]
    S_prev <- S_paths[, i]
    
    dW1 <- Z1 * sqrt_dt
    dW2 <- (rho_curr * Z1 + rho_prime * Z2) * sqrt_dt
    
    sqrt_V <- sqrt(V_prev)
    
    # check negative/zero variance
    
    if (method == "E") { # Euler-Maruyama Scheme
      
      V_next <- (
        V_prev 
        + kappa_curr * (theta_curr - V_prev) * dt     
        + sigma_curr * sqrt_V * dW2                  
      )
      
    } else if (method == "M") {
      V_euler_term <- (
        V_prev 
        + kappa_curr * (theta_curr - V_prev) * dt
        + sigma_curr * sqrt_V * dW2
      )
      
      V_milstein_term <- 0.25 * sigma_curr^2 * (dW2^2 - dt)
      
      V_next <- V_euler_term + V_milstein_term 
    }
    
    V_next <- pmax(V_next, 0)
    
    V_paths[, i + 1] <- V_next
    S_paths[, i + 1] <- (
      S_prev 
      + mu_curr * S_prev * dt                       
      + S_prev * sqrt_V * dW1                  
    )
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

  








set.seed(999)
N <- 250
v0 <- 0.04
S0 <- 100

Reg_chain <- simulate_Reg(series_length = N)
plot(Reg_chain)



Reg_param_theta <- matrix(
  c( # mu, kappa, theta, sigma, rho
    0.05, 3.0, 0.03, 0.15, -0.1, 0.05, 0.5, 0.2, 0.15, -0.1
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)



sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param_theta, T, N, M=2, method = "M")
S_simulated <- sim_series$S_paths[1,]
# S_mu_model <- fit_HMM(S_simulated, Reg_chain)
# summary(S_mu_model)
# S_mu_model


V_simulated <- sim_series$V_paths[1,]
# V_mu_model <- fit_HMM(V_simulated, Reg_chain)





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
source("LL_HMM_R.R")

# followed the example
series_control <- Heston_set_controls( 
  states      = 2,     # 2 state
  sdds        = "Heston",         
  date_column = "Date",
  file        = my_data_df, 
  data_column = "price",      
  logreturns  = FALSE,         
  from        = date_sequence[1],             
  to          = date_sequence[length(date_sequence)],
  runs = 100
)


data_hmm <- prepare_data(series_control)


model_hmm <- Heston_fit_model(data_hmm) # fit HMM


model_hmm <- decode_states(model_hmm)



Heston_fit_model(S_simulated, Reg_chain, "Heston")