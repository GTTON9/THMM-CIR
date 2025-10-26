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



