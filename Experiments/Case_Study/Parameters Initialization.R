initial_heuristic <- function(data, states, positive_mu) {
  
  # initial clustering
  cluster <- stats::kmeans(
    data[!is.na(data)], centers = states, iter.max = 100, nstart = 100
  )$cluster
  print("cluster")
  print(cluster)
  
  # Gamma <- matrix(0.1 / (states - 1), nrow = states, ncol = states)
  # diag(Gamma) <- 0.9
  
  
  Gamma <- matrix(NA, states, states)
  for (i in 1:states) {
    indices_i <- which(cluster[1:(length(cluster) - 1)] == i)
    if (length(indices_i) == 0) {
      next
    }
    
    N_i <- length(indices_i) 
    
    for (j in 1:states) {
      N_ij <- sum(cluster[indices_i + 1] == j)
      N_ij <- max(N_ij, 1e-3)
      Gamma[i, j] <- N_ij / N_i
    }
  }
  print(Gamma)
  
  kappa_est <- numeric(states)
  theta_est <- numeric(states)
  sigma_est <- numeric(states)
  
  
  dt <- 1/252
  
  Delta_Xt <- diff(data)
  Xt_lag <- data[-length(data)]
  
  cluster_lag <- cluster[-length(cluster)]
  
  for (s in seq_len(states)) {
    
    Xt_s <- data[cluster == s]
    
    # theta (numerical protection)
    theta_s <- mean(Xt_s, na.rm = TRUE)
    theta_est[s] <- max(theta_s, 1e-6)
    
    
    
    indices_s <- which(cluster_lag == s)
    
    if (length(indices_s) < 3) {
      warning(paste("State", s, "has too few observations for OLS. Using default kappa = 1.0"))
      kappa_est[s] <- 2
      next 
    }
    
    Delta_Xt_s <- Delta_Xt[indices_s]
    Xt_lag_s <- Xt_lag[indices_s]

    # sigma (numerical protection)
    sigma_s_base <- sd(Xt_s, na.rm = TRUE)
    sigma_est[s] <- sqrt(var(Delta_Xt_s, na.rm = TRUE) / (dt * mean(Xt_s, na.rm = TRUE)))
    
    
    # Fit the OLS model
    ols_model_s <- lm(Delta_Xt_s ~ Xt_lag_s)
    
    # Scatter plot
    plot(Xt_lag_s, Delta_Xt_s,
         xlab = "Xt_lag_s",
         ylab = "Delta_Xt_s",
         main = "Scatter Plot of ΔXt vs Xt_lag",
         pch = 19,
         col = "blue")
    
    # Add fitted regression line
    abline(ols_model_s, col = "red", lwd = 2)
    
    residuals_s <- residuals(ols_model_s)
    
    # --- Plot 1: Residuals vs Fitted values ---
    plot(fitted(ols_model_s), residuals_s,
         xlab = "Fitted values",
         ylab = "Residuals",
         main = "Residuals vs Fitted",
         pch = 19,
         col = "darkgreen")
    abline(h = 0, col = "red", lwd = 2)
    
    # --- Plot 2: Autocorrelation of residuals ---
    acf(residuals_s,
        main = "Autocorrelation of Residuals",
        col = "blue",
        lwd = 2)
    
    beta_hat <- coef(ols_model_s)["Xt_lag_s"]
    alpha_hat <- coef(ols_model_s)["(Intercept)"]
    
    
    # print(paste("kappa (by alpha):",alpha_hat/(dt*theta_s)))
    kappa_s_raw <- -beta_hat / dt
    print(paste("kappa (by beta):",kappa_s_raw))
    # numerical protection
    kappa_est[s] <- max(kappa_s_raw, 1e-3)
    # kappa_est[s] <- 10
  }
  
  list(
    "cluster" = cluster,
    "pars" = list(
      "kappa" = kappa_est, 
      "theta" = theta_est, 
      "sigma" = sigma_est, 
      "Gamma" = Gamma
    )
  )
}




set.seed(333)
N <- 250
v0 <- 0.03
S0 <- 100

Reg_chain <- simulate_Reg(series_length = N, Reg_tran = matrix(c(1, 0, 0, 1), 2, 2))




param_one_reg <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta, sigma, rho
    0.5,   10,     0.03,  0.05,  0, # calm
    0.5,   10,     0.03,  0.05,  0 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)


# 62-144 regime change
sim_series  <- simulate_heston(S0, v0, Reg_chain, param_one_reg, T, N, M=1, method = "M")
data <- sim_series$V_paths[1,]
plot(data, type = "l")
# plot(sim_series$S_paths[1,])
initial_heuristic(data, 1, TRUE)


# initial_heuristic(V_simulated, 2, TRUE)



