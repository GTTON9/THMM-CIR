get_init<- function(data, states) {
  
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
    
    # sigma (numerical protection)
    sigma_s_base <- sd(Xt_s, na.rm = TRUE)
    sigma_est[s] <- sigma_s_base
    
    
    indices_s <- which(cluster_lag == s)
    print(indices_s)
    if (length(indices_s) < 3) {
      warning(paste("State", s, "has too few observations for OLS. Using default kappa = 1.0"))
      kappa_est[s] <- 2
      next 
    }
    
    Delta_Xt_s <- Delta_Xt[indices_s]
    Xt_lag_s <- Xt_lag[indices_s]
    
    ols_model_s <- lm(Delta_Xt_s ~ Xt_lag_s)
    
    # question is, there are several sequences of OLS, shall we weigh them and get the estimate?
    
    beta_hat <- coef(ols_model_s)["Xt_lag_s"]
    alpha_hat <- coef(ols_model_s)["(Intercept)"]
    
    print(beta_hat)
    print(paste("kappa (by alpha):",alpha_hat/(dt*theta_s)))
    kappa_s_raw <- -beta_hat / dt
    print(paste("kappa (by beta):",kappa_s_raw))
    # numerical protection
    kappa_est[s] <- max(kappa_s_raw, 1e-3)
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



T =1
N = 250
v0 = 0.03
S0 =100


mu <- matrix(c(), nrow = 1)
kappa <- matrix()
theta <- matrix()
sigma <- matrix()
rho <- matrix()
 

Reg_chain <- simulate_Reg(series_length = N)
plot(Reg_chain)
Reg_param <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta, sigma, rho
    0.5,   10,    0.03,   0.05,  -0.1, # calm
    0.5,   5,     0.6,   0.05,  -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)


# 62-144 regime change
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T = 0.5, N, M=1, method = "M")
S_simulated <- sim_series$S_paths[1,]
plot(S_simulated,  type = "l")

V_simulated <- sim_series$V_paths[1,]
plot(V_simulated, type = "l")


plot(S_simulated, type = "l")
plot(diff(V_simulated), type = "l")
plot(diff(log(S_simulated)), type = "l")








log_returns <- diff(log(S_simulated)) 
volatility_proxy <- log_returns^2 * 252 
plot(V_simulated, type = 'l', ylim = c(0,1))
lines(volatility_proxy, col = "blue")

# volatility_proxy <- log(volatility_proxy)
states <- 2 

cluster_results <- stats::kmeans(
  volatility_proxy[!is.na(10 * volatility_proxy)], 
  centers = states, 
  iter.max = 100, 
  nstart = 100 
)
cluster_assignment <- cluster_results$cluster 
plot(Reg_chain+1, col = "blue")
lines(cluster_assignment)










# --- 1. Define Input Data and Time Step ---

# V_simulated: Your input volatility series (standard deviation)
# Assuming V_simulated is already loaded, e.g.,
# V_simulated <- sqrt(V_simulated_variance)

# Determine the time step (assuming daily data, t is in years)
Delta_t <- 1 / 252 

# --- 2. Calculate the Variance Proxy (v_t) ---

# The Heston model uses variance (v_t), which is (volatility)^2
v_t <- V_simulated

# --- 3. Prepare Lagged Data for Regression ---

# Dependent variable: v_{i+1}
v_i_plus_1 <- v_t[2:length(v_t)]

# Independent variable: v_i
v_i <- v_t[1:(length(v_t) - 1)]

# --- 4. Run the Linear Regression ---

# Regression model: v_{i+1} = alpha + beta * v_i
model <- lm(v_i_plus_1 ~ v_i)

# Extract coefficients
alpha_hat <- coef(model)[1] # Intercept (alpha)
beta_hat <- coef(model)[2]  # Slope (beta)

# --- 5. Calculate Heston Parameter Estimators ---

# Heston Kappa (\kappa) Estimate (Mean-reversion rate)
# Formula: \hat{\kappa} = (1 - \hat{\beta}) / \Delta t
kappa_hat <- (1 - beta_hat) / Delta_t

# Heston Theta (\theta) Estimate (Long-term mean variance)
# Formula: \hat{\theta} = \hat{\alpha} / (\hat{\kappa} * \Delta t)
# Check for stability before dividing by kappa_hat
if (kappa_hat > 1e-6) {
  theta_hat <- alpha_hat / (kappa_hat * Delta_t)
} else {
  # If kappa is near zero or negative, the process is unstable,
  # and theta estimate is unreliable.
  theta_hat <- alpha_hat / (1 - beta_hat) # Alternative division
}


# --- 6. Output Results and Sanity Check ---

cat("--- Heston Initial Estimator Results from Volatility Series ---\n")
cat(sprintf("Time Step (Delta_t): %.6f\n", Delta_t))
cat(sprintf("Regression Beta (\u03B2\u0302): %.6f\n", beta_hat))
cat(sprintf("Regression Alpha (\u03B1\u0302): %.10f\n", alpha_hat))
cat("\n")

if (kappa_hat > 0) {
  cat(sprintf("Initial Kappa Estimate (\u03BA\u0302): %.4f (per year)\n", kappa_hat))
  cat(sprintf("Initial Theta Estimate (\u03B8\u0302): %.8f (Variance, annualised)\n", theta_hat))
  cat(sprintf("Implied Long-Term Volatility (\u221A\u03B8\u0302): %.4f\n", sqrt(theta_hat)))
} else {
  cat(sprintf("Warning: Estimated \u03BA\u0302 (%.4f) is non-positive or near zero. \n", kappa_hat))
  cat("This suggests the mean-reversion assumption may not hold for this data.\n")
  # Provide a stable, albeit arbitrary, fallback for optimization initialization
  kappa_hat_fallback <- 0.5
  theta_hat_fallback <- mean(v_t)
  cat(sprintf("Using Fallback Initial Estimates: \u03BA = %.2f, \u03B8 = %.8f\n", 
              kappa_hat_fallback, theta_hat_fallback))
}












# checking the parameter sensitivity
# plot the transition distribution for different set of parameters


library(ggplot2)


source("Model_Simulation.R")

set.seed(333)
N <- 250
v0 <- 0.04
S0 <- 100

Reg_chain <- simulate_Reg(series_length = N)




Reg_param_theta <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta, sigma, rho
    0.5,   10,     0.03,  0.05,  -0.1, # calm
    0.5,   5,     0.15,  0.05,  -0.1 # turbulent
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



plot(S_simulated, type = "l")
plot(V_simulated, type = "l")
log_returns <- diff(log(S_simulated)) 
volatility_proxy <- log_returns^2 
plot(volatility_proxy)
lines(V_simulated^2)
# volatility_proxy <- log(volatility_proxy)
states <- 2 

cluster_results <- stats::kmeans(
  volatility_proxy[!is.na(10 * volatility_proxy)], 
  centers = states, 
  iter.max = 100, 
  nstart = 100 
)
cluster_assignment <- cluster_results$cluster 
plot(Reg_chain+1, col = "blue")
lines(cluster_assignment)

  
  
  
  
  