source("Heston_likelihood.R")
library(ggplot2)
library(tidyr)
library(plotly)



ln_d_Heston <- function(V_t, V_t_plus_k, k = 1, kappa, theta, sigma) {
  
  C_numerator <- 2 * kappa
  C_denominator <- (1 - exp(-kappa * k)) * (sigma^2)
  
  
  C <- C_numerator / C_denominator
  
  q_numerator <- 2 * kappa * theta
  q_denominator <- sigma^2
  
  q <- (q_numerator / q_denominator) - 1
  
  u <- C * V_t * exp(-kappa * k)
  v <- C * V_t_plus_k
  
  x <- 2 * sqrt(u * v)
  
  print(paste("x:",x))
  print(paste("q:",q))
  # e^(-x)I_v(x), so take log term
  
  if (x/q < 0.3) {
    # scaled_I_q <- (1 / sqrt(2*pi*q)) * (exp(1) * x / (2*q))^q
    scaled_I_q <- besselI(x, nu = q, expon.scaled = TRUE)
  } else {
    scaled_I_q <- besselI(x, nu = q, expon.scaled = TRUE)
  }
  
  # log(I_q)
  ln_I_q <- -1*x + log(scaled_I_q)
  
  # ln(d_Heston) = log(C) - u - v + (q / 2) * log(v / u) + ln_I_q
  ln_transition_density <- log(C) - u - v + (q / 2) * log(v / u) + ln_I_q
  
  return(ln_transition_density)
}


# --- 2. 定义参数网格和固定值 (V_T 固定, V_t 变化) ---
V_T_fixed   <- 0.03
V_t_grid    <- seq(0.0001, 0.06, length.out = 100) 
k_fixed     <- 1.0 


Kappa_vals <- c(0.2, 0.5, 1, 2, 5, 10) 
Theta_vals <- c(0.01, 0.02, 0.03, 0.04, 0.05)
Sigma_vals <- c(0.1, 0.3, 0.5, 1.0, 1.5)

kappa_base <- 5.0
theta_base <- 0.03
sigma_base <- 0.5


results_kappa <- data.frame(V_t_1 = V_t_grid)
results_theta <- data.frame(V_t_1 = V_t_grid)
results_sigma <- data.frame(V_t_1 = V_t_grid)



# calculate density
for (k in Kappa_vals) {
  
  log_dens <- sapply(V_t_grid, function(vt_1) 
    ln_d_Heston(V_t = V_T_fixed, V_t_plus_k = vt_1, k = k_fixed, 
                kappa = k, theta = theta_base, sigma = sigma_base))
  results_kappa[[paste0("k_", k)]] <- log_dens
}

for (t in Theta_vals) {
  log_dens <- sapply(V_t_grid, function(vt_1) 
    ln_d_Heston(V_t = V_T_fixed, V_t_plus_k = vt_1, k = k_fixed, 
                kappa = kappa_base, theta = t, sigma = sigma_base))
  results_theta[[paste0("t_", t)]] <- log_dens
}


for (s in Sigma_vals) {
  log_dens <- sapply(V_t_grid, function(vt_1) 
    ln_d_Heston(V_t = V_T_fixed, V_t_plus_k = vt_1, k = k_fixed, 
                kappa = kappa_base, theta = theta_base, sigma = s))
  results_sigma[[paste0("s_", s)]] <- log_dens
}



# --- 1. Kappa Effect Plot ---
df_kappa_long <- pivot_longer(results_kappa, -V_t_1, names_to = "Kappa", values_to = "LogDensity")

plot_kappa <- ggplot(df_kappa_long, aes(x = V_t_1, y = LogDensity, color = Kappa, group = Kappa)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression(paste("Effect of Kappa (", kappa, ") on Inverse Transition Density (V"[T] * " = 0.03)")),
    x = expression(paste("Initial Volatility (", V[t], ")")),
    y = expression(paste("Log Density ", log[e] * (p(V[T] * " = 0.03 | " * V[t])))),
    subtitle = bquote(theta == .(theta_base) ~ ", " ~ sigma == .(sigma_base))
  ) +
  theme_minimal() +
  scale_color_discrete(labels = paste("kappa =", Kappa_vals)) +
  geom_vline(xintercept = V_T_fixed, linetype = "dotted", color = "black")


# --- 2. Theta Effect Plot ---
df_theta_long <- pivot_longer(results_theta, -V_t_1, names_to = "Theta", values_to = "LogDensity")

plot_theta <- ggplot(df_theta_long, aes(x = V_t_1, y = LogDensity, color = Theta, group = Theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression(paste("Effect of Theta (", theta, ") on Inverse Transition Density (V"[T] * " = 0.03)")),
    x = expression(paste("Initial Volatility (", V[t], ")")),
    y = expression(paste("Log Density ", log[e] * (p(V[T] * " = 0.03 | " * V[t])))),
    subtitle = bquote(kappa == .(kappa_base) ~ ", " ~ sigma == .(sigma_base))
  ) +
  theme_minimal() +
  scale_color_discrete(labels = paste("theta =", Theta_vals)) +
  geom_vline(xintercept = V_T_fixed, linetype = "dotted", color = "black")


# --- 3. Sigma Effect Plot ---
df_sigma_long <- pivot_longer(results_sigma, -V_t_1, names_to = "Sigma", values_to = "LogDensity")

plot_sigma <- ggplot(df_sigma_long, aes(x = V_t_1, y = LogDensity, color = Sigma, group = Sigma)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression(paste("Effect of Sigma (", sigma, ") on Inverse Transition Density (V"[T] * " = 0.03)")),
    x = expression(paste("Initial Volatility (", V[t], ")")),
    y = expression(paste("Log Density ", log[e] * (p(V[T] * " = 0.03 | " * V[t])))),
    subtitle = bquote(kappa == .(kappa_base) ~ ", " ~ theta == .(theta_base))
  ) +
  theme_minimal() +
  scale_color_discrete(labels = paste("sigma =", Sigma_vals)) +
  geom_vline(xintercept = V_T_fixed, linetype = "dotted", color = "black")

# Print all charts
print(plot_kappa)
print(plot_theta)
print(plot_sigma)


















V_T_fixed    <- 0.03
k_fixed      <- 1.0  
kappa_base   <- 2.0 
theta_base   <- 0.03 
sigma_base   <- 0.5

V_t_grid_3d  <- seq(0.005, 0.055, length.out = 30) 



Theta_grid_3d <- seq(0.005, 0.055, length.out = 30)

LogDensity_matrix_theta <- matrix(NA, nrow = length(V_t_grid_3d), ncol = length(Theta_grid_3d))

for (i in 1:length(V_t_grid_3d)) {
  for (j in 1:length(Theta_grid_3d)) {
    vt_1 <- V_t_grid_3d[i]
    theta <- Theta_grid_3d[j]
    
    LogDensity_matrix_theta[i, j] <- ln_d_Heston(
      V_t = V_T_fixed, V_t_plus_k = vt_1, k = k_fixed,
      kappa = kappa_base, theta = theta, sigma = sigma_base
    )
  }
}

plot_theta_3d <- plot_ly(
  x = V_t_grid_3d, y = Theta_grid_3d, z = LogDensity_matrix_theta, type = "surface",
  colorscale = "Viridis" 
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Initial Volatility V[t]"),
      yaxis = list(title = "Long-Run Mean Theta (θ)"),
      zaxis = list(title = "Log Density")
    ),
    title = paste0("Joint Effect of V[t] and Theta on Log Density (V[T]=", V_T_fixed, ")")
  )

print(plot_theta_3d)





Kappa_grid_3d <- seq(0.1, 5.0, length.out = 30) 

LogDensity_matrix_kappa <- matrix(NA, nrow = length(V_t_grid_3d), ncol = length(Kappa_grid_3d))

for (i in 1:length(V_t_grid_3d)) {
  for (j in 1:length(Kappa_grid_3d)) {
    vt_1 <- V_t_grid_3d[i]
    kappa <- Kappa_grid_3d[j]
    
    LogDensity_matrix_kappa[i, j] <- ln_d_Heston(
      V_t = V_T_fixed, V_t_plus_k = vt_1, k = k_fixed,
      kappa = kappa, theta = theta_base, sigma = sigma_base
    )
  }
}

plot_kappa_3d <- plot_ly(
  x = V_t_grid_3d, y = Kappa_grid_3d, z = LogDensity_matrix_kappa, type = "surface",
  colorscale = "Plasma"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Initial Volatility V[t]"),
      yaxis = list(title = "Mean Reversion Rate Kappa (κ)"),
      zaxis = list(title = "Log Density")
    ),
    title = paste0("Joint Effect of V[t] and Kappa on Log Density (V[T]=", V_T_fixed, ")")
  )

print(plot_kappa_3d)



Sigma_grid_3d <- seq(0.1, 1.5, length.out = 30) 

LogDensity_matrix_sigma <- matrix(NA, nrow = length(V_t_grid_3d), ncol = length(Sigma_grid_3d))

for (i in 1:length(V_t_grid_3d)) {
  for (j in 1:length(Sigma_grid_3d)) {
    vt_1 <- V_t_grid_3d[i]
    sigma <- Sigma_grid_3d[j]
    
    LogDensity_matrix_sigma[i, j] <- ln_d_Heston(
      V_t = V_T_fixed, V_t_plus_k = vt_1, k = k_fixed,
      kappa = kappa_base, theta = theta_base, sigma = sigma
    )
  }
}

plot_sigma_3d <- plot_ly(
  x = V_t_grid_3d, y = Sigma_grid_3d, z = LogDensity_matrix_sigma, type = "surface",
  colorscale = "Cividis"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Initial Volatility V[t]"),
      yaxis = list(title = "Vol of Vol Sigma (σ)"),
      zaxis = list(title = "Log Density")
    ),
    title = paste0("Joint Effect of V[t] and Sigma on Log Density (V[T]=", V_T_fixed, ")")
  )

print(plot_sigma_3d)
