source("Heston_likelihood.R")
########################

library(ggplot2)
library(reshape2)
library(dplyr)

kappa_base <- c(2, 1)
theta_base <- c(0.1, 0.2)
sigma_base <- c(0.1, 0.1)


best_params <- expand.grid(
  kappa = c(2, 1),
  theta = c(0.1, 0.2),
  sigma = c(0.05, 0.05)
)


# 固定的 Gamma 矩阵 (用于 LL 计算)
Gamma_fix <- matrix(c(0.99, 0.01, 0.09, 0.91), nrow=2, byrow=TRUE)




theta_fix <- theta_base
sigma_fix <- sigma_base
kappa_fix <- kappa_base

# --- VARYING RANGE ---
kappa_range <- seq(-1, 5, length.out = 30)

ll_matrix_kappa <- matrix(NA_real_, length(kappa_range), length(kappa_range)) 


for (i in 1:length(kappa_range)) {
  for (j in 1:length(kappa_range)) {
    kappa1 <- kappa_range[i]
    kappa2 <- kappa_range[j]
    
    params_test <- list(
      kappa = c(kappa1, kappa2),
      theta = theta_fix,
      sigma = sigma_fix,
      Gamma = Gamma_fix
    )
    
    ll_matrix_kappa[i, j] <- Heston_loglik_natural(params_test, V_simulated)
  }
}

# --- PLOT ---
ll_df_kappa <- melt(ll_matrix_kappa)
names(ll_df_kappa) <- c("Kappa1_Index", "Kappa2_Index", "LogLikelihood")

kappa_plot_full_grid <- ll_df_kappa %>%
  mutate(
    Kappa1 = kappa_range[Kappa1_Index],
    Kappa2 = kappa_range[Kappa2_Index]
  ) %>%
  
  ggplot(aes(x = Kappa1, y = Kappa2, z = LogLikelihood)) +
  
  geom_tile(aes(fill = LogLikelihood)) +
  
  stat_contour(
    color = "white",
    alpha = 0.8,
    linewidth = 0.5,
    bins = 50
  ) +
  
  geom_point(
    x = best_params$kappa[1],
    y = best_params$kappa[2],
    fill = "red", 
    color = "white", 
    size = 4,
    shape = 21 
  ) +
  scale_fill_gradientn(
    colors = c("blue", "cyan", "green", "yellow", "red"),
    name = "LL"
  ) +
  labs(
    title = expression("Log-Likelihood Surface: " * kappa[1] * " vs " * kappa[2] * " (Full Grid)"),
    subtitle = paste0(
      "Fixed: ",
      "theta=(", round(theta_fix[1], 3), ", ", round(theta_fix[2], 3), "), ",
      "sigma=(", round(sigma_fix[1], 3), ", ", round(sigma_fix[2], 3), "), ",
      "p11=", round(Gamma_fix[1,1], 3), ", p22=", round(Gamma_fix[2,2], 3)
    ),
    x = expression(kappa[1] * " (Regime 1)"),
    y = expression(kappa[2] * " (Regime 2)")
  ) +
  theme_minimal() +
  coord_fixed() 

print(kappa_plot_full_grid)

















theta_range <- seq(0.00, 0.5, length.out = 30)
ll_matrix_theta <- matrix(NA_real_, length(theta_range), length(theta_range))

# --- LOOP ---
for (i in 1:length(theta_range)) {
  for (j in 1:length(theta_range)) {
    theta1 <- theta_range[i]
    theta2 <- theta_range[j]
    
    params_test <- list(
      kappa = kappa_fix, 
      theta = c(theta1, theta2), 
      sigma = sigma_fix,
      Gamma = Gamma_fix
    )
    ll_matrix_theta[i, j] <- Heston_loglik_natural(params_test, V_simulated)
  }
}

# --- PLOT ---
ll_df_theta <- melt(ll_matrix_theta)
names(ll_df_theta) <- c("Theta1_Index", "Theta2_Index", "LogLikelihood")

theta_plot_full_grid <- ll_df_theta %>%
  mutate(
    Theta1 = theta_range[Theta1_Index],
    Theta2 = theta_range[Theta2_Index]
  ) %>%
  ggplot(aes(x = Theta1, y = Theta2, z = LogLikelihood)) +
  geom_tile(aes(fill = LogLikelihood)) +
  stat_contour(
    color = "white", 
    alpha = 0.8, 
    linewidth = 0.5, 
    bins = 50
  ) +
  geom_point(
    x = best_params$theta[1], 
    y = best_params$theta[2], 
    fill = "red", 
    color = "white", 
    size = 4, 
    shape = 21 
  ) +
  scale_fill_gradientn(
    colors = c("blue", "cyan", "green", "yellow", "red"),
    name = "LL"
  ) +
  labs(
    title = expression("Log-Likelihood Surface: " * theta[1] * " vs " * theta[2] * " (Full Grid)"),
    subtitle = paste0(
      "Fixed: ",
      "kappa=(", round(kappa_fix[1], 2), ", ", round(kappa_fix[2], 2), "), ",
      "sigma=(", round(sigma_fix[1], 3), ", ", round(sigma_fix[2], 3), "), ",
      "p11=", round(Gamma_fix[1,1], 3), ", p22=", round(Gamma_fix[2,2], 3)
    ),
    x = expression(theta[1] * " (Regime 1)"),
    y = expression(theta[2] * " (Regime 2)")
  ) +
  theme_minimal() +
  coord_fixed()

print(theta_plot_full_grid)