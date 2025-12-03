

source("Heston_likelihood.R")
Heston_loglik_natural <- function(params_list, observations) {
  
  kappa <- params_list$kappa
  theta <- params_list$theta
  sigma <- params_list$sigma
  Gamma <- params_list$Gamma
  
  nstates <- length(kappa)
  T_obs <- length(observations)
  
  delta <- try({
    solve(t(diag(nstates) - Gamma + 1), rep(1, nstates))
  }, silent = TRUE)
  
  if (inherits(delta, "try-error")) {
    delta <- rep(1, nstates) / nstates
  }
  
  
  allprobs_log <- matrix(NA_real_, nstates, T_obs - 1)
  
  for (i in 1:nstates) {
    # 确保调用的是返回 LOG 值的函数
    # 注意：这里假设您已经有了之前定义的 get_transition_density_heston_ln
    # 且该函数返回的是 log 值
    allprobs_log[i, ] <- get_transition_density_heston_ln(
      observations, kappa[i], theta[i], sigma[i]
    )
  }
  
  # 3. 计算 HMM Log-Likelihood
  # 使用之前修正过的 LL_HMM_Log_R (接受 Log 输入)
  ll <- LL_HMM_R(allprobs_log, Gamma, delta)
  
  # 处理数值异常
  if (is.na(ll) || is.infinite(ll)) {
    return(-1e300) # 返回极小值 (因为这是最大化任务)
  }
  
  return(ll)
}



library(foreach)
library(doParallel)

Heston_fit_grid <- function(observations, grid_params, ncluster = 1) {
  
  message("--- Starting Grid Search for Heston HMM ---")
  
  
  regime_grid <- expand.grid(
    kappa = grid_params$kappa,
    theta = grid_params$theta,
    sigma = grid_params$sigma
  )
  

  feller_check <- (2 * regime_grid$kappa * regime_grid$theta) > (regime_grid$sigma^2)

  
  n_regime_combs <- nrow(regime_grid)
  message(paste("Unique regime parameter combinations:", n_regime_combs))
  
  # 3. 构建双状态模型网格 (2-State Model Grid)
  # 我们需要组合 Regime 1 和 Regime 2，以及转移概率
  # 为了减少计算量，我们强制 theta_2 >= theta_1 (假设状态2是高均值状态)
  
  model_grid <- list()
  counter <- 0
  
  diag_probs <- grid_params$diag_prob
  
  for (i in 1:n_regime_combs) {
    for (j in 1:n_regime_combs) {
      
      # 对称性破缺：只计算 theta[j] >= theta[i] 的情况
      # 这样可以减少约 50% 的计算量，避免 (State A, State B) 和 (State B, State A) 重复计算
      if (regime_grid$theta[j] < regime_grid$theta[i]) next
      
      # 如果 theta 相同，则要求 sigma[j] >= sigma[i] 以进一步去重
      if (regime_grid$theta[j] == regime_grid$theta[i] && 
          regime_grid$sigma[j] < regime_grid$sigma[i]) next
      
      # 遍历转移概率
      for (p11 in diag_probs) {
        for (p22 in diag_probs) {
          
          counter <- counter + 1
          
          # 构建 Gamma 矩阵
          Gamma <- matrix(c(p11, 1-p11, 
                            1-p22, p22), nrow=2, byrow=TRUE)
          
          model_grid[[counter]] <- list(
            kappa = c(regime_grid$kappa[i], regime_grid$kappa[j]),
            theta = c(regime_grid$theta[i], regime_grid$theta[j]),
            sigma = c(regime_grid$sigma[i], regime_grid$sigma[j]),
            Gamma = Gamma
          )
        }
      }
    }
  }
  
  message(paste("Total models to evaluate:", length(model_grid)))

    results <- list()
    pb <- txtProgressBar(min = 0, max = length(model_grid), style = 3)
    for (idx in 1:length(model_grid)) {
      params <- model_grid[[idx]]
      ll <- Heston_loglik_natural(params, observations)
      results[[idx]] <- c(ll, idx)
      setTxtProgressBar(pb, idx)
    }
    close(pb)
  
  
  # 5. 提取最佳结果
  # results 是一个 list of vectors
  res_mat <- do.call(rbind, results)
  best_idx_in_res <- which.max(res_mat[, 1])
  best_ll <- res_mat[best_idx_in_res, 1]
  best_grid_idx <- res_mat[best_idx_in_res, 2]
  
  best_params <- model_grid[[best_grid_idx]]
  
  message(paste("--- Grid Search Completed ---"))
  message(paste("Best Log-Likelihood:", round(best_ll, 4)))
  
  return(list(
    best_params = best_params,
    max_ll = best_ll,
    all_results = res_mat 
  ))
}





library(foreach)
library(doParallel)
library(doSNOW)
library(progress) 

#' @title Grid Search for Heston HMM Parameters (Denser Grid)
Heston_fit_grid_dense <- function(observations, grid_defs, ncluster = 1) {
  
  message("--- Starting Denser Grid Search for Heston HMM ---")
  
  
  required_params <- c("kappa", "theta", "sigma", "diag_prob")
  if (!all(required_params %in% names(grid_defs))) {
    stop("grid_defs 必须包含 'kappa', 'theta', 'sigma', 'diag_prob' 的定义。")
  }
  
  kappa_vec <- seq(grid_defs$kappa$min, grid_defs$kappa$max, length.out = grid_defs$kappa$n_points)
  theta_vec <- seq(grid_defs$theta$min, grid_defs$theta$max, length.out = grid_defs$theta$n_points)
  sigma_vec <- seq(grid_defs$sigma$min, grid_defs$sigma$max, length.out = grid_defs$sigma$n_points)
  diag_probs <- seq(grid_defs$diag_prob$min, grid_defs$diag_prob$max, length.out = grid_defs$diag_prob$n_points)
  
  
  regime_grid <- expand.grid(
    kappa = kappa_vec,
    theta = theta_vec,
    sigma = sigma_vec
  )
  
  n_regime_combs <- nrow(regime_grid)
  
  model_grid <- list()
  counter <- 0
  
  for (i in 1:n_regime_combs) {
    for (j in 1:n_regime_combs) {
      if (regime_grid$theta[j] < regime_grid$theta[i]) next
      if (regime_grid$theta[j] == regime_grid$theta[i] && 
          regime_grid$sigma[j] < regime_grid$sigma[i]) next
      
      for (p11 in diag_probs) {
        for (p22 in diag_probs) {
          counter <- counter + 1
          Gamma <- matrix(c(p11, 1-p11, 1-p22, p22), nrow=2, byrow=TRUE)
          model_grid[[counter]] <- list(
            kappa = c(regime_grid$kappa[i], regime_grid$kappa[j]),
            theta = c(regime_grid$theta[i], regime_grid$theta[j]),
            sigma = c(regime_grid$sigma[i], regime_grid$sigma[j]),
            Gamma = Gamma
          )
        }
      }
    }
  }
  
  total_models <- length(model_grid)
  message(paste("Total models to evaluate:", total_models))
  

  pb <- progress::progress_bar$new(
    format = "[:bar] :percent [:elapsedfull / :eta] | Model: :current",
    total = total_models, 
    width = 60
  )
  
  start_time <- Sys.time()
  
  if (ncluster > 1) {

    cl <- parallel::makeCluster(ncluster)

    doSNOW::registerDoSNOW(cl)
    
    message(paste("Parallel processing started on", ncluster, "cores..."))
    

    opts <- list(progress = function(n) pb$tick())
    
   
    results <- foreach::foreach(idx = 1:total_models, 
                                .packages = c('base'), 
                                .export = c("Heston_loglik_natural", 
                                            "LL_HMM_R", 
                                            "log_sum_exp", 
                                            "get_transition_density_heston_ln", 
                                            "ln_d_Heston"),
                                .options.snow = opts) %dopar% {
                                  
                                  params <- model_grid[[idx]]
                                  ll <- Heston_loglik_natural(params, observations)
                                  return(c(ll, idx))
                                }
    parallel::stopCluster(cl)
    
  } else {

    message("Serial processing started...")
    results <- list()
    for (idx in 1:total_models) {
      params <- model_grid[[idx]]
      ll <- Heston_loglik_natural(params, observations)
      results[[idx]] <- c(ll, idx)

      pb$tick(tokens = list(current = idx)) 
    }
  }
  
  end_time <- Sys.time()
  
  res_mat <- do.call(rbind, results)

  best_idx_in_res <- which.max(res_mat[, 1]) 
  best_ll <- res_mat[best_idx_in_res, 1]
  best_grid_idx <- res_mat[best_idx_in_res, 2]
  best_params <- model_grid[[best_grid_idx]]
  
  message(paste("\n--- Grid Search Completed ---"))
  message(paste("Time taken:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
  message(paste("Best Log-Likelihood:", round(best_ll, 4)))
  
  return(list(
    best_params = best_params,
    max_ll = best_ll,
    all_results = res_mat
  ))
}








source("Model_Simulation.R")

set.seed(999)
N <- 250
v0 <- 0.03
S0 <- 100

Reg_chain <- simulate_Reg(series_length = N)




Reg_param_theta <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta,  sigma,  rho
    0.5,   2,    0.1,   0.05,   -0.1, # calm
    0.5,   1,     0.2,    0.05,   -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)


# 62-144 regime change
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param_theta, T, N, M=1, method = "M")
S_simulated <- sim_series$S_paths
plot(S_simulated,  type = "l")
# S_mu_model <- fit_HMM(S_simulated, Reg_chain)
# summary(S_mu_model)
# S_mu_model


V_simulated <- sim_series$V_paths



# my_grids <- list(
#   kappa = c(1, 10),         
#   theta = c(0.02, 0.8),   
#   sigma = c(0.05, 0.2),     
#   diag_prob = c(0.99, 0.99) 
# )
# 
# 
# 
# grid_result <- Heston_fit_grid(V_simulated, my_grids, ncluster = 4)
# 
# print(grid_result$best_params)














Reg_param_theta <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta,  sigma,  rho
    0.5,   10,    0.03,   0.05,   -0.1, # calm
    0.5,   5,     0.05,    0.05,   -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)




dense_grids <- list(
  kappa = list(min = 0, max = 15, n_points = 15),
  theta = list(min = 0.00, max = 0.6, n_points =13),
  sigma = list(min = 0.03, max = 0.07, n_points = 5),

  diag_prob = list(min = 0.99, max = 0.99, n_points = 1)
)

grid_result <- Heston_fit_grid_dense(V_simulated, dense_grids, ncluster = 14)

print(grid_result$best_params)



best_params <- grid_result$best_params

kappa_fix <- best_params$kappa
sigma_fix <- best_params$sigma
Gamma_fix <- best_params$Gamma # 包含 p11 和 p22

library(ggplot2)
library(dplyr)
library(reshape2)


theta_range <- seq(0.00, 0.5, length.out = 30) # 30x30 = 900个点
theta_range <- seq(0.00, 0.5, length.out = 30) # 30x30 = 900个点

ll_matrix <- matrix(NA_real_, length(theta_range), length(theta_range))

for (i in 1:length(theta_range)) {
  for (j in 1:length(theta_range)) {
    theta1 <- theta_range[i]
    theta2 <- theta_range[j]
    
    # 保持其他参数固定
    params_test <- list(
      kappa = kappa_fix, 
      theta = c(theta1, theta2), 
      sigma = sigma_fix,
      Gamma = Gamma_fix
    )
    
    
    if (theta2 >= theta1) {
      ll_matrix[i, j] <- Heston_loglik_natural(params_test, V_simulated)
    } else {
      ll_matrix[i, j] <- NA_real_ # 不符合条件不计算
    }
  }
}





# 将矩阵转换为长格式数据框 (Long Format Data Frame)
ll_df <- melt(ll_matrix)
names(ll_df) <- c("Theta1_Index", "Theta2_Index", "LogLikelihood")

# 将索引映射回实际的 Theta 值
ll_df <- ll_df %>%
  mutate(
    Theta1 = theta_range[Theta1_Index],
    Theta2 = theta_range[Theta2_Index]
  ) %>%
  filter(!is.na(LogLikelihood)) # 移除未计算的点

# 绘制等高线图
contour_plot <- ggplot(ll_df, aes(x = Theta1, y = Theta2, z = LogLikelihood)) +
  # 1. 热力图层 (Heatmap)
  geom_tile(aes(fill = LogLikelihood)) +
  # 2. 等高线层 (Contour)
  stat_contour(
    color = "white", 
    alpha = 0.8,
    linewidth = 0.5,
    bins = 15 # 更多的线可以显示更多细节
  ) +
  # 3. 最佳参数点 (Grid Search 的结果)
  geom_point(
    x = best_params$theta[1], 
    y = best_params$theta[2], 
    color = "red", 
    size = 3, 
    shape = 21
  ) +
  scale_fill_gradientn(
    colors = c("blue", "cyan", "green", "yellow", "red"),
    name = "Log-Likelihood"
  ) +
  labs(
    title = expression("Heston HMM Log-Likelihood Surface: " * theta[1] * " vs " * theta[2]),
    subtitle = paste0(
      "Fixed: ",
      "kappa=(", round(kappa_fix[1], 2), ", ", round(kappa_fix[2], 2), "), ",
      "sigma=(", round(sigma_fix[1], 2), ", ", round(sigma_fix[2], 2), "), ",
      "p11=", round(Gamma_fix[1,1], 3), ", p22=", round(Gamma_fix[2,2], 3)
    ),
    x = expression(theta[1] * " (Regime 1)"),
    y = expression(theta[2] * " (Regime 2)")
  ) +
  theme_minimal() +
  coord_fixed() # 保持 X/Y 轴比例

print(contour_plot)















# ----------------------------------------------------
# 初始固定参数和依赖项 (请确保它们在您的环境中已定义)
# ----------------------------------------------------

# 初始固定参数组合 (基于您的要求)
kappa_base <- c(10, 5)
theta_base <- c(0.03, 0.6)
sigma_base <- c(0.05, 0.05)


best_params <- expand.grid(
  kappa = c(10, 5),
  theta = c(0.03, 0.6),
  sigma = c(0.05, 0.05)
)


# 固定的 Gamma 矩阵 (用于 LL 计算)
Gamma_fix <- matrix(c(0.99, 0.01, 0.09, 0.91), nrow=2, byrow=TRUE)

# 假设 V_simulated 和 best_params 已经定义
library(ggplot2)
library(reshape2)
library(dplyr) 





























########################

library(ggplot2)
library(reshape2)
library(dplyr)


theta_fix <- theta_base
sigma_fix <- sigma_base

# --- VARYING RANGE ---
kappa_range <- seq(-20, 20, length.out = 60)

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
    bins = 15
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




















# --- FIXED PARAMETERS ---
kappa_fix <- c(10, 5)
sigma_fix <- c(0.05, 0.05)
# 假设 Gamma_fix 已定义，例如：
# Gamma_fix <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2, byrow=TRUE)

# --- VARYING RANGE ---
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
    bins = 15
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