library(ggplot2)
library(tidyr)
library(dplyr)


forward_filter <- function(observations, nstates, Gamma, kappa, theta, sigma) {
  
  T <- length(observations) - 1
  logGamma <- log(Gamma)
  delta <- log(oeli::stationary_distribution(Gamma))
  
  # emission log-prob
  allprobs <- matrix(0, nstates, T)
  for (i in 1:nstates) {
    allprobs[i,] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
  }
  
  logalpha <- matrix(-Inf, nstates, T)
  
  # initialization
  logalpha[,1] <- delta + allprobs[,1]
  
  # recursion
  for (t in 2:T) {
    for (j in 1:nstates) {
      logalpha[j,t] <- logsumexp(logalpha[,t-1] + logGamma[,j]) + allprobs[j,t]
    }
  }
  
  # normalized posterior: p(s_t=j | y_1:t)
  posterior <- matrix(0, nstates, T)
  reg <- matrix(0, 1, T)
  for (t in 1:T) {
    c <- logsumexp(logalpha[,t])
    posterior[,t] <- exp(logalpha[,t] - c)
    reg[,t] <- which.max(posterior[,t])
  }
  
  return(list(posterior = posterior,
              states_estimate = as.vector(reg)))
}

logsumexp <- function(x) {
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
  
}

plot_viterbi <- function(V_simulated, nstates, Gamma, kappa, theta, sigma, Reg_chain){
  
  
  prob <- forward_filter(V_simulated, nstates, Gamma, kappa, theta, sigma)
  states_estimate <- prob$states_estimate
  a <- prob$posterior

  time_index <- 1:ncol(a) # 假设时间索引从 1 到 N
  
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = a[1,],
    Prob_State2 = 1 - a[1,]
  )
 
  
  state_data <- data.frame(
    Time = time_index,
    Estimated_State = states_estimate -1, 
    True_State = as.factor(Reg_chain[1:ncol(a)])        
  )
  
  
print("np")
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  


  max_time <- max(state_data$Time) 
  regime_start_times <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>% 
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>% 
    rename(Time_Start = Time)

  
plot <- ggplot() +
    
    geom_rect(
      data = regime_start_times,
      aes(
        xmin = Time_Start, 
        xmax = Time_End, 
        ymin = -0.05, 
        ymax = 1.05,
        fill = True_State # True_State 映射到背景色图例
      ),
      alpha = 0.2, 
      inherit.aes = FALSE
    ) +
    
  
    geom_line(
      data = prob_data_long, 
      aes(x = Time, y = Probability, color = Probability_Type), # Probability_Type 映射到颜色图例
      linewidth = 0.9
    ) +
    
    # III. 机制估计点 (现在使用 shape 映射来创建图例)
    geom_point(
      data = state_data,
      # 映射 shape = "Estimated State" 以创建图例
      aes(x = Time, y = Estimated_State, shape = "Estimated State"), 
      size = 1,
      color = "black", # 强制颜色为黑色
      alpha = 0.8
    ) +
    

    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"), 
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)"),
      name = "True Regime Background"
    ) +
  scale_color_manual(
    values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
    labels = c("Prob_State1" = expression(P(X[1:t], S[t] == 1 * " | " * bold(theta))),
               "Prob_State2" = expression(P(X[1:t], S[t] == 2 * " | " * bold(theta)))),
               name = "Forward Probability" 
    )+
    scale_shape_manual(
      # 定义 'Estimated State' 对应的形状为 'x' (pch=4)
      values = c("Estimated State" = 1), 
      name = "Regime Estimate" # 新的图例标题
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      name = "Posterior Probability / State"
    ) +
    labs(
      title = "HMM Regime Switching Heston: Estimation",
      x = "Time Step",
    ) +
    theme_minimal() +
    # 统一图例顺序
    guides(
      fill = guide_legend(order = 1),       
      color = guide_legend(order = 2),
      shape = guide_legend(order = 3) # 确保 shape 图例显示
    )
  return(plot)
}








# 简化版：只需要真实波动率、参数和状态来构建置信区间
# plot_cir_confidence_simple <- function(true_vol, 
#                                        param, 
#                                        states_estimate, 
#                                        dt = 1/252, 
#                                        interval_step = 10,
#                                        V_grid = NULL) {
#   
#   # 确保长度一致
#   min_length <- min(length(true_vol), length(states_estimate))
#   true_vol <- true_vol[1:min_length]
#   states_estimate <- states_estimate[1:min_length]
#   time_points <- 1:min_length
#   
#   cat("数据长度:", min_length, "\n")
#   
#   # 自动生成V_grid
#   if (is.null(V_grid)) {
#     V_min <- max(0.0001, min(true_vol) * 0.1)
#     V_max <- max(true_vol) * 3
#     V_grid <- seq(V_min, V_max, length.out = 500)
#   }
#   
#   # 找到regime转换点
#   changepoints <- c(1, which(diff(states_estimate) != 0) + 1, length(states_estimate) + 1)
#   
#   # 存储所有置信区间数据
#   all_ci_data <- data.frame()
#   
#   # 为每个regime段计算置信区间
#   for (i in 1:(length(changepoints) - 1)) {
#     start_idx <- changepoints[i]
#     end_idx <- changepoints[i + 1] - 1
#     
#     if (end_idx - start_idx < interval_step) next
#     
#     current_regime <- states_estimate[start_idx]
#     kappa_reg <- param$kappa[current_regime]
#     theta_reg <- param$theta[current_regime]
#     sigma_reg <- param$sigma[current_regime]
#     
#     # 使用regime的第一个点作为参考
#     v_ref <- true_vol[start_idx]
#     
#     # 在该regime内每隔interval_step计算一次置信区间
#     for (j in seq(interval_step, (end_idx - start_idx + 1), by = interval_step)) {
#       t_idx <- start_idx + j - 1
#       if (t_idx > end_idx) break
#       
#       k <- j * dt  # 从regime起点开始的时间间隔
#       
#       # 使用 ln_d_Heston 计算密度
#       log_densities <- sapply(V_grid, function(v) {
#         ln_d_Heston(V_t = v_ref, 
#                     V_t_plus_k = v, 
#                     k = k, 
#                     kappa = kappa_reg, 
#                     theta = theta_reg, 
#                     sigma = sigma_reg)
#       })
#       
#       # 转换为密度
#       densities <- exp(log_densities)
#       
#       # 检查有效性
#       if (any(is.na(densities)) || all(densities == 0)) {
#         warning(sprintf("Invalid densities at time %d, regime %d", t_idx, current_regime))
#         next
#       }
#       
#       # 数值积分计算CDF
#       dV <- diff(V_grid)[1]
#       cdf <- cumsum(densities * dV)
#       cdf <- cdf / max(cdf)  # 归一化
#       
#       # 找到5th和95th百分位数
#       idx_lower <- which.min(abs(cdf - 0.05))
#       idx_upper <- which.min(abs(cdf - 0.95))
#       
#       lower_val <- V_grid[idx_lower]
#       upper_val <- V_grid[idx_upper]
#       
#       # 验证合理性
#       if (!is.na(lower_val) && !is.na(upper_val) && 
#           lower_val > 0 && upper_val > lower_val) {
#         
#         new_row <- data.frame(
#           time = t_idx,
#           lower = lower_val,
#           upper = upper_val,
#           regime = current_regime
#         )
#         all_ci_data <- rbind(all_ci_data, new_row)
#       }
#     }
#   }
#   
#   if (nrow(all_ci_data) == 0) {
#     warning("No confidence intervals computed")
#     return(NULL)
#   }
#   
#   # 准备绘图数据
#   vol_data <- data.frame(
#     time = time_points, 
#     true = true_vol
#   )
#   
#   cat("\n置信区间样本:\n")
#   print(head(all_ci_data, 10))
#   
#   # 绘图
#   p <- ggplot() +
#     # 真实波动率线
#     geom_line(data = vol_data, 
#               aes(x = time, y = true), 
#               color = "red", size = 0.5, alpha = 0.8) +
#     
#     # 置信区间点
#     geom_point(data = all_ci_data, 
#                aes(x = time, y = lower, fill = factor(regime)), 
#                shape = 21, size = 0.7, alpha = 0.7, stroke = 1) +
#     geom_point(data = all_ci_data, 
#                aes(x = time, y = upper, fill = factor(regime)), 
#                shape = 21, size = 0.7, alpha = 0.7, stroke = 1) +
#     
#     # 连接上下界的竖线
#     geom_segment(data = all_ci_data,
#                  aes(x = time, xend = time, y = lower, yend = upper, 
#                      color = factor(regime)),
#                  alpha = 0.3, size = 0.8) +
#     
#     scale_color_manual(
#       values = c("1" = "blue", "2" = "green"),
#       name = "Regime"
#     ) +
#     scale_fill_manual(
#       values = c("1" = "blue", "2" = "green"),
#       name = "Regime (90% CI)"
#     ) +
#     labs(
#       title = "CIR Model Confidence Intervals (via ln_d_Heston)",
#       subtitle = sprintf("90%% CI | Every %d days | Red line: True volatility", interval_step),
#       x = "Time (Days)", 
#       y = "Volatility"
#     ) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(size = 14, face = "bold"),
#       legend.position = "bottom"
#     )
#  
#   
#   return(list(plot = p, ci_data = all_ci_data))
# }





plot_cir_confidence_simple <- function(true_vol, 
                                       param, 
                                       states_estimate, 
                                       dt = 1/252, 
                                       interval_step = 10,
                                       V_grid = NULL) {
  
  min_length <- min(length(true_vol), length(states_estimate))
  true_vol <- true_vol[1:min_length]
  states_estimate <- states_estimate[1:min_length]
  time_points <- 1:min_length
  
  
  if (is.null(V_grid)) {
    V_min <- max(0.0001, min(true_vol) * 0.1)
    V_max <- max(true_vol) * 3
    V_grid <- seq(V_min, V_max, length.out = 500)
  }
  
  
  changepoints <- c(1, which(diff(states_estimate) != 0) + 1, length(states_estimate) + 1)
  
  
  all_ci_data <- data.frame()
  
  
  for (i in 1:(length(changepoints) - 1)) {
    start_idx <- changepoints[i]
    end_idx <- changepoints[i + 1] - 1
    
    if (end_idx - start_idx < interval_step) next
    
    current_regime <- states_estimate[start_idx]
    kappa_reg <- param$kappa[current_regime]
    theta_reg <- param$theta[current_regime]
    sigma_reg <- param$sigma[current_regime]
    
    v_ref <- true_vol[start_idx]
    
    
    for (j in seq(interval_step, (end_idx - start_idx + 1), by = interval_step)) {
      t_idx <- start_idx + j - 1
      if (t_idx > end_idx) break
      
      k <- j * dt # time steps from the regime change point
      
      log_densities <- sapply(V_grid, function(v) {
        ln_d_Heston(V_t = v_ref, 
                    V_t_plus_k = v, 
                    k = k, 
                    kappa = kappa_reg, 
                    theta = theta_reg, 
                    sigma = sigma_reg)
      })
      
      densities <- exp(log_densities)
      
      if (any(is.na(densities)) || all(densities == 0)) {
        warning(sprintf("Invalid densities at time %d, regime %d", t_idx, current_regime))
        next
      }
      
      dV <- diff(V_grid)[1]
      cdf <- cumsum(densities * dV)
      cdf <- cdf / max(cdf) 
      
      
      idx_lower <- which.min(abs(cdf - 0.025))
      idx_mean <- which.min(abs(cdf - 0.5))
      idx_upper <- which.min(abs(cdf - 0.975))
      
      lower_val <- V_grid[idx_lower]
      mean_val <- V_grid[idx_mean]
      upper_val <- V_grid[idx_upper]
      
      # 验证合理性
      if (!is.na(lower_val) && !is.na(upper_val) && 
          lower_val > 0 && upper_val > lower_val) {
        
        new_row <- data.frame(
          time = t_idx,
          lower = lower_val,
          mean = mean_val,
          upper = upper_val,
          regime = current_regime
        )
        all_ci_data <- rbind(all_ci_data, new_row)
      }
    }
  }
  
  if (nrow(all_ci_data) == 0) {
    warning("No confidence intervals computed")
    return(NULL)
  }
  
  # 准备绘图数据
  vol_data <- data.frame(
    time = time_points, 
    true = true_vol
  )
  
  cat("\n置信区间样本:\n")
  print(head(all_ci_data, 10))
  
  # 绘图
  p <- ggplot() +
    # 真实波动率线
    geom_line(data = vol_data, 
              aes(x = time, y = true), 
              color = "black", size = 0.5, alpha = 0.8) +
    
    # 置信区间点
    geom_point(data = all_ci_data, 
               aes(x = time, y = lower, fill = factor(regime)), 
               shape = 21, size = 0.7, alpha = 0.7, stroke = 1) +
    geom_point(data = all_ci_data, 
               aes(x = time, y = mean, fill = factor(regime)), 
               shape = 21, size = 0.7, alpha = 0.7, stroke = 1)+
    geom_point(data = all_ci_data, 
               aes(x = time, y = upper, fill = factor(regime)), 
               shape = 21, size = 0.7, alpha = 0.7, stroke = 1) +
    
    # 连接上下界的竖线
    geom_segment(data = all_ci_data,
                 aes(x = time, xend = time, y = lower, yend = upper, 
                     color = factor(regime)),
                 alpha = 0.3, size = 0.8) +
    
    scale_color_manual(
      values = c("1" = "blue", "2" = "red"),
      name = "Regime"
    ) +
    scale_fill_manual(
      values = c("1" = "blue", "2" = "red"),
      name = "Regime (95% CI)"
    ) +
    labs(
      title = "CIR Model Confidence Intervals",
      subtitle = sprintf("95%% CI | Every %d days | black line: True volatility", interval_step),
      x = "Time (Days)", 
      y = "Volatility"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  
  return(list(plot = p, ci_data = all_ci_data))
}













# ================================================================
# Plot Viterbi Results with O-M Functional (Same Format as Original)
# ================================================================
plot_viterbi_om <- function(batch_results_complete,
                            nstates, 
                            Gamma, 
                            kappa, 
                            theta, 
                            sigma,
                            Reg_chain,
                            lambda_om = 1.0,
                            use_upper = TRUE,
                            use_lower = TRUE,
                            normalize_method = "log",
                            show_om_contribution = FALSE) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # ===== Step 1: Prepare O-M matrices =====
  cat("Preparing O-M matrices...\n")
  
  # Normalize O-M values if needed
  if (normalize_method == "log") {
    OM_upper <- log(batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] + 1)
    OM_lower <- log(batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] + 1)
  } else if (normalize_method == "zscore") {
    OM_all <- c(
      batch_results_complete$OM_upper_R1,
      batch_results_complete$OM_upper_R2,
      batch_results_complete$OM_lower_R1,
      batch_results_complete$OM_lower_R2
    )
    OM_finite <- OM_all[is.finite(OM_all)]
    mean_om <- mean(OM_finite)
    sd_om <- sd(OM_finite)
    
    OM_upper <- (batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] - mean_om) / sd_om
    OM_lower <- (batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] - mean_om) / sd_om
  } else {
    # Raw values
    OM_upper <- batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")]
    OM_lower <- batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")]
  }
  
  OM_upper <- as.matrix(OM_upper)
  OM_lower <- as.matrix(OM_lower)
  
  # ===== Step 2: Run Viterbi Algorithm (FULL: Forward + Backward) =====
  cat("Running Viterbi Algorithm (Global Decoding)...\n")
  
  viterbi_res <- viterbi_om_pure(
    OM_upper = OM_upper,
    OM_lower = OM_lower,
    nstates = nstates,
    Gamma = Gamma,
    lambda_om = lambda_om,
    use_upper = use_upper,
    use_lower = use_lower
  )
  
  states_estimate <- viterbi_res$states_estimate
  viterbi_score <- viterbi_res$viterbi_score  # [nstates x T]
  
  # ===== Step 2.5: Compute Delta (Normalized Path Probabilities) =====
  # Delta represents confidence in the Viterbi path at each time
  T_len <- ncol(viterbi_score)
  delta <- matrix(0, nstates, T_len)
  
  for (t in 1:T_len) {
    # Softmax normalization of Viterbi scores
    max_score <- max(viterbi_score[, t])
    exp_scores <- exp(viterbi_score[, t] - max_score)
    delta[, t] <- exp_scores / sum(exp_scores)
  }
  
  # ===== Step 3: Prepare data for plotting =====
  
  time_index <- 1:T_len
  
  # Delta probability data
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = delta[1, ],
    Prob_State2 = delta[2, ]
  )
  
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  # State data
  state_data <- data.frame(
    Time = time_index,
    Estimated_State = states_estimate - 1,
    True_State = as.factor(Reg_chain[1:T_len])
  )
  
  # Regime background rectangles
  max_time <- max(state_data$Time)
  regime_start_times <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  
  # ===== Step 4: Create main plot =====
  
  main_plot <- ggplot() +
    # Background shading for true regimes
    geom_rect(
      data = regime_start_times,
      aes(
        xmin = Time_Start,
        xmax = Time_End,
        ymin = -0.05,
        ymax = 1.05,
        fill = True_State
      ),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    
    # Delta probability lines (path confidence)
    geom_line(
      data = prob_data_long,
      aes(x = Time, y = Probability, color = Probability_Type),
      linewidth = 0.9
    ) +
    
    # Viterbi optimal path (step line instead of points)
    geom_step(
      data = state_data,
      aes(x = Time, y = Estimated_State, linetype = "Viterbi Path"),
      color = "black",
      linewidth = 1.2,
      alpha = 0.8
    ) +
    
    # Scales
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)"),
      name = "True Regime Background"
    ) +
    scale_color_manual(
      values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
      labels = c(
        "Prob_State1" = expression(delta[1](t) ~ "(R1 path score)"),
        "Prob_State2" = expression(delta[2](t) ~ "(R2 path score)")
      ),
      name = "Viterbi Delta (Path Confidence)"
    ) +
    scale_linetype_manual(
      values = c("Viterbi Path" = "solid"),
      name = "Optimal Path"
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      name = "Delta Probability / State"
    ) +
    labs(
      title = sprintf("Viterbi HMM with O-M Functional (norm=%s)", 
                      normalize_method),
      subtitle = sprintf("Global Optimal Path | Accuracy: %.2f%%", 
                         mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100),
      x = "Time Step"
    ) +
    theme_minimal() +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2),
      linetype = guide_legend(order = 3)
    )
  
  # ===== Step 5: Optional O-M contribution plot (PROBABILITY SCALE) =====
  
  if (show_om_contribution) {
    
    # Calculate O-M contribution as PROBABILITY (not log)
    # Emission probability: exp(-λ × L[φ])
    
    # Total O-M functional per regime
    L_total_R1 <- OM_upper[, 1] + OM_lower[, 1]
    L_total_R2 <- OM_upper[, 2] + OM_lower[, 2]
    
    # Convert to probability: P ∝ exp(-λ × L)
    prob_R1 <- exp(-lambda_om * L_total_R1)
    prob_R2 <- exp(-lambda_om * L_total_R2)
    
    # Normalize to make probabilities comparable
    # (Each regime's probability at each time point)
    total_prob <- prob_R1 + prob_R2
    prob_R1_norm <- prob_R1 / total_prob
    prob_R2_norm <- prob_R2 / total_prob
    
    om_contrib_data <- data.frame(
      Time = 1:nrow(batch_results_complete),
      Prob_R1 = prob_R1_norm,
      Prob_R2 = prob_R2_norm
    )
    
    om_contrib_long <- pivot_longer(
      om_contrib_data,
      cols = starts_with("Prob_"),
      names_to = "Regime",
      values_to = "OM_Probability"
    )
    
    om_plot <- ggplot(om_contrib_long, aes(x = Time, y = OM_Probability, color = Regime)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
      scale_color_manual(
        values = c("Prob_R1" = "blue", "Prob_R2" = "red"),
        labels = c("Prob_R1" = "Regime 1 (Calm)", "Prob_R2" = "Regime 2 (Turbulent)"),
        name = "Regime"
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      ) +
      labs(
        title = "O-M Based Path Probability (Normalized)",
        subtitle = expression(
          "P(path | regime) = " * exp(-lambda %.% L[phi])
        ),
        x = "Time Step",
        y = "Probability"
      ) +
      theme_minimal() +
      theme(
        plot.subtitle = element_text(size = 9, color = "gray30")
      )
    
    # Combine plots
    library(gridExtra)
    combined_plot <- grid.arrange(main_plot, om_plot, ncol = 1, heights = c(2, 1))
    
    return(list(
      plot = combined_plot,
      main_plot = main_plot,
      om_plot = om_plot,
      states_estimate = states_estimate,
      delta = delta,  # Instead of posterior
      viterbi_score = viterbi_score,
      om_probabilities = om_contrib_data,
      accuracy = mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State)))
    ))
    
  } else {
    return(list(
      plot = main_plot,
      states_estimate = states_estimate,
      delta = delta,  # Instead of posterior
      viterbi_score = viterbi_score,
      accuracy = mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State)))
    ))
  }
}
