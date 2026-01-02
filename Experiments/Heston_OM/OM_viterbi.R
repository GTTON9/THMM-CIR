

#' Forward Filter for Pure O-M Functional (THMM)
#' 
#' Following PDF Section 11.3.1 - Forward Probabilities
#' Formula: α_k(j) = [Σ_i α_{k-1}(i) P_ij] exp(-L_j[φ_k])
#'
#' @param OM_upper Matrix [T x nstates] of O-M functionals for upper paths
#' @param OM_lower Matrix [T x nstates] of O-M functionals for lower paths
#' @param nstates Number of hidden states
#' @param Gamma Transition probability matrix
#' @param lambda_om Weight for O-M contribution (default 1.0)
#' @param use_upper Whether to use upper path (default TRUE)
#' @param use_lower Whether to use lower path (default TRUE)
#'
#' @return List containing:
#'   - posterior: Matrix [nstates x T] of P(s_t = i | φ_{1:t})
#'   - alpha: Matrix [nstates x T] of forward probabilities
#'   - states_estimate: Vector [T] of most likely states at each time
#'   - L_total: Matrix [nstates x T] of total O-M costs

forward_filter_om_pure <- function(OM_upper,
                                   OM_lower,
                                   nstates,
                                   Gamma,
                                   lambda_om = 1.0,
                                   use_upper = TRUE,
                                   use_lower = TRUE) {
  
  T <- nrow(OM_upper)
  
  cat("=================================================================\n")
  cat("Pure O-M Forward Filter (THMM)\n")
  cat("=================================================================\n\n")
  
  # ===== Step 1: Compute Total O-M Functional for Each Path =====
  # This is L_i[φ_k] in PDF formula
  
  cat("Computing O-M functionals...\n")
  
  L_total <- matrix(0, T, nstates)
  
  for (t in 1:T) {
    for (i in 1:nstates) {
      
      om_value <- 0
      
      # Upper path
      if (use_upper && !is.na(OM_upper[t, i]) && is.finite(OM_upper[t, i])) {
        om_value <- om_value + lambda_om * OM_upper[t, i]
      }
      
      # Lower path
      if (use_lower && !is.na(OM_lower[t, i]) && is.finite(OM_lower[t, i])) {
        om_value <- om_value + lambda_om * OM_lower[t, i]
      }
      
      L_total[t, i] <- om_value
    }
  }
  
  cat("  L_total statistics:\n")
  cat("    Mean:", mean(L_total), "\n")
  cat("    Range: [", min(L_total), ",", max(L_total), "]\n\n")
  
  # ===== Step 2: Initial Distribution =====
  # π_i = stationary distribution
  
  delta <- oeli::stationary_distribution(Gamma)
  
  cat("Initial distribution (π):\n")
  for (i in 1:nstates) {
    cat(sprintf("  State %d: %.4f\n", i, delta[i]))
  }
  cat("\n")
  
  # ===== Step 3: Forward Algorithm (Log-Space) =====
  
  # log_alpha[i, t] = log P(φ_1,...,φ_t, s_t = i)
  log_alpha <- matrix(-Inf, nstates, T)
  
  # Initialization (k=1)
  # Formula: α_1(i) = π_i × exp(-L_i[φ_1])
  # Log form: log α_1(i) = log π_i - L_i[φ_1]
  
  for (i in 1:nstates) {
    log_alpha[i, 1] <- log(delta[i]) - L_total[1, i]
  }
  
  cat("Initialization:\n")
  cat("  log_alpha[, 1] =", log_alpha[, 1], "\n\n")
  
  # Recursion (k = 2, ..., T)
  # Formula: α_k(j) = [Σ_i α_{k-1}(i) P_ij] × exp(-L_j[φ_k])
  # Log form: log α_k(j) = log[Σ_i exp(log α_{k-1}(i) + log P_ij)] - L_j[φ_k]
  
  cat("Forward recursion...\n")
  pb <- txtProgressBar(min = 2, max = T, style = 3)
  
  log_Gamma <- log(Gamma)
  
  for (t in 2:T) {
    for (j in 1:nstates) {
      
      # Compute Σ_i α_{k-1}(i) P_ij in log space
      # log Σ_i exp(x_i) = logsumexp(x)
      
      log_terms <- numeric(nstates)
      for (i in 1:nstates) {
        log_terms[i] <- log_alpha[i, t-1] + log_Gamma[i, j]
      }
      
      # logsumexp for numerical stability
      log_sum <- logsumexp(log_terms)
      
      # Add emission probability: exp(-L_j[φ_k])
      log_alpha[j, t] <- log_sum - L_total[t, j]
    }
    
    setTxtProgressBar(pb, t)
  }
  
  close(pb)
  
  # ===== Step 4: Compute Posterior Probabilities =====
  # P(s_t = i | φ_{1:t}) = α_t(i) / Σ_j α_t(j)
  
  cat("\n\nComputing posterior probabilities...\n")
  
  posterior <- matrix(0, nstates, T)
  
  for (t in 1:T) {
    # Normalize: P(s_t = i | φ_{1:t}) ∝ α_t(i)
    log_norm <- logsumexp(log_alpha[, t])
    
    for (i in 1:nstates) {
      posterior[i, t] <- exp(log_alpha[i, t] - log_norm)
    }
  }
  
  # ===== Step 5: State Estimation =====
  # Most likely state at each time: argmax_i P(s_t = i | φ_{1:t})
  
  states_estimate <- apply(posterior, 2, which.max)
  
  # ===== Step 6: Summary Statistics =====
  
  transitions <- diff(states_estimate)
  n_switches <- sum(transitions != 0)
  
  cat("\n=== Forward Filter Results (Pure O-M) ===\n")
  cat("Total time windows:", T, "\n")
  cat("Number of regime switches:", n_switches, "\n")
  
  for (i in 1:nstates) {
    n_in_regime <- sum(states_estimate == i)
    pct <- 100 * n_in_regime / T
    avg_posterior <- mean(posterior[i, states_estimate == i], na.rm = TRUE)
    
    cat(sprintf("  Regime %d:\n", i))
    cat(sprintf("    Frequency: %d windows (%.1f%%)\n", n_in_regime, pct))
    cat(sprintf("    Avg posterior when assigned: %.4f\n", avg_posterior))
  }
  
  cat("\nDetected switches:\n")
  switch_idx <- which(transitions != 0)
  if (length(switch_idx) > 0) {
    for (idx in switch_idx) {
      cat(sprintf("  t=%d: R%d → R%d (prob: %.3f → %.3f)\n", 
                  idx, 
                  states_estimate[idx], 
                  states_estimate[idx+1],
                  posterior[states_estimate[idx], idx],
                  posterior[states_estimate[idx+1], idx+1]))
    }
  } else {
    cat("  (No regime switches detected)\n")
  }
  
  # Compute average posterior confidence
  avg_confidence <- mean(apply(posterior, 2, max))
  cat(sprintf("\nAverage posterior confidence: %.4f\n", avg_confidence))
  
  cat("=================================================================\n")
  
  return(list(
    posterior = posterior,
    alpha = exp(log_alpha),
    log_alpha = log_alpha,
    states_estimate = states_estimate,
    L_total = L_total,
    n_switches = n_switches,
    switch_points = switch_idx,
    avg_confidence = avg_confidence
  ))
}


#' Helper: Log-Sum-Exp for Numerical Stability
logsumexp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}





viterbi_om_pure <- function(OM_upper,
                            OM_lower,
                            nstates,
                            Gamma,
                            lambda_om = 1.0,
                            use_upper = TRUE,
                            use_lower = TRUE) {
  
  T_len <- nrow(OM_upper)
  log_Gamma <- log(Gamma)
  
  # ===== Step 1: 预处理 O-M 成本 (Emission costs) =====
  L_total <- matrix(0, T_len, nstates)
  for (t in 1:T_len) {
    for (i in 1:nstates) {
      val <- 0
      if (use_upper) val <- val + lambda_om * OM_upper[t, i]
      if (use_lower) val <- val + lambda_om * OM_lower[t, i]
      L_total[t, i] <- val
    }
  }
  
  # ===== Step 2: 初始化 Viterbi 矩阵 =====
  # v[i, t] 存储到时刻 t 状态为 i 的最大累积对数概率
  v <- matrix(-Inf, nstates, T_len)
  # psi[i, t] 存储回溯路径（哪个前置状态导致了当前的最大概率）
  psi <- matrix(0, nstates, T_len)
  
  # 初始分布 (Stationary)
  delta <- oeli::stationary_distribution(Gamma)
  
  for (i in 1:nstates) {
    v[i, 1] <- log(delta[i]) - L_total[1, i]
  }
  
  # ===== Step 3: 前向递归 (Recursion) =====
  for (t in 2:T_len) {
    for (j in 1:nstates) {
      # 寻找最大跳转概率: max_i [ v_{t-1}(i) + log(P_ij) ]
      trans_probs <- v[, t-1] + log_Gamma[, j]
      best_prev_state <- which.max(trans_probs)
      
      v[j, t] <- trans_probs[best_prev_state] - L_total[t, j]
      psi[j, t] <- best_prev_state
    }
  }
  
  # ===== Step 4: 回溯路径 (Backward Tracing) =====
  states_viterbi <- numeric(T_len)
  
  # 终点时刻的最优状态
  states_viterbi[T_len] <- which.max(v[, T_len])
  
  # 从 T-1 到 1 回溯
  for (t in (T_len - 1):1) {
    states_viterbi[t] <- psi[states_viterbi[t + 1], t + 1]
  }
  
  cat("Viterbi decoding complete.\n")
  return(list(
    states_estimate = states_viterbi,
    viterbi_score = v,
    path_matrix = psi
  ))
}











# ================================================================
# Plot Viterbi Results with O-M Functional (Probability Scale)
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
  
  # ===== Step 2: Run Forward Filter with O-M =====
  cat("Running Forward Filter with O-M...\n")
  
  observations <- batch_results_complete$v
  
  prob <- forward_filter_om(
    observations = observations,
    OM_upper = OM_upper,
    OM_lower = OM_lower,
    nstates = nstates,
    Gamma = Gamma,
    kappa = kappa,
    theta = theta,
    sigma = sigma,
    lambda_om = lambda_om,
    use_upper = use_upper,
    use_lower = use_lower
  )
  
  states_estimate <- prob$states_estimate
  posterior <- prob$posterior
  
  # ===== Step 3: Prepare data for plotting =====
  
  time_index <- 1:ncol(posterior)
  
  # Probability data
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = posterior[1, ],
    Prob_State2 = posterior[2, ]
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
    True_State = as.factor(Reg_chain[1:ncol(posterior)])
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
    
    # Posterior probability lines
    geom_line(
      data = prob_data_long,
      aes(x = Time, y = Probability, color = Probability_Type),
      linewidth = 0.9
    ) +
    
    # Estimated states
    geom_point(
      data = state_data,
      aes(x = Time, y = Estimated_State, shape = "Estimated State"),
      size = 1,
      color = "black",
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
        "Prob_State1" = expression(P(s[t] == 1 * " | " * bold(y)[1:t], bold(phi)[1:t])),
        "Prob_State2" = expression(P(s[t] == 2 * " | " * bold(y)[1:t], bold(phi)[1:t]))
      ),
      name = "Forward Probability (O-M Enhanced)"
    ) +
    scale_shape_manual(
      values = c("Estimated State" = 1),
      name = "Regime Estimate"
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      name = "Posterior Probability / State"
    ) +
    labs(
      title = sprintf("HMM with O-M Functional (norm=%s)", 
                      normalize_method),
      x = "Time Step"
    ) +
    theme_minimal() +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2),
      shape = guide_legend(order = 3)
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
          "P(path | regime)"
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
      posterior = posterior,
      om_probabilities = om_contrib_data  # 返回概率数据
    ))
    
  } else {
    return(list(
      plot = main_plot,
      states_estimate = states_estimate,
      posterior = posterior
    ))
  }
}



plot_viterbi_om_path <- function(batch_results_complete,
                                 nstates, 
                                 Gamma, 
                                 Reg_chain,
                                 lambda_om = 1.0,
                                 use_upper = TRUE,
                                 use_lower = TRUE,
                                 normalize_method = "log") {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # ===== Step 1: 准备 O-M 矩阵 (保持原有逻辑) =====
  cat("Preparing O-M matrices...\n")
  
  if (normalize_method == "log") {
    OM_upper <- log(batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] + 1)
    OM_lower <- log(batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] + 1)
  } else if (normalize_method == "zscore") {
    OM_all <- c(batch_results_complete$OM_upper_R1, batch_results_complete$OM_upper_R2,
                batch_results_complete$OM_lower_R1, batch_results_complete$OM_lower_R2)
    OM_finite <- OM_all[is.finite(OM_all)]
    mean_om <- mean(OM_finite); sd_om <- sd(OM_finite)
    OM_upper <- (batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] - mean_om) / sd_om
    OM_lower <- (batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] - mean_om) / sd_om
  } else {
    OM_upper <- batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")]
    OM_lower <- batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")]
  }
  
  OM_upper <- as.matrix(OM_upper)
  OM_lower <- as.matrix(OM_lower)
  
  # ===== Step 2: 运行 Viterbi 算法 (替换原有的 Forward Filter) =====
  cat("Running Viterbi Algorithm with O-M...\n")
  
  # 调用之前定义的 viterbi_om_pure 函数
  viterbi_results <- viterbi_om_pure(
    OM_upper = OM_upper,
    OM_lower = OM_lower,
    nstates = nstates,
    Gamma = Gamma,
    lambda_om = lambda_om,
    use_upper = use_upper,
    use_lower = use_lower
  )
  
  # Viterbi 得到的是 1, 2 序列，转为 0, 1 以匹配 Reg_chain
  viterbi_path <- viterbi_results$states_estimate - 1
  
  # ===== Step 3: 准备绘图数据 =====
  time_index <- 1:length(viterbi_path)
  
  state_data <- data.frame(
    Time = time_index,
    Viterbi_State = viterbi_path,
    True_State = as.factor(Reg_chain[1:length(viterbi_path)])
  )
  
  # 背景矩形数据 (保持原有逻辑)
  regime_background <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max(time_index) + 1)) %>%
    rename(Time_Start = Time)
  
  # ===== Step 4: 创建 Viterbi 路径图 =====
  
  p <- ggplot() +
    # 背景着色：代表真实的 Regime
    geom_rect(
      data = regime_background,
      aes(xmin = Time_Start, xmax = Time_End, ymin = -0.05, ymax = 1.05, fill = True_State),
      alpha = 0.15
    ) +
    
    # 真实的 Regime 阶梯线 (虚线)
    geom_step(
      data = state_data,
      aes(x = Time, y = as.numeric(as.character(True_State)), linetype = "True Path"),
      color = "gray40", linewidth = 0.8
    ) +
    
    # Viterbi 识别出的最可能路径 (粗实线)
    geom_step(
      data = state_data,
      aes(x = Time, y = Viterbi_State, color = "Viterbi Path"),
      linewidth = 1.2
    ) +
    
    # 格式设置
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "True Regime 1", "1" = "True Regime 2"),
      name = "Background (Ground Truth)"
    ) +
    scale_color_manual(
      values = c("Viterbi Path" = "darkblue"),
      name = "Inference"
    ) +
    scale_linetype_manual(
      values = c("True Path" = "dotted"),
      name = "Reference"
    ) +
    scale_y_continuous(
      breaks = c(0, 1),
      labels = c("Regime 1", "Regime 2"),
      limits = c(-0.1, 1.1)
    ) +
    labs(
      title = "Regime Detection: Viterbi Global Decoding with O-M Functional",
      subtitle = sprintf("Lambda_OM = %.2f | Normalization = %s", lambda_om, normalize_method),
      x = "Time Step",
      y = "Decoded State"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(
    plot = p,
    viterbi_path = viterbi_path,
    accuracy = mean(viterbi_path == as.numeric(as.character(state_data$True_State)))
  ))
}

# ================================================================
# Usage Example
# ================================================================

plot_result <- plot_viterbi_om(
  batch_results_complete = batch_results_complete,
  nstates = 2,
  Gamma = Gamma,
  kappa = c(10, 5),
  theta = c(0.1, 0.5),
  sigma = c(0.1, 0.1),
  Reg_chain = Reg_chain_year,
  lambda_om = 1,
  normalize_method = "log",
  show_om_contribution = TRUE
)











