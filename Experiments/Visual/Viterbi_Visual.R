
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


plot_viterbi <- function(V_simulated, nstates, Gamma, kappa, theta, sigma, Reg_chain){
  
  
  prob <- forward_filter(V_simulated, nstates, Gamma, kappa, theta, sigma)
  states_estimate <- as.vector(prob$states_estimate)
  a <- prob$posterior
  # states_estimate <- viterbi(V_simulated, nstates, Gamma, kappa, theta, sigma)
  
  time_index <- 1:ncol(a) # 假设时间索引从 1 到 N
  
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = a[1,],
    Prob_State2 = 1 - a[1,]
  )
  
  
  state_data <- data.frame(
    Time = time_index,
    Estimated_State = states_estimate - 1, 
    True_State = Reg_chain - 1        
  )

  
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  
  
  
  
  time_index <- 1:ncol(a) # 假设时间索引从 1 到 N
  
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = a[1,],
    Prob_State2 = 1 - a[1,]
  )
  
  # --- 步骤 2: 创建状态链数据框 ---
  state_data <- data.frame(
    Time = time_index,
    Estimated_State = states_estimate - 1, # 确保状态为 0/1
    True_State = Reg_chain - 1             # 确保状态为 0/1
  )
  
  # --- 步骤 3: 转换概率数据为长格式 (用于 ggplot 的多线绘制) ---
  library(tidyr)
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  
  
plot <-   ggplot() +
    
    # I. 背景色区域 (突出 True Regime)
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
    
    # II. 后验概率线 (主要信息)
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
    
    # --- 调整和标签 ---
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"), 
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)"),
      name = "True Regime Background"
    ) +
    scale_color_manual(
      values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
      labels = c("Prob_State1" = "P(State 1)", "Prob_State2" = "P(State 2)"),
      name = "Posterior Probability"
    ) +
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
