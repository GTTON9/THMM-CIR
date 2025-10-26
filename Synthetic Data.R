set.seed(999) 

simulate_Reg <- function(
    series_length = 252,
    initial_state = 0,
    Reg_tran = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2),
    plot = TRUE
) {
  
  
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
         main = "Regime ",
         xlab = "Time Step",
         ylab = "Regime",
         yaxt = 'n' 
    )
    # 自定义Y轴标签为 0 和 1
    axis(2, at = c(0, 1), labels = c("Calm", "Turbulent"))
  }
  
  
  invisible(Reg_chain) 
}


# Reg_chain <- simulate_Reg(Reg_tran = matrix(c(1, 0, 0, 1), 2, 2))
Reg_chain <- simulate_Reg()

# Euler-Method



simulate_heston <- function(S0, v0, mu, Reg_series, Reg_param, rho, T, N, M, method = "E") {
  
  # M: number of path simulated by MC

  if (2 * Reg_param[1,1] * Reg_param[1,2] < Reg_param[1,3]^2 || 
      2 * Reg_param[2,1] * Reg_param[2,2] < Reg_param[2,3]^2) {
    message("Warning: Feller condition (2*kappa*theta > sigma^2) NOT satisfied.")
  }

  
  dt <- T / N           
  sqrt_dt <- sqrt(dt)
  
  
  if (rho < -1 || rho > 1) {
    stop("相关系数 rho 必须在 [-1, 1] 之间。")
  }
  
 
  method <- toupper(method) 
  if (!(method %in% c("E", "M"))) {
    stop("method 参数必须是 'E' (Euler) 或 'M' (Milstein)。")
  }
  

  rho_prime <- sqrt(1 - rho^2)
  
  # matrix for path (Monte Carlo)
  S_paths <- matrix(0, nrow = M, ncol = N + 1)
  V_paths <- matrix(0, nrow = M, ncol = N + 1)
  
  S_paths[, 1] <- S0
  V_paths[, 1] <- v0
  
  for (i in 1:N) {
    
    # stochastic term
    Z1 <- rnorm(M)  
    Z2 <- rnorm(M)  
    
    V_prev <- V_paths[, i]
    S_prev <- S_paths[, i]
    
    dW1 <- Z1 * sqrt_dt
    dW2 <- (rho * Z1 + rho_prime * Z2) * sqrt_dt
    
    sqrt_V <- sqrt(V_prev)
    
    print(Reg_series[i])
    if(Reg_series[i] == 0){
      kappa <- Reg_param[1,1]
      theta <- Reg_param[1,2]
      sigma <- Reg_param[1,3]
      
    }else{
      
      kappa_2 <- Reg_param[2,1]
      theta_2 <- Reg_param[2,2]
      sigma_2 <- Reg_param[2,3]
    }
    
    if (method == "E") { # Euler-Maruyama Scheme
      
      # V_t 更新：只包含欧拉项
      V_next <- (
        V_prev 
        + kappa * (theta - V_prev) * dt     # 漂移项
        + xi * sqrt_V * dW2            # 扩散项
      )
      
    } else if (method == "M") {
      # === 米尔斯坦修正法 (Milstein Scheme) for V_t ===
      
      # 欧拉项: f_v * dt + g_v * dW2
      V_euler_term <- (
        V_prev 
        + kappa * (theta - V_prev) * dt
        + xi * sqrt_V * dW2
      )
      
      # 米尔斯坦修正项: 1/4 * xi^2 * ((dW2)^2 - dt)
      V_milstein_term <- 0.25 * xi^2 * (dW2^2 - dt)
      
      V_next <- V_euler_term + V_milstein_term
    }
    
    # --- 方差 V_t 路径存储 (两种方法都需要最终截断) ---
    V_paths[, i + 1] <- pmax(V_next, 0)
    
    # --- 资产价格 S_t 的更新 (两种方法都使用欧拉项，以简化) ---
    # S_t 仍使用欧拉公式：
    S_paths[, i + 1] <- (
      S_prev 
      + mu * S_prev * dt                 # 漂移项
      + S_prev * sqrt_V * dW1       # 扩散项
    )
  }
  
  # 返回结果（以列表形式）
  return(list(S_paths = S_paths, V_paths = V_paths, method_used = method))
}

results_E <- simulate_heston(S0, v0, mu,Reg_chain, Reg_param, rho, T, N, M, method = "E")


kappa_1 <- 2.0; 
theta_1 <- 0.1; 
sigma_1 <- 0.5;

# turbulent regime
kappa_2 <- 2.0; 
theta_2 <- 0.2; 
sigma_2 <- 0.8;

Reg_param <- matrix(c(kappa_1, theta_1, sigma_1,
                      kappa_2, theta_2, sigma_2), nrow = 2)


S0 <- 100.0; v0 <- 0.04; mu <- 0.05;  rho <- -0.7; T <- 1.0; 
N <- 252; M <- 10000;

# 1. 使用欧拉法
results_E <- simulate_heston(S0, v0, mu,Reg_chain, Reg_param, rho, T, N, M, method = "E")
print("=== 欧拉法模拟结果 (E) ===")
print(paste0("最终价格 E[S_T]: ", round(mean(results_E$S_paths[, N + 1]), 4)))
print(paste0("最终方差 E[V_T]: ", round(mean(results_E$V_paths[, N + 1]), 4)))
matplot(t(results_E$S_paths[1:10, ]), type = "l", 
        main = "Heston Model Asset Price Paths", 
        xlab = "Time Step", ylab = "Asset Price S_t")

  
# 2. 使用米尔斯坦修正法 
results_M <- simulate_heston(S0, v0, mu,Reg_chain, Reg_param, rho, T, N, M, method = "M")
print("=== 米尔斯坦修正法模拟结果 (M) ===")
print(paste0("最终价格 E[S_T]: ", round(mean(results_M$S_paths[, N + 1]), 4)))
print(paste0("最终方差 E[V_T]: ", round(mean(results_M$V_paths[, N + 1]), 4)))


matplot(t(results_M$S_paths[1:10, ]), type = "l", 
        main = "Heston Model Asset Price Paths", 
        xlab = "Time Step", ylab = "Asset Price S_t")







# --- 1. Heston 模拟数据的准备 ---

# 假设 N, S_paths, results_E 已经存在于环境中
N <- 252 
S_simulated <- results_E$S_paths[1, ] 

# 计算对数收益率
simulated_returns <- diff(log(S_simulated)) 
series_length <- length(simulated_returns)

# 创建假设的日期序列
# 从 2024-01-01 开始创建 series_length 个连续日期
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)

# 将收益率和日期放入 data.frame
my_data_df <- data.frame(
  Date = date_sequence,
  returns = simulated_returns
)

# --- 2. 设置 fHMM_controls 参数 ---

# 使用 set_controls 正确地创建 fHMM_controls 对象
controls_simulated <- set_controls( 
  states      = 2,              # 拟合 2 个隐藏状态
  sdds        = "t",            # 使用 t 分布 (Student's t-distribution)
  file        = my_data_df,     # 传入数据框
  data_column = "returns",      # 指定收益率数据列名
  logreturns  = FALSE,          # 数据已经是收益率
  
  # 使用索引而非日期列
  from        = 1,              # 序列起始索引
  to          = series_length   # 序列结束索引 (即 N)
)

# --- 3. 准备数据并拟合模型 ---

# 准备 fHMM 需要的数据结构
data_hmm <- prepare_data(controls_simulated)

# 检查数据摘要
print("--- 准备数据摘要 ---")
summary(data_hmm) 

# 拟合 HMM 模型
model_hmm <- fit_model(data_hmm)

# --- 4. 结果分析与可视化 ---

# 解码状态链
model_hmm <- decode_states(model_hmm)

# 打印模型拟合结果（包括转移矩阵和状态参数）
print("--- HMM 拟合模型摘要 ---")
summary(model_hmm)

# 绘制时间序列图和状态图 (假设您要看第一条路径的状态切换)
plot(model_hmm, plot_type = c("ts", "sdds"))

library(fHMM)
dax <- download_data(symbol = "^GDAXI")
syn_data <- results_M
controls <- set_controls( 
  states      = 3,
  sdds        = "t",
  file        = dax,
  date_column = "Date",
  data_column = "Close",
  logreturns  = TRUE,
  from        = "2000-01-01",
  to          = "2022-12-31"
)
?set_controls
data <- prepare_data(controls)
summary(data)
model <- fit_model(data)
model <- decode_states(model)
model <- compute_residuals(model)
summary(model)

events <- fHMM_events(
  list(dates = c("2001-09-11", "2008-09-15", "2020-01-27"),
       labels = c("9/11 terrorist attack", "Bankruptcy Lehman Brothers", "First COVID-19 case Germany"))
)
plot(model, plot_type = c("sdds","ts"), events = events)
plot(model, plot_type = "pr")

# Calm Regime/Turbuletn Regime
# reversion rate similar, 
# theta_1 should be lower than theta_2
# sigma_1 < sigma_2
fixed paramters, 
# 欧拉法
# start with OU-process, no discretization problems
