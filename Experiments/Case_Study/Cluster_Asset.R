
source("Model_Simulation.R")
T =1 
N = 1000
v0 = 0.3
S0 =100


Reg_chain <- simulate_Reg(series_length = N)
plot(Reg_chain)
Reg_param <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta, sigma, rho
    0.5,   10,    0.3,   0.05,  -0.1, # calm
    0.5,   5,     0.6,   0.05,  -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)


# 62-144 regime change
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T = 0.5, N, M=1, method = "E")
S_simulated <- sim_series$S_paths
plot(S_simulated,  type = "l")

V_simulated <- sim_series$V_paths
plot(V_simulated, type = "l")




log_returns <- diff(log(S_simulated)) 
volatility_proxy <- abs(log_returns) * sqrt(252)


plot(V_simulated, type = 'l', ylim = c(0,1))
lines(volatility_proxy, col = "blue")
lines(lowess(volatility_proxy, f = 0.2), col = "red")
lowess_proxy <- lowess(volatility_proxy, f = 0.2)
# volatility_proxy <- lowess(volatility_proxy, f = 0.001)
# volatility_proxy <- log(volatility_proxy)
states <- 2 

cluster_results <- stats::kmeans(
  volatility_proxy[!is.na(lowess_proxy)], 
  centers = states, 
  iter.max = 100, 
  nstart = 100 
)
cluster_assignment <- cluster_results$cluster 
plot(cluster_assignment, col = "blue")
lines(Reg_chain+1)












obs_per_day <- 10
total_obs <- length(S_simulated)

# 对数收益率的数量 (比价格观测少 1)
total_returns <- total_obs - 1

# 确定完整的天数 (即收益率数量可以被 obs_per_day 整除的天数)
num_full_days <- floor(total_returns / obs_per_day)


# 裁剪收益率向量，使其长度恰好可以被 obs_per_day 整除
required_returns_length <- num_full_days * obs_per_day

# 计算对数收益率: r[t] = log(S[t]) - log(S[t-1])
log_returns <- diff(log(S_simulated))

# 裁剪对数收益率向量
returns_trimmed <- log_returns[1:required_returns_length]

# --- 3. 计算已实现方差 (Realized Variance) ---

# 将收益率重塑为矩阵，每列代表一天的收益率
# 矩阵有 obs_per_day 行，num_full_days 列。
returns_matrix <- matrix(
  data = returns_trimmed,
  nrow = obs_per_day, 
  ncol = num_full_days,
  byrow = FALSE # 确保每列是一个交易日内的 10 个收益率
)

# 计算每日已实现方差 (RV^2)
# RV^2 = sum(r_j^2) for r_j within the day
RV_squared_daily <- colSums(returns_matrix^2)

# --- 4. 计算已实现波动率 (Realized Volatility) ---

RV_daily <- sqrt(RV_squared_daily)

# --- 5. 结果展示 ---

cat("每日观测次数 (obs_per_day):", obs_per_day, "\n")
cat("总共计算了多少天的 RV:", length(RV_daily), "\n\n")
cat("已实现波动率 (RV_daily) 估计值:\n")
print(RV_daily)

plot(V_simulated, type = 'l', ylim = c(0,1))
lines(rep(sqrt(RV_daily), each = 10), col = "blue")



