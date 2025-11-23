# 假设 simulate_heston 已返回 paths
Reg_series <- simulate_Reg(series_length = N)
res <- simulate_heston(S0 = S0, v0 = v0, Reg_series = Reg_series, 
                       Reg_param = cbind(mu, kappa, theta, sigma, rho),
                       T = 1, N = 250, M = 1, method = "E_C", seed = 999, min_var = 1e-12)

V <- as.numeric(res$V_paths[1, ])   # 取第一条路径（单样本情况）
N <- length(V) - 1
dt <- 1 / 250  # 你之前 T/N, adjust if different

# compute regime-specific q, C for each time step
qs <- numeric(N)
Cvals <- numeric(N)
u_vals <- numeric(N)
v_vals <- numeric(N)
z_vals <- numeric(N)
z_over_q <- numeric(N)

for (i in 1:N) {
  reg_index <- Reg_series[i] + 1
  kappa_curr <- kappa[reg_index]
  theta_curr <- theta[reg_index]
  sigma_curr <- sigma[reg_index]
  exp_k_dt <- exp(-kappa_curr * dt)
  C_val <- (2 * kappa_curr) / ((1 - exp_k_dt) * sigma_curr^2)
  q_val <- (2 * kappa_curr * theta_curr) / (sigma_curr^2) - 1
  
  u <- C_val * V[i] * exp_k_dt
  v <- C_val * V[i+1]
  z <- 2 * sqrt(pmax(0, u * v))
  
  qs[i] <- q_val
  Cvals[i] <- C_val
  u_vals[i] <- u
  v_vals[i] <- v
  z_vals[i] <- z
  z_over_q[i] <- ifelse(abs(q_val) > 0, z / abs(q_val), NA)
}

summary_df <- data.frame(q = qs, C = Cvals, u = u_vals, v = v_vals, z = z_vals, z_over_q = z_over_q)
print(summary(summary_df))

# 查看极端值
cat("Percentiles of z/q:\n")
print(quantile(z_over_q, probs = c(0, 0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.99, 0.999, 1), na.rm=TRUE))

# 查看 if any C, u, v blow up
cat("C range, u range, v range:\n")
print(range(Cvals, na.rm=TRUE))
print(range(u_vals, na.rm=TRUE))
print(range(v_vals, na.rm=TRUE))
