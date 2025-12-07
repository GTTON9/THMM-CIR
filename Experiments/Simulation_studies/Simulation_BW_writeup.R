source("Gen_fit.R")



# easy case
N <- 252
v0 <- 0.03
S0 <- 100
nstates = 2
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(10, 5)
theta <- c(0.03, 0.6)
sigma <- c(0.2, 0.2)
rho <- c(-0.1, -0.2)

BW_result <- Gen_fit(Gamma, mu, kappa, theta, sigma, rho, V_ = TRUE , plot_path = TRUE)
BW_result$param
BW_result$states_estimate
BW_result$fisher_inverse


variance <- BW_result$fisher_inverse
theta_hat <- c(
  BW_result$param$Gamma[1, 1], # Gamma_11
  BW_result$param$Gamma[2, 2], # Gamma_22
  BW_result$param$kappa[1],    # kappa_1
  BW_result$param$kappa[2],    # kappa_2
  BW_result$param$theta[1],    # theta_1
  BW_result$param$theta[2],    # theta_2
  BW_result$param$sigma[1],    # sigma_1
  BW_result$param$sigma[2]     # sigma_2
)
Z_score <- qnorm(0.975) 

# 计算标准误差 (Standard Errors)
standard_errors <- sqrt(variance)

# 计算误差边际 (Margin of Error)
margin_of_error <- Z_score * standard_errors

# 计算 95% 置信区间的上下限
lower_bound <- theta_hat - margin_of_error
upper_bound <- theta_hat + margin_of_error

# -----------------
# 4. 整理结果并输出
# -----------------

# 参数名称向量 (与 theta_hat 顺序一致)
param_names <- c(
  "Gamma_11", "Gamma_22", 
  "kappa_1", "kappa_2", 
  "theta_1", "theta_2", 
  "sigma_1", "sigma_2"
)

# 组合成数据框
confidence_intervals <- data.frame(
  Parameter = param_names,
  Estimate = theta_hat,
  Std_Error = standard_errors,
  Z_score = Z_score,
  CI_95_Lower = lower_bound,
  CI_95_Upper = upper_bound
)

# 打印结果
print("参数的 95% Wald 置信区间:")
print(confidence_intervals)











# hard case


N <- 252
v0 <- 1
S0 <- 100
nstates = 2
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(2, 1)
theta <- c(0.1, 0.2)
sigma <- c(0.05, 0.05)
rho <- c(-0.1, -0.1)


source("Gen_fit.R")
BW_result <- Gen_fit(Gamma, mu, kappa, theta, sigma, rho, V_ = TRUE , plot_path = TRUE)
BW_result$param
BW_result$states_estimate



variance <- BW_result$fisher_inverse
theta_hat <- c(
  BW_result$param$Gamma[1, 1], # Gamma_11
  BW_result$param$Gamma[2, 2], # Gamma_22
  BW_result$param$kappa[1],    # kappa_1
  BW_result$param$kappa[2],    # kappa_2
  BW_result$param$theta[1],    # theta_1
  BW_result$param$theta[2],    # theta_2
  BW_result$param$sigma[1],    # sigma_1
  BW_result$param$sigma[2]     # sigma_2
)
Z_score <- qnorm(0.975) 

# 计算标准误差 (Standard Errors)
standard_errors <- sqrt(variance)

# 计算误差边际 (Margin of Error)
margin_of_error <- Z_score * standard_errors

# 计算 95% 置信区间的上下限
lower_bound <- theta_hat - margin_of_error
upper_bound <- theta_hat + margin_of_error

# -----------------
# 4. 整理结果并输出
# -----------------

# 参数名称向量 (与 theta_hat 顺序一致)
param_names <- c(
  "Gamma_11", "Gamma_22", 
  "kappa_1", "kappa_2", 
  "theta_1", "theta_2", 
  "sigma_1", "sigma_2"
)

# 组合成数据框
confidence_intervals <- data.frame(
  Parameter = param_names,
  Estimate = theta_hat,
  Std_Error = standard_errors,
  Z_score = Z_score,
  CI_95_Lower = lower_bound,
  CI_95_Upper = upper_bound
)

# 打印结果
print("参数的 95% Wald 置信区间:")
print(confidence_intervals)





