source("Model_Simulation.R")
source("Gen_fit.R")

v0 <- 0.03
S0 <- 100
nstates = 2
T = 10
n_days <- 2520
n_intraday <- 2400
Gamma <- matrix(c(0.999, 0.001, 0.001, 0.999), 2, 2)
mu <- c(0.1, 0.3)
kappa <- c(10, 5)
theta <- c(0.03, 0.6)
sigma <- c(0.1, 0.1)
rho <- c(-0.1, -0.1)


N <- 2520 * n_intraday
Reg_chain_year <- simulate_Reg(series_length = n_days, Reg_tran = Gamma)
Reg_chain <- rep(Reg_chain_year, each = n_intraday)
Reg_param <- cbind(mu, kappa, theta, sigma, rho)
sim_series_10_year  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E",interp = F, seed = 999)
S_simulated <- sim_series_10_year$S_paths
# plot(S_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)], type = "l")
V_simulated <- sim_series_10_year$V_paths


par(mfrow = c(1, 2))

plot(
  S_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)],
  type = "l",
  main = "S path",
  col = "blue")

plot(
  V_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)],
  type = "l",
  main = "V path",
  col = "red"
)

par(mfrow = c(1, 1))




S <- S_simulated


RV_V <- numeric(n_days)

for (t in 1:n_days) {
  idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
  S_day <- S[idx]
  r_day <- diff(log(S_day))      # length 3
  RV_V[t] <- sum(r_day^2) * 252 
}


plot(RV_V, col = "grey")
lines(V_simulated[seq(1, length(V_simulated), by = n_intraday)], col = "blue")
# lines(lowess(RV_V, f = 0.01), col = "red")



# RV_V <- lowess(RV_V, f = 0.03)$y




series_length <- length(RV_V)
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)







source("Gen_fit.R")
my_data_df <- data.frame(
  Date = date_sequence,
  Var = RV_V
)
colnames(my_data_df) <- c("Date", "Var")



# followed the example
series_control <- Heston_set_controls( 
  states      = 2,     # 2 state
  sdds        = "Heston",         
  date_column = "Date",
  file        = my_data_df, 
  data_column = "Var",      
  logreturns  = FALSE,         
  from        = date_sequence[1],             
  to          = date_sequence[length(date_sequence)],
  runs = 10
)


data_hmm <- prepare_data(series_control)
model_hmm <- Heston_fit_model(data_hmm) 

final_model <- decode_states_heston(model_hmm) 
states_estimate <- final_model$decoding
states_estimate

plot(states_estimate, col = "blue")
lines(Reg_chain_year+1)
param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)


p <- plot_viterbi(RV_V, nstates, param$Gamma, 
                  param$kappa, param$theta, 
                  param$sigma, Reg_chain_year)
plot(p)
param







calculate_cir_densities_by_regime <- function(time_points, estimated_vol, true_vol, 
                                              param, states_estimate, 
                                              dt = 1/252, 
                                              V_grid = seq(0.001, 1, length.out = 100)) {
  
  # 确保所有输入长度一致
  min_length <- min(length(time_points), length(estimated_vol), 
                    length(true_vol), length(states_estimate))
  
  time_points <- time_points[1:min_length]
  estimated_vol <- estimated_vol[1:min_length]
  true_vol <- true_vol[1:min_length]
  states_estimate <- states_estimate[1:min_length]
  
  cat("使用数据长度:", min_length, "\n")
  
  # 找到regime转换点
  changepoints <- c(1, which(diff(states_estimate) != 0) + 1, length(states_estimate) + 1)
  
  # 存储所有密度数据
  all_density_data <- list()
  
  # 为每个regime段计算密度
  for (i in 1:(length(changepoints) - 1)) {
    start_idx <- changepoints[i]
    end_idx <- changepoints[i + 1] - 1
    
    current_regime <- states_estimate[start_idx]
    kappa_reg <- param$kappa[current_regime]
    theta_reg <- param$theta[current_regime]
    sigma_reg <- param$sigma[current_regime]
    
    regime_times <- time_points[start_idx:end_idx]
    regime_estimated <- estimated_vol[start_idx:end_idx]
    
    # 使用regime的第一个点作为参考
    v_ref <- regime_estimated[1]
    
    # 为这个regime计算不同时间点的密度
    for (j in 1:length(regime_estimated)) {
      t_idx <- regime_times[j]
      k <- j * dt  # 从regime起点开始的时间间隔
      
      # 计算整个V_grid上的密度
      log_densities <- sapply(V_grid, function(v) {
        ln_d_Heston(V_t = v_ref, 
                    V_t_plus_k = v, 
                    k = k, 
                    kappa = kappa_reg, 
                    theta = theta_reg, 
                    sigma = sigma_reg)
      })
      
      densities <- exp(log_densities)
      
      # 存储结果
      density_result <- data.frame(
        time = t_idx,
        V_grid = V_grid,
        density = densities,
        log_density = log_densities,
        regime = current_regime,
        true_vol = true_vol[t_idx],
        estimated_vol = estimated_vol[t_idx]
      )
      
      all_density_data[[length(all_density_data) + 1]] <- density_result
    }
  }
  
  # 合并所有结果
  density_df <- do.call(rbind, all_density_data)
  
  cat("\n密度计算完成，共", nrow(density_df), "个数据点\n")
  
  return(list(
    density_data = density_df,
    changepoints = changepoints,
    time_points = time_points
  ))
}












library(ggplot2)




true_vol <- V_simulated[seq(1, length(V_simulated), by = n_intraday)]


result <- plot_cir_confidence_simple(
  true_vol = true_vol,
  param = param,
  states_estimate = states_estimate,
  dt = 1/252,
  interval_step = 1
)


print(result$plot)



print(summary(result$ci_data))


param


S_daily <- S_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)]



result_mle <- estimate_mu_simple_mle(
  S_daily = S_daily,
  S_intraday = S_simulated,
  n_intraday = n_intraday,
  states_estimate = states_estimate,
  nstates = 2,
  dt = 1/252
)

result_mle$mu



