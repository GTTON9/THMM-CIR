



library(ggplot2)

# ============ 各种波动率估计器 ============

# 1. Realized Volatility (RV) - 已有的
calculate_RV <- function(S, n_days, n_intraday) {
  RV <- numeric(n_days)
  for (t in 1:n_days) {
    idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
    S_day <- S[idx]
    r_day <- diff(log(S_day))
    RV[t] <- sum(r_day^2) * 252
  }
  return(RV)
}

# 2. Parkinson Volatility (High-Low)
calculate_Parkinson <- function(S, n_days, n_intraday) {
  PK <- numeric(n_days)
  for (t in 1:n_days) {
    idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
    S_day <- S[idx]
    
    H <- max(S_day)  # High
    L <- min(S_day)  # Low
    
    # Parkinson公式
    PK[t] <- (1 / (4 * log(2))) * (log(H / L))^2 * 252
  }
  return(PK)
}

# 3. Garman-Klass Volatility
calculate_GarmanKlass <- function(S, n_days, n_intraday) {
  GK <- numeric(n_days)
  for (t in 1:n_days) {
    idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
    S_day <- S[idx]
    
    O <- S_day[1]                    # Open
    H <- max(S_day)                  # High
    L <- min(S_day)                  # Low
    C <- S_day[length(S_day)]        # Close
    
    # Garman-Klass公式
    term1 <- 0.5 * (log(H / L))^2
    term2 <- -(2 * log(2) - 1) * (log(C / O))^2
    
    GK[t] <- (term1 + term2) * 252
  }
  return(GK)
}

# 4. Rogers-Satchell Volatility
calculate_RogersSatchell <- function(S, n_days, n_intraday) {
  RS <- numeric(n_days)
  for (t in 1:n_days) {
    idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
    S_day <- S[idx]
    
    O <- S_day[1]                    # Open
    H <- max(S_day)                  # High
    L <- min(S_day)                  # Low
    C <- S_day[length(S_day)]        # Close
    
    # Rogers-Satchell公式
    RS[t] <- sqrt(log(H / C) * log(H / O) + log(L / C) * log(L / O)) * sqrt(252)
  }
  return(RS^2)
}

# 5. Yang-Zhang Volatility (最复杂但最精确)
calculate_YangZhang <- function(S, n_days, n_intraday) {
  YZ <- numeric(n_days)
  for (t in 1:n_days) {
    idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
    S_day <- S[idx]
    
    O <- S_day[1]
    H <- max(S_day)
    L <- min(S_day)
    C <- S_day[length(S_day)]
    
    # Overnight volatility
    if (t > 1) {
      prev_idx <- ((t - 2) * n_intraday + 1):((t - 1) * n_intraday)
      C_prev <- S[prev_idx][length(S[prev_idx])]
      sigma_o <- (log(O / C_prev))^2
    } else {
      sigma_o <- 0
    }
    
    # Close-to-close volatility
    sigma_c <- (log(C / O))^2
    
    # Rogers-Satchell component
    sigma_rs <- log(H / C) * log(H / O) + log(L / C) * log(L / O)
    
    # Yang-Zhang combines these
    k <- 0.34 / (1.34 + (n_intraday + 1) / (n_intraday - 1))
    YZ[t] <- sqrt(sigma_o + k * sigma_c + (1 - k) * sigma_rs) * sqrt(252)
  }
  return(YZ^2)
}

# ============ 计算所有估计器 ============

S <- S_simulated
n_days <- 252
n_intraday <- 2400

cat("计算各种波动率估计器...\n")
plot(S[seq(1, length(S_simulated), by = n_intraday)], type = 'l', col = "black")
# 计算所有估计器
RV_raw <- calculate_RV(S, n_days, n_intraday)
PK_raw <- calculate_Parkinson(S, n_days, n_intraday)
GK_raw <- calculate_GarmanKlass(S, n_days, n_intraday)
RS_raw <- calculate_RogersSatchell(S, n_days, n_intraday)
YZ_raw <- calculate_YangZhang(S, n_days, n_intraday)

# 真实波动率
true_vol <- V_simulated[seq(1, length(V_simulated), by = n_intraday)]
smooth_f <- 0.1
RV_smooth <- lowess(RV_raw, f = smooth_f)$y
PK_smooth <- lowess(PK_raw, f = smooth_f)$y
GK_smooth <- lowess(GK_raw, f = smooth_f)$y
RS_smooth <- lowess(RS_raw, f = smooth_f)$y
YZ_smooth <- lowess(YZ_raw, f = smooth_f)$y


min_length <- min(length(true_vol), length(RV_smooth), length(PK_smooth), 
                  length(GK_smooth), length(RS_smooth), length(YZ_smooth))

cat("各向量长度:\n")
cat("true_vol:", length(true_vol), "\n")
cat("RV_smooth:", length(RV_smooth), "\n")
cat("PK_smooth:", length(PK_smooth), "\n")
cat("GK_smooth:", length(GK_smooth), "\n")
cat("RS_smooth:", length(RS_smooth), "\n")
cat("YZ_smooth:", length(YZ_smooth), "\n")
cat("使用长度:", min_length, "\n\n")

# 截取到相同长度
true_vol <- true_vol[1:min_length]
RV_smooth <- RV_smooth[1:min_length]
PK_smooth <- PK_smooth[1:min_length]
GK_smooth <- GK_smooth[1:min_length]
RS_smooth <- RS_smooth[1:min_length]
YZ_smooth <- YZ_smooth[1:min_length]
time_points <- 1:min_length
plot(V_simulated)
# 创建数据框
vol_comparison <- data.frame(
  time = rep(time_points, 6),
  volatility = c(true_vol, RV_smooth, PK_smooth, GK_smooth, RS_smooth, YZ_smooth),
  estimator = factor(rep(c("True", "RV", "Parkinson", "Garman-Klass", 
                           "Rogers-Satchell", "Yang-Zhang"), 
                         each = min_length),
                     levels = c("True", "RV", "Parkinson", "Garman-Klass", 
                                "Rogers-Satchell", "Yang-Zhang"))
)
# 绘制对比图
p1 <- ggplot(vol_comparison, aes(x = time, y = volatility, color = estimator)) +
  geom_line(size = 1, alpha = 0.8) +
  scale_color_manual(
    values = c("True" = "black", 
               "RV" = "blue", 
               "Parkinson" = "red",
               "Garman-Klass" = "green",
               "Rogers-Satchell" = "purple",
               "Yang-Zhang" = "orange")
  ) +
  labs(
    title = "Comparison of Volatility Estimators",
    subtitle = "All estimators smoothed with lowess (f=0.1)",
    x = "Time (Days)",
    y = "Volatility",
    color = "Estimator"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

print(p1)

# ============ 计算估计误差 ============

calculate_metrics <- function(estimated, true_val) {
  mae <- mean(abs(estimated - true_val))
  rmse <- sqrt(mean((estimated - true_val)^2))
  correlation <- cor(estimated, true_val)
  bias <- mean(estimated - true_val)
  
  return(c(MAE = mae, RMSE = rmse, Correlation = correlation, Bias = bias))
}

cat("\n=== 估计器性能指标 ===\n")
metrics_df <- data.frame(
  RV = calculate_metrics(RV_smooth, true_vol),
  Parkinson = calculate_metrics(PK_smooth, true_vol),
  GarmanKlass = calculate_metrics(GK_smooth, true_vol),
  RogersSatchell = calculate_metrics(RS_smooth, true_vol),
  YangZhang = calculate_metrics(YZ_smooth, true_vol)
)

print(round(metrics_df, 4))

# 保存指标
write.csv(metrics_df, "volatility_estimator_metrics.csv")

# ============ 分regime对比 ============

# 为每个估计器运行置信区间分析
estimators_list <- list(
  RV = RV_smooth,
  Parkinson = PK_smooth,
  GarmanKlass = GK_smooth,
  RogersSatchell = RS_smooth,
  YangZhang = YZ_smooth
)

all_results <- list()

for (est_name in names(estimators_list)) {
  cat(sprintf("\n计算 %s 的置信区间...\n", est_name))
  
  result <- plot_cir_confidence_simple(
    true_vol = estimators_list[[est_name]],
    param = param,
    states_estimate = states_estimate,
    dt = 1/252,
    interval_step = 10
  )
  
  all_results[[est_name]] <- result
  
  # 修改标题
  result$plot <- result$plot + 
    labs(title = sprintf("CIR Confidence Intervals - %s Estimator", est_name))
  
  print(result$plot)
  ggsave(sprintf("cir_ci_%s.png", tolower(est_name)), 
         result$plot, width = 12, height = 6, dpi = 300)
}

# ============ 组合对比图 ============

# 将所有置信区间数据合并
all_ci_combined <- data.frame()
for (est_name in names(all_results)) {
  ci_data <- all_results[[est_name]]$ci_data
  ci_data$estimator <- est_name
  all_ci_combined <- rbind(all_ci_combined, ci_data)
}

# 创建facet图
p2 <- ggplot() +
  # 真实波动率（黑色）
  geom_line(data = data.frame(time = time_points, true = true_vol),
            aes(x = time, y = true), 
            color = "black", size = 1, alpha = 0.8) +
  
  # 置信区间
  geom_point(data = all_ci_combined,
             aes(x = time, y = lower, fill = factor(regime)),
             shape = 21, size = 2, alpha = 0.6) +
  geom_point(data = all_ci_combined,
             aes(x = time, y = upper, fill = factor(regime)),
             shape = 21, size = 2, alpha = 0.6) +
  geom_segment(data = all_ci_combined,
               aes(x = time, xend = time, y = lower, yend = upper,
                   color = factor(regime)),
               alpha = 0.3, size = 0.5) +
  
  facet_wrap(~ estimator, ncol = 2, scales = "free_y") +
  
  scale_color_manual(values = c("1" = "blue", "2" = "darkgreen")) +
  scale_fill_manual(values = c("1" = "blue", "2" = "darkgreen")) +
  
  labs(
    title = "CIR Confidence Intervals Comparison Across Estimators",
    subtitle = "Black line: True volatility | Circles: 90% CI bounds",
    x = "Time (Days)",
    y = "Volatility",
    fill = "Regime"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p2)


cat("\n=== 置信区间覆盖率分析 ===\n")

coverage_analysis <- function(ci_data, true_vol) {
  # 检查真实值是否在置信区间内
  in_ci <- sapply(1:nrow(ci_data), function(i) {
    t_idx <- ci_data$time[i]
    true_val <- true_vol[t_idx]
    true_val >= ci_data$lower[i] && true_val <= ci_data$upper[i]
  })
  
  coverage <- mean(in_ci) * 100
  return(coverage)
}

coverage_results <- data.frame(
  Estimator = names(all_results),
  Coverage = sapply(names(all_results), function(est_name) {
    coverage_analysis(all_results[[est_name]]$ci_data, true_vol)
  })
)

print(coverage_results)




