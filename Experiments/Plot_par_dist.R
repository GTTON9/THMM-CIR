# 1. 加载必要的库
library(tidyverse)

# --- 2. 定义输入数据 ---
# ... (此部分保持不变) ...

# Fisher 逆矩阵的对角线元素 (方差)
FIM_Variances <- BW_result$fisher_inverse
names(FIM_Variances) <- c("Gamma_11", "Gamma_22", "kappa_1", "kappa_2", "theta_1", "theta_2", "sigma_1", "sigma_2")

# 点估计值 
Estimates <- c(
  BW_result$param$Gamma[1, 1], # Gamma_11
  BW_result$param$Gamma[2, 2], # Gamma_22
  BW_result$param$kappa[1],    # kappa_1
  BW_result$param$kappa[2],    # kappa_2
  BW_result$param$theta[1],    # theta_1
  BW_result$param$theta[2],    # theta_2
  BW_result$param$sigma[1],    # sigma_1
  BW_result$param$sigma[2]     # sigma_2
)
names(Estimates) <- names(FIM_Variances)

# 真实参数值 (以您更新后的值为准)
True_Values <- c(0.99, 0.99, 10, 5, 0.03, 0.6, 0.05, 0.05)

# True_Values <- c(0.99, 0.99, 1, 2, 0.1, 0.2, 0.05, 0.05)

names(True_Values) <- names(FIM_Variances)

# --- 3. 计算渐近 95% CI 并创建绘图数据框 ---

# 计算标准误差 (SE = sqrt(方差))
SE <- sqrt(FIM_Variances)

# 计算 95% 置信区间 (使用 Z 分数 1.96)
Z_score <- 1.96
CI_95_Lower_FIM <- Estimates - Z_score * SE
CI_95_Upper_FIM <- Estimates + Z_score * SE

# 修正：创建最终绘图数据框，**包含 SE 列**
df_plot_fim <- data.frame(
  Parameter = names(Estimates),
  Estimate = Estimates,
  SE = SE, 
  CI_lower = CI_95_Lower_FIM,
  CI_upper = CI_95_Upper_FIM,
  True_Value = True_Values
)

param_order <- names(FIM_Variances)
df_plot_fim$Parameter <- factor(df_plot_fim$Parameter, levels = param_order)

# --- 4. 绘制 2x4 渐近正态分布图 ---

# 1. 生成正态分布曲线数据：(函数体保持不变，确保了稳定性)
generate_normal_data <- function(df) {
  df$Parameter <- as.character(df$Parameter) 
  
  data_list <- lapply(1:nrow(df), function(i) {
    param_data <- df[i, ]
    estimate_val <- param_data$Estimate[[1]]
    se_val <- param_data$SE[[1]]
    param_name <- param_data$Parameter[[1]]
    
    x_range <- seq(estimate_val - 4 * se_val, 
                   estimate_val + 4 * se_val, 
                   length.out = 200)
    
    data.frame(
      Parameter = param_name,
      x = x_range,
      y = dnorm(x_range, mean = estimate_val, sd = se_val)
    )
  })
  return(bind_rows(data_list))
}

df_distribution <- generate_normal_data(df_plot_fim)

# 确保分布数据中的参数顺序正确
df_distribution$Parameter <- factor(df_distribution$Parameter, levels = param_order)


# **新增步骤：计算每个参数的最大密度值，用于定位标签的 Y 轴高度**
df_y_max <- df_distribution %>%
  group_by(Parameter) %>%
  summarise(y_max = max(y, na.rm = TRUE))

df_plot_labels <- df_plot_fim %>%
  left_join(df_y_max, by = "Parameter")


# 2. 绘制图表
plot_fim_distribution <- ggplot(df_distribution, aes(x = x, y = y)) +
  
  # 绘制正态密度曲线
  geom_line(color = "darkblue", linewidth = 1) +
  geom_area(fill = "lightblue", alpha = 0.6) +
  
  # 添加 True Value (红色实线)
  geom_vline(data = df_plot_fim, 
             aes(xintercept = True_Value), 
             color = "red", 
             linetype = "solid", 
             linewidth = 1) +
  
  # 添加 95% CI 边界 (黑色虚线)
  geom_vline(data = df_plot_fim, 
             aes(xintercept = CI_lower), 
             color = "black", 
             linetype = "dashed", 
             linewidth = 0.5) +
  geom_vline(data = df_plot_fim, 
             aes(xintercept = CI_upper), 
             color = "black", 
             linetype = "dashed", 
             linewidth = 0.5) +
  
  # **新增部分：添加真实值标签**
  geom_text(data = df_plot_labels,
            aes(x = True_Value,
                y = 0.95 * y_max, # 将标签放置在接近最大密度处
                label = sprintf("%.2f", True_Value)), # 格式化为 4 位小数
            hjust = -0.1, # 将标签稍微向右推，避免重叠
            color = "red",
            size = 3.5,
            fontface = "bold") +
  
  # 分面为 2 行 4 列，scales="free" 确保每个子图的轴适应其参数范围
  facet_wrap(~ Parameter, ncol = 4, scales = "free") +
  
  # 样式和标签设置
  labs(
    title = "8 CIR paramter Estimation (Asymptotic Normal Distribution)",
    x = "Parameters",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

# 打印绘图结果
print(plot_fim_distribution)





