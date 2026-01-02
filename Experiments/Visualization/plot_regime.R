# 确保数据是数值型
Reg_chain_year <- as.numeric(Reg_chain_year)
n <- length(Reg_chain_year)
time_index <- 1:n

# 1. 初始化绘图窗口
# 为了美观，我们将 0 映射到 y=1，将 1 映射到 y=2
plot(time_index, Reg_chain_year + 1, type = "n", 
     yaxt = "n", 
     ylab = "True Regime", 
     xlab = "Time (1 year)",
     main = "True Regimes Over Time",
     ylim = c(0.5, 2.5))

# 2. 自定义 Y 轴刻度
# 逻辑：0+1=1 (Regime 1), 1+1=2 (Regime 2)
axis(2, at = c(1, 2), labels = c("Reg 1", "Reg 2"), las = 1)

# 3. 循环画线
for(i in 1:(n-1)) {
  curr <- Reg_chain_year[i]
  next_v <- Reg_chain_year[i+1]
  
  # 颜色逻辑：0 (Regime 1) 用红色(火砖色)，1 (Regime 2) 用蓝色
  line_color <- ifelse(curr == 0, "firebrick", "royalblue")
  
  # 画水平线（注意高度加 1）
  segments(i, curr + 1, i + 1, curr + 1, col = line_color, lwd = 3)
  
  # 状态切换时的垂直虚线
  if(curr != next_v) {
    segments(i + 1, curr + 1, i + 1, next_v + 1, col = "gray80", lty = 3)
  }
}

# 4. 图例
legend("topright", 
       legend = c("Regime 1 (Turbulent)", "Regime 2 (Calm)"), 
       col = c("firebrick", "royalblue"), 
       lwd = 3, bg = "white", cex = 0.8)

grid(nx = NULL, ny = NA, col = "gray90")