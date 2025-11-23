if (q > 50 && x < 100) { # based on the application， stirling approximation
  
  scaled_I_q <- (1 / sqrt(2*pi*q)) * (exp(1) * x / (2*q))^q
  ln_I_q <- log(scaled_I_q)
  # besselI(100, 50)
  
} else if(q > 50 && x > 100){
  # need reference, depends on ratio (Airy functions)
  print("numeric warning")
  
} else { 
  scaled_I_q <- besselI(x, nu = q, expon.scaled = TRUE)
  ln_I_q <- 1*x + log(scaled_I_q)
}



## x < 100 && q > 50
x = 80
q =50
scaled_I_q <- (1 / sqrt(2*pi*q)) * (exp(1) * x / (2*q))^q
ln_I_q <- log(scaled_I_q)
ln_I_q

scaled_I_q <- besselI(x, nu = q, expon.scaled = TRUE)
ln_I_q <- 1*x + log(scaled_I_q)


## 0 < x < 100 && 0 < q < 50 
x = 5000
q = 25000

scaled_I_q <- besselI(x, nu = q, expon.scaled = TRUE)
ln_I_q <- 1*x + log(scaled_I_q)
ln_I_q
 

# z >1
Debye_Bessel <- function(x, q){
  z <- x/q
  eta <- sqrt(1 + z^2) + log(z / (1 + sqrt(1 + z^2)))
  nom <- exp(q * eta )
  denom <- sqrt(2 * pi * q) * (1 + z^2)^(1/4)
  return(nom/denom)
}
Debye_Bessel(x,q)


z <- x/q
eta <- sqrt(1 + z^2) + log(z / (1 + sqrt(1 + z^2)))
nom <- exp(q * eta )
denom <- sqrt(2 * pi * q) * (1 + z^2)^(1/4)
nom/denom




logI_debye <- function(x, nu) {
  z <- x / nu
  sqrt1z2 <- sqrt(1 + z^2)
  eta <- sqrt1z2 + log(z / (1 + sqrt1z2))
  log_pref <- nu * eta - 0.5 * log(2 * pi * nu) - 0.25 * log(1 + z^2)
  return(log_pref)
}

logI_debye()






# 函数: 计算 I_nu(x) 的对数渐近近似值
# x: 自变量 (Argument)
# nu: 阶数 (Order, 必须很大)
log_I_asymptotic_small_z <- function(x, nu) {
  
  if (nu <= 0 || x <= 0) {
    stop("Parameters 'x' and 'nu' must be positive.")
  }
  
  # lgamma(nu + 1) 计算 ln(Gamma(nu + 1))
  # nu * log(x/2) 计算 ln((x/2)^nu)
  log_approx <- nu * log(x / 2) - lgamma(nu + 1)
  
  return(log_approx)
}

# --- 示例用法 ---
nu_val <- 25000  # 大阶数
x_val <- 5000    # 相对较小的自变量 (z = 5/20 = 0.25 << 1)

# 计算 log(I_nu(x)) 的近似值
log_I_approx_result <- log_I_asymptotic_small_z(x_val, nu_val)
cat(sprintf("log(I_%.1f(%.1f)) 的渐近近似值: %.4f\n", 
            nu_val, x_val, log_I_approx_result))

# 计算 I_nu(x) 的原始近似值
I_approx_result <- exp(log_I_approx_result)
cat(sprintf("I_%.1f(%.1f) 的近似值: %.4e\n", 
            nu_val, x_val, I_approx_result))

# 对比 R 内置精确值
# 注意：I_true 的值非常小，使用 log(besselI(...)) 更稳定
log_I_true_result <- log(besselI(x_val, nu_val, expon.scaled = TRUE)) + x_val
cat(sprintf("R 内置函数 log(I_%.1f(%.1f)) 的精确值: %.4f\n", 
            nu_val, x_val, log_I_true_result))









besselI_asymp <- function(v, x) {
  # Large-order asymptotic approximation for I_v(x)
  # valid when v, x are large and z = x/v is fixed (not close to 1)
  
  z <- x / v
  psi <- sqrt(1 + z^2)
  eta <- psi + log(z / (1 + psi))
  
  # log-scale computation to avoid overflow/underflow
  logIv <- -0.5 * log(2 * pi * v) - 0.25 * log(1 + z^2) + v * eta
  
  # raw value (may underflow if extremely small)
  # Iv_approx <- exp(logIv)
  
  return(logIv)
}
besselI_asymp(56, 552.100989303939)
