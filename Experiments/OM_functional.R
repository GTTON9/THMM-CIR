Heston_prepare_test_function <- function(
    observations, 
    vix_data = NULL,
    method = c("vix", "realized_ma", "constant"),
    smooth_bandwidth = 5
) {
  
  method <- match.arg(method)
  T <- length(observations)
  
  if (method == "vix") {
    # VIX作为市场benchmark
    if (is.null(vix_data)) {
      stop("vix_data must be provided when method='vix'")
    }
    
    # VIX通常是百分比形式的年化波动率，转换为方差
    vix_variance <- (vix_data / 100)^2
    
    # 对齐长度
    if (length(vix_variance) != T) {
      warning("VIX data length differs from observations. Truncating to match.")
      vix_variance <- vix_variance[1:min(T, length(vix_variance))]
      if (length(vix_variance) < T) {
        # 如果VIX数据不够长，用最后一个值填充
        vix_variance <- c(vix_variance, rep(tail(vix_variance, 1), T - length(vix_variance)))
      }
    }
    
    # 平滑处理
    phi <- stats::filter(
      vix_variance, 
      rep(1/smooth_bandwidth, smooth_bandwidth), 
      sides = 2
    )
    
    # 处理NA (边界效应)
    phi[is.na(phi)] <- vix_variance[is.na(phi)]
    
  } else if (method == "realized_ma") {
    # 用观测数据自身的移动平均
    
    # 滚动窗口计算实现方差的移动平均
    phi <- stats::filter(
      observations, 
      rep(1/smooth_bandwidth, smooth_bandwidth), 
      sides = 2
    )
    
    # 处理NA
    phi[is.na(phi)] <- observations[is.na(phi)]
    
  } else if (method == "constant") {
    # 简单的常数测试函数（用于测试）
    phi <- rep(mean(observations, na.rm = TRUE), T)
  }
  
  # 数值保护：确保phi严格为正
  phi <- pmax(phi, 1e-6)
  
  return(as.numeric(phi))
}









#' Compute Onsager-Machlup functional for CIR/Heston variance process
#'
#' @description
#' Computes the OM functional value for a given path segment.
#' This quantifies how "likely" the path is under the regime's dynamics.
#'
#' Lower I[φ] → more probable path
#' Higher I[φ] → less probable path
#'
#' @param phi_segment Test function values over the time segment [t-1, t]
#' @param kappa Mean reversion speed
#' @param theta Long-run mean (variance)
#' @param sigma Volatility of volatility (vol-of-vol)
#' @param dt Time step (default 1/252 for daily data)
#'
#' @return Scalar: OM functional value I[φ]
#'
#' @details
#' For CIR process: dv = κ(θ - v)dt + ξ√v dW
#' 
#' OM functional: I[φ] = (1/2) ∫ [(φ̇ - κ(θ - φ))² / (ξ²φ)] dt - (1/2)κT
#'
#' @keywords internal

compute_OM_functional_segment <- function(
    phi_segment,      # 测试函数在[t-1, t]的值
    kappa,           
    theta,           
    sigma,           
    dt = 1/252
) {
  
  # 数值保护
  kappa <- max(kappa, 1e-6)
  theta <- max(theta, 1e-6)
  sigma <- max(sigma, 1e-6)
  phi_segment <- pmax(phi_segment, 1e-6)
  
  # 计算 φ 的导数（数值微分）
  n_points <- length(phi_segment)
  
  if (n_points < 2) {
    # 如果只有一个点，无法计算导数
    return(0)
  }
  
  # 用中心差分或前向差分
  phi_dot <- diff(phi_segment) / dt
  
  # 对应的phi值（取中点）
  phi_mid <- (phi_segment[-n_points] + phi_segment[-1]) / 2
  
  # CIR drift: b(φ) = κ(θ - φ)
  drift <- kappa * (theta - phi_mid)
  
  # CIR diffusion: σ²(φ) = ξ² * φ
  diffusion_sq <- sigma^2 * phi_mid
  
  # OM integrand: [(φ̇ - drift)²] / diffusion²
  integrand <- ((phi_dot - drift)^2) / diffusion_sq
  
  # 第一项：路径偏离项
  term1 <- 0.5 * sum(integrand) * dt
  
  # 第二项：散度修正
  # 对CIR: ∇·b = -κ
  T_total <- (n_points - 1) * dt
  term2 <- -0.5 * kappa * T_total
  
  # 总OM functional
  I <- term1 + term2
  
  # 数值保护：防止极端值
  I <- max(I, -1e10)
  I <- min(I, 1e10)
  
  return(I)
}


#' Compute OM functional for entire time series
#'
#' @description
#' Wrapper function to compute OM functional for all time segments
#'
#' @param phi_test Complete test function φ(t)
#' @param kappa,theta,sigma Regime parameters
#' @param dt Time step
#'
#' @return Vector of OM functional values for each time transition
#'
#' @keywords internal

compute_OM_functional_full <- function(
    phi_test,
    kappa,
    theta,
    sigma,
    dt = 1/252
) {
  
  T <- length(phi_test)
  OM_values <- numeric(T - 1)
  
  for (t in 1:(T - 1)) {
    # 取相邻两个点
    phi_segment <- phi_test[t:(t + 1)]
    
    OM_values[t] <- compute_OM_functional_segment(
      phi_segment = phi_segment,
      kappa = kappa,
      theta = theta,
      sigma = sigma,
      dt = dt
    )
  }
  
  return(OM_values)
}







library(quantmod)
vix <- getSymbols("^VIX", src = "yahoo", auto.assign = FALSE,
                  from = "2000-01-01")  # adjust date range as needed
head(vix)



