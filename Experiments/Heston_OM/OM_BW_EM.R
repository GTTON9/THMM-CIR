# ================================================================
# 0. Helper: Log-Sum-Exp for numerical stability
# ================================================================
log_sum_exp <- function(x) {
  max_x <- max(x)
  if (!is.finite(max_x)) return(max_x)
  max_x + log(sum(exp(x - max_x)))
}

# ================================================================
# 1. E-Step: Forward-Backward Algorithm (Log Space)
# ================================================================
#' Run Forward-Backward Algorithm in Log Space
#' @param log_allprobs Matrix [N_states x T] of log emission probs
#' @param log_Gamma Matrix [N_states x N_states] of log transition probs
#' @param log_delta Vector [N_states] of log initial probs
#' @return List containing smoothed probabilities (gamma) and transition pairs (xi)
forward_backward_log <- function(log_allprobs, log_Gamma, log_delta) {
  nstates <- nrow(log_allprobs)
  T_periods <- ncol(log_allprobs)
  
  # --- Forward Pass (alpha) ---
  log_alpha <- matrix(-Inf, nrow = nstates, ncol = T_periods)
  
  # Initialization (t=1)
  log_alpha[, 1] <- log_delta + log_allprobs[, 1]
  
  # Induction
  for (t in 2:T_periods) {
    for (j in 1:nstates) {
      # alpha[t, j] = log_sum_exp(alpha[t-1, i] + log_Gamma[i, j]) + log_prob[t, j]
      log_alpha[j, t] <- log_sum_exp(log_alpha[, t-1] + log_Gamma[, j]) + log_allprobs[j, t]
    }
  }
  
  # --- Backward Pass (beta) ---
  log_beta <- matrix(-Inf, nrow = nstates, ncol = T_periods)
  
  # Initialization (t=T): beta is 1 (log beta is 0)
  log_beta[, T_periods] <- 0
  
  # Induction
  for (t in (T_periods - 1):1) {
    for (i in 1:nstates) {
      # beta[t, i] = log_sum_exp(log_Gamma[i, j] + log_prob[t+1, j] + beta[t+1, j])
      term <- log_Gamma[i, ] + log_allprobs[, t+1] + log_beta[, t+1]
      log_beta[i, t] <- log_sum_exp(term)
    }
  }
  
  # --- Smoothing (Gamma and Xi) ---
  
  # 1. Compute State Probabilities (gamma)
  # log_gamma = alpha + beta - total_log_likelihood
  log_likelihood <- log_sum_exp(log_alpha[, T_periods])
  log_gamma <- log_alpha + log_beta - log_likelihood
  gamma <- exp(log_gamma)
  
  # Normalize gamma numerically to ensure sum to 1 (handling tiny errors)
  gamma <- apply(gamma, 2, function(x) x / sum(x))
  
  # 2. Compute Transition Expected Counts (xi_sum)
  # We need the sum over t=1 to T-1 of xi_{t}(i,j)
  xi_sum <- matrix(0, nrow = nstates, ncol = nstates)
  
  for (t in 1:(T_periods - 1)) {
    for (i in 1:nstates) {
      for (j in 1:nstates) {
        log_xi_t <- log_alpha[i, t] + log_Gamma[i, j] + 
          log_allprobs[j, t+1] + log_beta[j, t+1] - log_likelihood
        xi_sum[i, j] <- xi_sum[i, j] + exp(log_xi_t)
      }
    }
  }
  
  return(list(
    gamma = gamma,         # [N_states x T]
    xi_sum = xi_sum,       # [N_states x N_states] expected transitions
    log_lik = log_likelihood
  ))
}







# ================================================================
# 2. M-Step Helper: Optimize Heston Parameters for one state
# ================================================================

#' Objective function for Heston parameters in M-step
#' Minimizes: Sum_{t} weight_t * OM_Cost(path_t | params)
M_step_objective_one_state <- function(par, weights, result_with_paths, H, N, lambda_om) {
  # Decode params
  kappa <- par[1]
  theta <- par[2]
  sigma <- par[3]
  
  # Constraints check (return huge cost if violated)
  if(kappa <= 0 || theta <= 0 || sigma <= 0) return(1e10)
  
  total_weighted_cost <- 0
  T_periods <- length(weights)
  
  # Loop over time (can be vectorized or done in C++ for speed, but keeping R loop for clarity)
  # To speed up, we only compute for t where weight > 1e-4
  active_indices <- which(weights > 1e-6)
  
  if(length(active_indices) == 0) return(0) # No weight for this state
  
  for (t in active_indices) {
    w <- weights[t]
    
    # Upper Path
    cost_upper <- 0
    if (!is.na(result_with_paths$path_upper_a[t])) {
      val <- calculate_OM_single(
        a = result_with_paths$path_upper_a[t],
        b = result_with_paths$path_upper_b[t],
        kappa = kappa, theta = theta, sigma = sigma,
        H = H, N = N
      )
      if(is.finite(val)) cost_upper <- val else cost_upper <- 1e6
    }
    
    # Lower Path
    cost_lower <- 0
    if (!is.na(result_with_paths$path_lower_a[t])) {
      val <- calculate_OM_single(
        a = result_with_paths$path_lower_a[t],
        b = result_with_paths$path_lower_b[t],
        kappa = kappa, theta = theta, sigma = sigma,
        H = H, N = N
      )
      if(is.finite(val)) cost_lower <- val else cost_lower <- 1e6
    }
    
    # Contribution to Q-function (min cost = max prob)
    total_weighted_cost <- total_weighted_cost + w * (cost_upper + cost_lower)
  }
  
  return(total_weighted_cost)
}



#' Objective function for Heston parameters in M-step
M_step_objective_one_state <- function(par, weights, result_with_paths, H, N, lambda_om) {
  kappa <- par[1]
  theta <- par[2]
  sigma <- par[3]
  
  if(kappa <= 0 || theta <= 0 || sigma <= 0) return(1e10)
  
  total_weighted_cost <- 0
  total_weight_sum <- 0
  
  # 获取非零权重的索引
  active_indices <- which(weights > 1e-6)
  if(length(active_indices) == 0) return(0)
  
  for (t in active_indices) {
    w <- weights[t]
    
    # 1. 计算 OM 能量 (这一部分你已经有了)
    cost_upper <- 0
    if (!is.na(result_with_paths$path_upper_a[t])) {
      val <- calculate_OM_single(..., sigma = sigma, ...) # 你的原始调用
      if(is.finite(val)) cost_upper <- val else cost_upper <- 1e6
    }
    
    cost_lower <- 0
    if (!is.na(result_with_paths$path_lower_a[t])) {
      val <- calculate_OM_single(..., sigma = sigma, ...) # 你的原始调用
      if(is.finite(val)) cost_lower <- val else cost_lower <- 1e6
    }
    
    om_term <- cost_upper + cost_lower
    
    # 2. 关键修复：加入 Normalization 惩罚
    # 对于每个路径点，都需要一个 log(sigma) 惩罚。
    # N 是积分步数。实际上这是对路径测度 Jacobian 的近似。
    # 简单理解：为了让概率密度积分为1，sigma越大，峰值必须越低。
    # 系数 2 是因为我们有 Upper 和 Lower 两条路径
    norm_term <- 2 * log(sigma) 
    
    # 这里的 scale 需要注意：OM通常是积分类似项，norm_term 是对数项。
    # 如果 calculate_OM_single 已经包含了 1/sigma^2，那么这里必须加上 log(sigma)
    # 通常建议形式： Minimize: N * log(sigma) + lambda * OM_Integral
    
    # 组合 Cost: 我们希望最大化概率，即最小化 (-log P)
    # -log P ~ (Normalization) + lambda * (Energy)
    
    step_cost <- norm_term + lambda_om * om_term
    
    total_weighted_cost <- total_weighted_cost + w * step_cost
  }
  
  return(total_weighted_cost)
}


# ================================================================
# 3. Main EM Algorithm Function
# ================================================================

Heston_fit_model_EM <- function(result_with_paths,
                                controls,
                                max_iter = 50,
                                tol = 1e-4,
                                lambda_om = 1.0,
                                H = 5/252,
                                N_integration = 100,
                                seed = NULL,
                                verbose = TRUE) {
  
  set.seed(seed)
  nstates <- controls$states[1]
  T_periods <- nrow(result_with_paths)
  
  # --- Initialization ---
  # Use simple heuristics or random values to start
  # Or use the 'Heston_get_init' logic if available. Here we assume simple start.
  pseudo_data <- list(
    data = result_with_paths$v_start,
    controls = controls
  )
  class(pseudo_data) <- "Heston_data"
  
  initial_values <- Heston_get_init(
    data = pseudo_data,
    ncluster = 1,
    seed = 999,
    verbose = TRUE,
    initial_estimate = NULL
  )

  par <- tryCatch({
    parUncon2par_heston(initial_values[[1]], controls, FALSE, numerical_safeguard = TRUE)
  }, error = function(e) {
    return(NULL)
  })
  print(par)
  # Extract parameters
  Gamma <- par[["Gamma"]]
  kappa <- par[["kappa"]]

  theta <- par[["theta"]]
  sigma <- par[["sigma"]]
  
  
  log_lik_hist <- c()
  
  # EM Loop
  for (iter in 1:max_iter) {
    
    # ==========================
    # 1. E-STEP
    # ==========================
    
    # Calculate initial distribution (stationary)
    delta <- tryCatch({
      solve(t(diag(nstates) - Gamma + 1), rep(1, nstates))
    }, error = function(e) rep(1, nstates)/nstates)
    
    # Calculate Log Emission Probabilities Matrix [N_states x T]
    log_allprobs <- matrix(0, nstates, T_periods)
    
    # Pre-calculate OM costs for current parameters
    # Note: We reuse calculate_OM_for_likelihood but modify it to return log probs directly
    # Or just loop here:
    for (k in 1:nstates) {
      for (t in 1:T_periods) {
        # Calculate cost for upper + lower
        c_up <- 0; c_low <- 0;
        
        if(!is.na(result_with_paths$path_upper_a[t])) {
          v <- calculate_OM_single(result_with_paths$path_upper_a[t], 
                                   result_with_paths$path_upper_b[t], 
                                   kappa[k], theta[k], sigma[k], H, N_integration)
          c_up <- ifelse(is.finite(v), v, 1e6)
        }
        
        if(!is.na(result_with_paths$path_lower_a[t])) {
          v <- calculate_OM_single(result_with_paths$path_lower_a[t], 
                                   result_with_paths$path_lower_b[t], 
                                   kappa[k], theta[k], sigma[k], H, N_integration)
          c_low <- ifelse(is.finite(v), v, 1e6)
        }
        
        # log P = -lambda * (Cost_up + Cost_low)
        log_allprobs[k, t] <- -lambda_om * (c_up + c_low)
      }
    }
 
    # Run Forward-Backward
    log_Gamma <- log(Gamma + 1e-100) # Safety
    log_delta <- log(delta + 1e-100)
    
    fb_res <- forward_backward_log(log_allprobs, log_Gamma, log_delta)
    
    curr_ll <- fb_res$log_lik
    log_lik_hist <- c(log_lik_hist, curr_ll)
    
    if (verbose) {
      cat(sprintf("Iter %d: Log-Likelihood = %.4f\n", iter, curr_ll))
    }
    
    # Check Convergence
    if (iter > 1 && abs(curr_ll - log_lik_hist[iter-1]) < tol) {
      if (verbose) cat("Converged.\n")
      break
    }
    
    # ==========================
    # 2. M-STEP
    # ==========================
    
    # --- Update Transition Matrix (Gamma) ---
    # Gamma_ij = sum_t xi_t(i,j) / sum_t gamma_t(i)
    # xi_sum is sum over t=1..T-1
    # gamma sum should also be over t=1..T-1 for denominator
    
    denom <- rowSums(fb_res$gamma[, 1:(T_periods-1), drop=FALSE])
    new_Gamma <- fb_res$xi_sum / (denom + 1e-10) # Vector recycling
    
    # Normalize rows
    new_Gamma <- sweep(new_Gamma, 1, rowSums(new_Gamma), "/")
    Gamma <- new_Gamma
    
    # --- Update Emission Parameters (kappa, theta, sigma) for each state ---
    for (k in 1:nstates) {
      # Weights for this state
      weights_k <- fb_res$gamma[k, ]
      
      # Starting guess
      p0 <- c(kappa[k], theta[k], sigma[k])
      print(p0)
      # Optimize weighted OM cost
      # We use 'optim' with L-BFGS-B to enforce bounds (avoid negative params!)
      opt_res <- try(optim(
        par = p0,
        fn = M_step_objective_one_state,
        weights = weights_k,
        result_with_paths = result_with_paths,
        H = H, N = N_integration, lambda_om = lambda_om,
        method = "L-BFGS-B",
        lower = c(0.1, 0.01, 0.05), # Bounds: kappa, theta, sigma
        upper = c(50.0, 1.0, 0.8),
        control = list(factr = 1e2) # Loose convergence for inner loop speed
      ), silent = TRUE)
      
      if (!inherits(opt_res, "try-error")) {
        kappa[k] <- opt_res$par[1]
        theta[k] <- opt_res$par[2]
        sigma[k] <- opt_res$par[3]
      }
    }
    
    # Reorder parameters to keep states identifiable (e.g., by theta)
    # This prevents label switching mid-optimization
    ord <- order(theta)
    theta <- theta[ord]
    kappa <- kappa[ord]
    sigma <- sigma[ord]
    Gamma <- Gamma[ord, ord]
    
  } # End EM Loop
  
  # Format output to match previous style
  estimate_list <- list(kappa=kappa, theta=theta, sigma=sigma, Gamma=Gamma)
  
  out <- list(
    estimate = estimate_list, # Not parUncon anymore, but direct values
    ll = tail(log_lik_hist, 1),
    lls = log_lik_hist,
    gamma = fb_res$gamma,
    n_iter = iter
  )
  
  return(out)
}




# ================================================================
# 4. Usage Script
# ================================================================

# 1. 准备数据和环境 (假设之前代码已运行，result_with_paths 已存在)
# ...

# 2. 设置控制参数
controls_em <- Heston_set_controls(
  states = 2,
  sdds = "Heston",
  horizon = nrow(result_with_paths)
)

cat("========================================\n")
cat("Starting EM Algorithm for O-M THMM...\n")
cat("========================================\n")

# 3. 运行 EM 拟合
# 注意：EM 算法对初始值不那么敏感，但也依赖。这里内部使用随机初始化。
# 如果想跑多次取最优，可以在外层加个循环。

start_time <- Sys.time()

model_em <- Heston_fit_model_EM(
  result_with_paths = result_with_paths,
  controls = controls_em,
  max_iter = 50,          # 最大迭代次数
  tol = 1e-4,             # 收敛阈值
  lambda_om = 1.0,        # OM 权重
  H = 5/252,
  N_integration = 100,
  seed = 123,
  verbose = TRUE
)

end_time <- Sys.time()
cat("EM Fitting time:", round(difftime(end_time, start_time, units="mins"), 2), "mins\n")

# 4. 提取结果
cat("\n--- Estimated Parameters (EM) ---\n")
print(model_em$estimate)

# 5. 解码状态 (可以直接使用 EM 输出的 gamma，或者再跑一次 Viterbi)
# Viterbi is for "Most Likely Path", Gamma is for "Posterior Probability"
# Let's decode based on max posterior probability from EM

decoded_states_em <- apply(model_em$gamma, 2, which.max)
plot(decoded_states_em)


Gamma_BW <- model_em$estimate$Gamma

regime_params_BW <- list(
  list(kappa = model_em$estimate$kappa[1], 
       theta = model_em$estimate$theta[1], 
       sigma = model_em$estimate$sigma[1]),
  list(kappa = model_em$estimate$kappa[2], 
       theta = model_em$estimate$theta[2], 
       sigma = model_em$estimate$sigma[2])
)


result_complete <- calculate_OM_batch(result_with_paths, regime_params_BW)
OM_viterbi_result <- viterbi_om_pure(OM_upper= result_complete[, c("OM_upper_R1", "OM_upper_R2")],
                                     OM_lower= result_complete[, c("OM_lower_R1", "OM_lower_R2")],
                                     nstates,
                                     Gamma_BW,
                                     lambda_om = 1.0,
                                     use_upper = TRUE,
                                     use_lower = TRUE)

plot(OM_viterbi_result$states_estimate, col="blue")
lines(Reg_chain_year+1)

plot2 <- plot_viterbi_om(
  batch_results_complete = result_complete,
  nstates = 2,
  Gamma = Gamma_BW,
  kappa = c(regime_params_BW[[1]]$kappa, regime_params_BW[[2]]$kappa),
  theta = c(regime_params_BW[[1]]$theta, regime_params_BW[[2]]$theta),
  sigma = c(regime_params_BW[[1]]$sigma, regime_params_BW[[2]]$sigma),
  Reg_chain = Reg_chain_year ,
  lambda_om = 1,
  normalize_method = "log",
  show_om_contribution = TRUE
)

print(model_em$estimate)


