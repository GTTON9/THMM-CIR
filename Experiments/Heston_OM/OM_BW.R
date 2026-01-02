

# ================================================================
# O-M THMM Complete Implementation for Heston Parameter Fitting
# ================================================================

# Source your existing files first:
# source("Heston_get_init.R")
# source("fHMM_model.R")
# source("Heston_likelihood.R")

# ================================================================
# 1. Calculate O-M Emission Probabilities
# ================================================================

#' Calculate O-M Functional for Likelihood
#'
#' @param result_with_paths Data frame with path parameters (from construct_linear_paths)
#' @param kappa Vector of mean reversion speeds for each state
#' @param theta Vector of long-term means for each state
#' @param sigma Vector of volatilities for each state
#' @param H Time horizon (default 5/252)
#' @param N Integration steps (default 100)
#' @param lambda_om Weight for O-M contribution (default 1.0)
#' @param use_upper Use upper path (default TRUE)
#' @param use_lower Use lower path (default TRUE)
#'
#' @return Matrix [nstates x T] of log emission probabilities

calculate_OM_for_likelihood <- function(result_with_paths,
                                        kappa,
                                        theta,
                                        sigma,
                                        H = 5/252,
                                        N = 100,
                                        lambda_om = 1.0,
                                        use_upper = TRUE,
                                        use_lower = TRUE) {
  
  T_periods <- nrow(result_with_paths)

  nstates <- length(kappa)
  # Initialize emission probability matrix
  allprobs <- matrix(0, nrow = nstates, ncol = T_periods)
  
  for (t in 1:T_periods) {

    # Extract path parameters for time t

    a_upper <- result_with_paths$path_upper_a[t]
    b_upper <- result_with_paths$path_upper_b[t]
    a_lower <- result_with_paths$path_lower_a[t]
    b_lower <- result_with_paths$path_lower_b[t]

    for (i in 1:nstates) {
      
      # Calculate O-M functional for each state
      om_cost <- 0
      
      # Upper path contribution
      if (use_upper && !is.na(a_upper) && !is.na(b_upper)) {
        om_upper <- calculate_OM_single(
          a = a_upper,
          b = b_upper,
          kappa = kappa[i],
          theta = theta[i],
          sigma = sigma[i],
          H = H,
          N = N
        )
        # a_val <- a_upper; b_val <- b_upper; 

        
       
        
        
        if (is.finite(om_upper)) {
          om_cost <- om_cost + om_upper
        } else {

          om_cost <- om_cost + 1e6  # Penalty for invalid paths
        }
      }
      
      # Lower path contribution
      if (use_lower && !is.na(a_lower) && !is.na(b_lower)) {
        om_lower <- calculate_OM_single(
          a = a_lower,
          b = b_lower,
          kappa = kappa[i],
          theta = theta[i],
          sigma = sigma[i],
          H = H,
          N = N
        )
        
        
        if (is.finite(om_lower)) {
          om_cost <- om_cost + om_lower
        } else {
          om_cost <- om_cost + 1e6
        }
      }
      
      # Emission log-probability: log P(φ_t | state=i) = -λ × L[φ_t]
      allprobs[i, t] <- -lambda_om * om_cost

    }
  }

  
  return(allprobs)
}


calculate_OM_single <- function(a, b, kappa, theta, sigma, 
                                H = 5/252, N = 100) {
  # 1. 基础合法性检查
  # 如果起点就在 0 附近，或者终点 a+b*H 导致负数，直接判定为极低概率路径

  
  
  end_val <- a + b * H

  # 2. 向量化准备
  delta_tau <- H / N
  tau <- seq(0, H, length.out = N + 1)
  phi <- a + b * tau

  # 3. 核心计算
  phi_dot <- b
  drift <- kappa * (theta - phi)

  # 使用一个稍微大一点的 epsilon 或者加上一个小常数防止除零
  # 也可以考虑使用 pmax(phi, 1e-5)
  denom <- sigma^2 * phi


  
  f <- (phi_dot - drift)^2 / denom
  
  # 4. 梯形积分
  integral <- (sum(f) - 0.5 * (f[1] + f[N + 1])) * delta_tau
  L <- 0.5 * integral
  return(L)
}



# ================================================================
# 2. O-M Negative Log-Likelihood Function
# ================================================================

#' Negative Log-Likelihood Function for O-M THMM
#'
#' @param parUncon Unconstrained parameter vector
#' @param result_with_paths Path data (from construct_linear_paths)
#' @param controls Heston controls object
#' @param lambda_om Weight for O-M functional (default 1.0)
#' @param H Time horizon (default 5/252)
#' @param N_integration Integration steps (default 100)
#'
#' @return Negative log-likelihood value

Heston_nLL_hmm_om <- function(parUncon,
                              result_with_paths,
                              controls,
                              lambda_om = 1.0,
                              H = 5/252,
                              N_integration = 100) {
  
  # Set class for parameter vector
  class(parUncon) <- "parUncon"
  
  T_periods <- nrow(result_with_paths)
  nstates <- controls[["states"]][1]
  
  # Convert unconstrained parameters to constrained parameters
  par <- tryCatch({
    parUncon2par_heston(parUncon, controls, FALSE, numerical_safeguard = TRUE)
  }, error = function(e) {
    return(NULL)
  })
  
  # Extract parameters
  Gamma <- par[["Gamma"]]
  kappa <- par[["kappa"]]
  theta <- par[["theta"]]
  sigma <- par[["sigma"]]
  
  # Calculate initial distribution (stationary distribution)
  delta <- tryCatch({
    solve(t(diag(nstates) - Gamma + 1), rep(1, nstates))
  }, error = function(e) {
    rep(1, nstates) / nstates
  })
  
  # Calculate emission probabilities using O-M functional
  allprobs <- tryCatch({
    calculate_OM_for_likelihood(
      result_with_paths = result_with_paths,
      kappa = kappa,
      theta = theta,
      sigma = sigma,
      H = H,
      N = N_integration,
      lambda_om = lambda_om,
      use_upper = TRUE,
      use_lower = TRUE
    )
  }, error = function(e) {
    return(NULL)
  })
  
  
  
  print("-------------------------------------------------")
  print(Gamma)
  print(theta)
  print(kappa)
  print(sigma)
  print(allprobs)
  
  
  # Check if emission calculation failed
  if (is.null(allprobs) || any(!is.finite(allprobs))) {
    return(1e10)
  }

  # Calculate log-likelihood using forward algorithm
  ll <- - LL_HMM_R(allprobs, Gamma, delta)
  print(ll)
  # Return negative log-likelihood
  if (!is.finite(ll) || is.na(ll)) {
    return(1e10)
  }

  return(ll)
}

# ================================================================
# 3. O-M Model Fitting Function
# ================================================================

#' Fit Heston Model using O-M THMM
#'
#' @param result_with_paths Path data from construct_linear_paths
#' @param controls Heston controls (from Heston_set_controls)
#' @param fit Fitting options
#' @param runs Number of optimization runs
#' @param lambda_om Weight for O-M functional
#' @param H Time horizon for paths
#' @param N_integration Integration steps for O-M
#' @param ncluster Number of clusters for parallel processing
#' @param seed Random seed
#' @param verbose Print progress messages
#' @param initial_estimate Initial parameter estimate
#'
#' @return fHMM_model object with fitted parameters

Heston_fit_model_om <- function(result_with_paths,
                                controls,
                                fit = list(),
                                runs = 10,
                                origin = FALSE,
                                accept = 1:3,
                                gradtol = 0.01,
                                iterlim = 1000,
                                print.level = 0,
                                steptol = 0.01,
                                lambda_om = 1.0,
                                H = 5/252,
                                N_integration = 100,
                                ncluster = 1,
                                seed = NULL,
                                verbose = TRUE,
                                initial_estimate = NULL) {
  
  set.seed(seed)
  
  # Input validation
  if (!checkmate::test_count(ncluster, positive = TRUE)) {
    stop("'ncluster' must be a positive integer.", call. = FALSE)
  }
  
  if (!checkmate::test_flag(verbose)) {
    stop("'verbose' must be either TRUE or FALSE.", call. = FALSE)
  }
  
  # Update controls with fitting options
  controls <- Heston_set_controls(
    controls = controls,
    hierarchy = controls[["hierarchy"]],
    states = controls[["states"]],
    sdds = controls[["sdds"]],
    horizon = nrow(result_with_paths),
    period = controls[["period"]],
    data = controls[["data"]],
    fit = fit,
    runs = runs,
    origin = origin,
    accept = accept,
    gradtol = gradtol,
    iterlim = iterlim,
    print.level = print.level,
    steptol = steptol
  )
  
  # Create pseudo data object for initialization
  # Use variance data for initial parameter estimation
  pseudo_data <- list(
    data = result_with_paths$v_start,
    controls = controls
  )
  class(pseudo_data) <- "Heston_data"
  # Get initial values
  if (verbose) {
    message("Generating initial values for O-M THMM...")
  }
  
  initial_values <- Heston_get_init(
    data = pseudo_data,
    ncluster = ncluster,
    seed = seed,
    verbose = verbose,
    initial_estimate = initial_estimate
  )
  
  runs <- length(initial_values)
  
  # Define target function
  target <- function(parUncon) {
    Heston_nLL_hmm_om(
      parUncon = parUncon,
      result_with_paths = result_with_paths,
      controls = controls,
      lambda_om = lambda_om,
      H = H,
      N_integration = N_integration
    )
  }
  
  # Progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent, :eta ETA",
      total = runs,
      width = 45,
      clear = TRUE,
      complete = "=",
      incomplete = "-",
      current = ">"
    )
    pb$message("Maximizing O-M likelihood...")
  }
  
  # Optimization
  start_time <- Sys.time()
  
  if (ncluster == 1) {
    mods <- list()
    
    for (run in seq_len(runs)) {
      if (verbose) pb$tick(0)
      
      suppressWarnings({
        # mod <- try(
        #   stats::nlm(
        #     f = target,
        #     p = initial_values[[run]],
        #     iterlim = controls[["fit"]][["iterlim"]],
        #     steptol = controls[["fit"]][["steptol"]],
        #     gradtol = controls[["fit"]][["gradtol"]],
        #     print.level = 2,# controls[["fit"]][["print.level"]],
        #     hessian = FALSE
        #   ),
        # 
        #   silent = TRUE
        # )
        
        mod <- try({
          optim(
            par = p0,
            fn = target,          # 这里的 target 已经包含了 result_with_paths 等参数
            method = "L-BFGS-B",
            lower = c(0.1, 0.01, 0.05), # 顺序必须对应 parUncon 中的 kappa, theta, sigma
            upper = c(50.0, 1.0, 0.2),
            control = list(
              factr = 1e9,        # 迭代精度限制
              pgtol = 0, 
              maxit = controls[["fit"]][["iterlim"]]
            ),
            hessian = FALSE
          )
        }, silent = TRUE)
        
        opt_res <- try(optim(
        par = p0,
        fn = M_step_objective_one_state,
        weights = weights_k,
        result_with_paths = result_with_paths,
        H = H, N = N_integration, lambda_om = lambda_om,
        method = "L-BFGS-B",
        lower = c(0.1, 0.01, 0.05), # Bounds: kappa, theta, sigma
        upper = c(50.0, 1.0, 0.2),
        control = list(factr = 1e9) # Loose convergence for inner loop speed
      ), silent = TRUE)
      })
      
      accept_run <- !inherits(mod, "try-error") &&
        mod[["code"]] %in% controls[["fit"]][["accept"]]
      
      if (accept_run) {
        mods[[run]] <- mod
      } else {
        mods[[run]] <- NA
      }
      
      if (verbose) pb$tick()
    }
  }
  
  end_time <- Sys.time()
  
  # Extract log-likelihoods
  lls <- -unlist(sapply(mods, `[`, "minimum"), use.names = FALSE)
  
  if (all(is.na(lls))) {
    stop("None of the estimation runs ended successfully.\n",
         "Adapt 'accept' or increase 'runs' in 'controls'.",
         call. = FALSE)
  }
  
  # Approximate Hessian
  if (verbose) {
    message("Approximating Hessian...")
  }
  
  best_estimate <- mods[[which.max(lls)]][["estimate"]]
  
  fisher <- pracma::hessdiag(
    f = target,
    x = best_estimate
  )
  
  if (all(fisher > 0)) {
    inverse_fisher <- 1/fisher
  } else {
    hessian <- suppressWarnings(
      pracma::hessian(f = target, x0 = best_estimate)
    )
    inverse_fisher <- diag(MASS::ginv(hessian))
  }
  
  if (verbose) {
    message("Fitting completed!")
  }
  
  # Best model
  mod <- mods[[which.max(lls)]]
  ll <- -mod[["minimum"]]
  estimate <- mod[["estimate"]]
  class(estimate) <- "parUncon"
  
  estimation_time <- ceiling(difftime(end_time, start_time, units = "mins"))
  
  # Create model object
  out <- fHMM_model(
    data = pseudo_data,
    estimate = estimate,
    nlm_output = mod,
    estimation_time = estimation_time,
    ll = ll,
    lls = lls,
    gradient = mod$gradient,
    inverse_fisher = inverse_fisher,
    decoding = NULL
  )
  
  # Reorder states by mean
  out <- Heston_reorder_states(out, state_order = "mean")
  
  return(out)
}

# ================================================================
# 4. Decode States using O-M THMM
# ================================================================

#' Decode States for O-M THMM
#'
#' @param model Fitted model from Heston_fit_model_om
#' @param result_with_paths Path data
#' @param lambda_om Weight for O-M functional
#' @param H Time horizon
#' @param N_integration Integration steps
#'
#' @return Model with decoded states

decode_states_heston_om <- function(model,
                                    result_with_paths,
                                    lambda_om = 1.0,
                                    H = 5/252,
                                    N_integration = 100) {
  
  # Extract parameters
  par <- parUncon2par_heston(
    model$estimate,
    model$data$controls,
    FALSE,
    numerical_safeguard = TRUE
  )
  
  # Calculate emission probabilities
  allprobs <- calculate_OM_for_likelihood(
    result_with_paths = result_with_paths,
    kappa = par$kappa,
    theta = par$theta,
    sigma = par$sigma,
    H = H,
    N = N_integration,
    lambda_om = lambda_om
  )
  
  # Run Viterbi
  viterbi_result <- viterbi_om_pure(
    OM_upper = matrix(allprobs[1, ], ncol = 2),  # Dummy structure
    OM_lower = matrix(allprobs[2, ], ncol = 2),
    nstates = length(par$kappa),
    Gamma = par$Gamma,
    lambda_om = 1.0
  )
  
  # Add decoding to model
  model$decoding <- viterbi_result$path
  
  return(model)
}

# ================================================================
# 5. Complete Usage Example
# ================================================================

#' Complete O-M THMM Analysis Workflow
#'
#' @examples
#' # Step 1: Simulate or load data
#' S_daily <- ...  # Daily stock prices
#' V_daily <- ...  # Daily variance
#' 
#' # Step 2: Prepare paths
#' batch_results <- calculate_expected_moves_batch(S_daily, V_daily)
#' batch_results_with_IV <- calculate_IV_with_dividend_skew(
#'   batch_results, r = 0.02, H = 5/252
#' )
#' result_with_paths <- construct_linear_paths(batch_results_with_IV, H = 5/252)
#' 
#' # Step 3: Set controls
#' controls <- Heston_set_controls(
#'   states = 2,
#'   sdds = "Heston",
#'   horizon = nrow(result_with_paths),
#'   runs = 10
#' )
#' 
#' # Step 4: Fit model
#' model_om <- Heston_fit_model_om(
#'   result_with_paths = result_with_paths,
#'   controls = controls,
#'   runs = 10,
#'   lambda_om = 1.0,
#'   H = 5/252,
#'   N_integration = 100,
#'   seed = 999,
#'   verbose = TRUE
#' )
#' 
#' # Step 5: Extract parameters
#' param_om <- parUncon2par_heston(
#'   model_om$estimate,
#'   controls,
#'   FALSE,
#'   numerical_safeguard = TRUE
#' )
#' 
#' print(param_om$kappa)
#' print(param_om$theta)
#' print(param_om$sigma)
#' print(param_om$Gamma)
#' 
#' # Step 6: Decode states
#' model_om <- decode_states_heston_om(
#'   model_om,
#'   result_with_paths,
#'   lambda_om = 1.0
#' )
#' 
#' # Step 7: Compare with true regime
#' accuracy <- mean(model_om$decoding == (Reg_chain_year + 1))
#' cat(sprintf("Detection accuracy: %.2f%%\n", accuracy * 100))

run_om_thmm_workflow <- function(S_daily,
                                 V_daily,
                                 Reg_chain_year = NULL,
                                 true_params = NULL,
                                 Gamma_true = NULL) {
  
  cat("\n")
  cat("=" %R% 70)
  cat("\nO-M THMM COMPLETE WORKFLOW\n")
  cat("=" %R% 70)
  cat("\n\n")
  
  # Step 1: Calculate expected moves
  cat("Step 1: Calculating expected moves...\n")
  batch_results <- calculate_expected_moves_batch(S_daily, V_daily)
  
  # Step 2: Calculate IV
  cat("Step 2: Calculating implied volatilities...\n")
  batch_results_with_IV <- calculate_IV_with_dividend_skew(
    batch_results,
    r = 0.02,
    H = 5/252
  )
  
  # Step 3: Construct paths
  cat("Step 3: Constructing linear paths...\n")
  result_with_paths <- construct_linear_paths(batch_results_with_IV, H = 5/252)
  
  # Step 4: Set controls
  cat("Step 4: Setting up controls...\n")
  controls <- Heston_set_controls(
    states = 2,
    sdds = "Heston",
    horizon = nrow(result_with_paths),
    runs = 10
  )
  
  # Step 5: Fit model
  cat("Step 5: Fitting O-M THMM model...\n\n")
  model_om <- Heston_fit_model_om(
    result_with_paths = result_with_paths,
    controls = controls,
    runs = 10,
    lambda_om = 1.0,
    H = 5/252,
    N_integration = 100,
    seed = 999,
    verbose = TRUE
  )
  
  # Step 6: Extract parameters
  param_om <- parUncon2par_heston(
    model_om$estimate,
    controls,
    FALSE,
    numerical_safeguard = TRUE
  )
  
  # Step 7: Decode states
  model_om <- decode_states_heston_om(
    model_om,
    result_with_paths,
    lambda_om = 1.0
  )
  
  # Step 8: Results
  cat("\n")
  cat("=" %R% 70)
  cat("\nRESULTS\n")
  cat("=" %R% 70)
  cat("\n\n")
  
  cat("Estimated Parameters:\n")
  for (i in 1:2) {
    cat(sprintf("  Regime %d: κ=%.4f, θ=%.4f, σ=%.4f\n",
                i, param_om$kappa[i], param_om$theta[i], param_om$sigma[i]))
  }
  
  cat("\nTransition Matrix:\n")
  print(round(param_om$Gamma, 4))
  
  if (!is.null(Reg_chain_year)) {
    accuracy <- mean(model_om$decoding == (Reg_chain_year + 1))
    cat(sprintf("\nRegime Detection Accuracy: %.2f%%\n", accuracy * 100))
  }
  
  return(list(
    model = model_om,
    parameters = param_om,
    result_with_paths = result_with_paths
  ))
}

































# ================================================================
# METHOD 2: O-M THMM (New method)
# ================================================================



# # Step 1: Calculate expected moves
# cat("Calculating expected moves...\n")
# batch_results <- calculate_expected_moves_batch(S_daily, V_daily)
# 
# # Step 2: Calculate IV
# cat("Calculating implied volatilities...\n")
# batch_results_with_IV <- calculate_IV_with_dividend_skew(
#   batch_results,
#   r = 0.02,
#   H = 5/252,
#   d_call = 0,
#   d_put = 0.03,
#   d_iv = 0
# )
# 
# # Step 3: Construct linear paths
# cat("Constructing linear paths...\n")
# result_with_paths <- construct_linear_paths(batch_results_with_IV, H = 5/252)

# Step 4: Set controls for O-M THMM
cat("Setting up O-M THMM controls...\n")
controls_om <- Heston_set_controls(
  states = 2,
  sdds = "Heston",
  horizon = nrow(result_with_paths),
  runs = 10
)

# Step 5: Fit O-M THMM model
cat("Fitting O-M THMM model...\n\n")




model_om <- Heston_fit_model_om(
  result_with_paths = result_with_paths,
  controls = controls_om,
  runs = 100,
  lambda_om = 1.0,
  H = 5/252,
  N_integration = 100,
  seed = 999,
  verbose = TRUE
)




# Step 6: Extract parameters
param_om <- parUncon2par_heston(
  model_om$estimate,
  controls_om,
  FALSE,
  numerical_safeguard = TRUE
)
param_om
# Step 7: Decode states
model_om <- decode_states_heston_om(
  model_om,
  result_with_paths,
  lambda_om = 100,
  H = 5/252,
  N_integration = 100
)

plot(model_om$decoding[1,])
states_estimate_om <- model_om$decoding

cat("\nO-M THMM Results:\n")
cat("κ:", param_om$kappa, "\n")
cat("θ:", param_om$theta, "\n")
cat("σ:", param_om$sigma, "\n")
cat("Γ:\n")
print(param_om$Gamma)

accuracy_om <- mean(states_estimate_om == (Reg_chain_year + 1))
cat(sprintf("\nAccuracy: %.2f%%\n", accuracy_om * 100))

# ================================================================
# COMPARISON AND VISUALIZATION
# ================================================================

cat("\n")
cat("#" %R% 70)
cat("\nCOMPARISON: Traditional HMM vs O-M THMM\n")
cat("#" %R% 70)
cat("\n\n")

# Create comparison table
comparison_df <- data.frame(
  Parameter = c("κ₁", "θ₁", "σ₁", "P₁₁", "P₁₂",
                "κ₂", "θ₂", "σ₂", "P₂₁", "P₂₂"),
  True = c(
    kappa[1], theta[1], sigma[1], Gamma[1,1], Gamma[1,2],
    kappa[2], theta[2], sigma[2], Gamma[2,1], Gamma[2,2]
  ),
  Traditional_HMM = c(
    param_traditional$kappa[1],
    param_traditional$theta[1],
    param_traditional$sigma[1],
    param_traditional$Gamma[1,1],
    param_traditional$Gamma[1,2],
    param_traditional$kappa[2],
    param_traditional$theta[2],
    param_traditional$sigma[2],
    param_traditional$Gamma[2,1],
    param_traditional$Gamma[2,2]
  ),
  OM_THMM = c(
    param_om$kappa[1],
    param_om$theta[1],
    param_om$sigma[1],
    param_om$Gamma[1,1],
    param_om$Gamma[1,2],
    param_om$kappa[2],
    param_om$theta[2],
    param_om$sigma[2],
    param_om$Gamma[2,1],
    param_om$Gamma[2,2]
  )
)

# Calculate errors
comparison_df$Error_Traditional <- abs(comparison_df$Traditional_HMM - comparison_df$True)
comparison_df$Error_OM <- abs(comparison_df$OM_THMM - comparison_df$True)
comparison_df$Pct_Error_Traditional <- 
  comparison_df$Error_Traditional / comparison_df$True * 100
comparison_df$Pct_Error_OM <- 
  comparison_df$Error_OM / comparison_df$True * 100

print(comparison_df)

# Summary statistics
cat("\n\nSummary:\n")
cat(sprintf("Traditional HMM Accuracy: %.2f%%\n", accuracy_traditional * 100))
cat(sprintf("O-M THMM Accuracy: %.2f%%\n", accuracy_om * 100))
cat(sprintf("Traditional Mean |Error|: %.4f\n", mean(comparison_df$Error_Traditional)))
cat(sprintf("O-M THMM Mean |Error|: %.4f\n", mean(comparison_df$Error_OM)))

# ================================================================
# VISUALIZATION
# ================================================================

param_om




