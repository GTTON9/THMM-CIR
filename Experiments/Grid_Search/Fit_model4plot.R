

#' Negative log-likelihood function of an HMM
#'
#' @description
#' This function computes the negative log-likelihood of an HMM.
#'
#' @param parUncon
#' An object of class \code{parUncon}.
#' @param observations
#' The vector of the simulated or empirical data used for estimation.
#' @param controls
#' An object of class \code{Heston_controls}.
#'
#' @return
#' The negative log-likelihood value.
#'
#' @keywords internal

# install.packages("Bessel")
# library(Bessel)
Heston_nLL_hmm <- function(parUncon, observations, controls) {
  
  class(parUncon) <- "parUncon"
  T <- length(observations)
  
  nstates <- controls[["states"]][1]
  
  par <- parUncon2par_heston(parUncon, controls, FALSE, numerical_safeguard = TRUE)
  Heston_trace$par_list[[length(Heston_trace$par_list) + 1]] <- par
  
  # print("----- likelihood par-----")
  # print(par)
  sdd <- controls[["sdds"]][[1]]$name
  
  Gamma <- par[["Gamma"]]
  
  delta <- try(
    solve(t(diag(nstates) - Gamma + 1), rep(1, nstates)),
    silent = TRUE
  )
  
  if (inherits(delta, "try-error")) {
    delta <- rep(1, nstates) / nstates
  }
  
  kappa <- par[["kappa"]]
  theta <- par[["theta"]]
  sigma <- par[["sigma"]]
  # kappa <- c(10,5)
  
  allprobs <- matrix(NA_real_, nstates, T-1)
  for (i in 1:nstates) {
    if (sdd == "Heston") {
      
      allprobs[i, ] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
    } 
    
    else {
      stop("Unknown state-dependent distribution", call. = FALSE)
    }
  }
  # print(allprobs)
  
  ll <- - LL_HMM_R(allprobs, Gamma, delta)
  
  # print(paste("negative_log_likelihood:",ll))
  
  return(ll)
  # -LL_HMM_Rcpp(allprobs, Gamma, delta, nstates, T)
  
}
get_transition_density_Gaussian <- function(observations, kappa, theta, sigma){
  T_obs <- length(observations)
  all_ln_densities <- numeric(T_obs - 1) 
  
  for (t in 1:(T_obs - 1)) {
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1]
    transition_density <- get_gaussian_approx_density(
      V_t, V_t_plus_k, 
      kappa = kappa, theta = theta, sigma = sigma, 
      dt = 1/252)
    
    all_ln_densities[t] <- transition_density
  }
  return(all_ln_densities)
}


get_transition_density_heston_ln <- function(observations, kappa, theta, sigma) {
  T_obs <- length(observations)
  
  all_ln_densities <- numeric(T_obs - 1) 
  
  for (t in 1:(T_obs - 1)) {
    
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1]
    
    ln_transition_density <- ln_d_Heston(
      V_t = V_t, 
      V_t_plus_k = V_t_plus_k, 
      k = 1/252, kappa = kappa, theta = theta, sigma = sigma
    )
    
    all_ln_densities[t] <- (ln_transition_density)
  }
  
  return(all_ln_densities)
}


# 
# ln_d_Heston <- function(V_t, V_t_plus_k, k = 1/252, kappa, theta, sigma) {
#   
#   # --- 1. 参数计算 ---
#   sigma <- max(sigma, 1e-8)
#   exp_k <- exp(-kappa * k)
#   den_factor <- (1 - exp_k) * (sigma^2)
#   
#   if(abs(den_factor) < 1e-10) return(-Inf) 
#   
#   C <- 2 * kappa / den_factor
#   q <- (2 * kappa * theta) / (sigma^2) - 1 # 
#   u <- C * V_t * exp_k
#   v <- C * V_t_plus_k
#   z_arg <- 2 * sqrt(u * v) 
#   
#   bessel_val <- tryCatch({
#     BesselI(z_arg, nu = q, expon.scaled = TRUE) 
#   }, error = function(e) {
#     besselI(z_arg, nu = q, expon.scaled = TRUE)
#   })
#   
#   
#   if (is.finite(bessel_val) && bessel_val > 0) {
#     ln_I_q <- z_arg + log(bessel_val)
#     
#   } else {
#     
#     
#     q_abs <- abs(q)
#     
#     if (z_arg > 100 && q_abs < z_arg / 5) {
#       # besselIasym for large z regular q
#       ln_I_q <- tryCatch({
#         besselIasym(z_arg, nu = q, log = TRUE)
#       }, error = function(e) {
#         # DLMF 10.40.1: I_q(z) ≈ exp(z) / sqrt(2*pi*z) (z -> ∞)
#         warning(paste("besselIasym failed for z =", z_arg, 
#                       ". Switching to large-z asymptotic approximation (DLMF 10.40.1). Original error:", 
#                       conditionMessage(e)))
#         z_arg - 0.5 * log(2 * pi * z_arg)
#       })
#       
#       
#     } else if (q_abs > 10 && q_abs > z_arg * 1.2) {
#       # for large q regular z
#       # besselI.nuAsym (nu >= 0)
#       ln_I_q <- tryCatch({
#         besselI.nuAsym(z_arg, nu = q_abs, log = TRUE, k.max = 5)
#       }, error = function(e) {
#         # q * (1 + log(z / (2 * q))) - 0.5 * log(2 * pi * q)
#         # DLMF 10.41.1
#         warning(paste("besselI.nuAsym failed for q =", q_abs, 
#                       ". Switching to asymptotic approximation (DLMF 10.41.1). Original error:", 
#                       conditionMessage(e)))
#         
#         q_abs * (1 + log(z_arg / (2 * q_abs))) - 0.5 * log(2 * pi * q_abs)
#       })
#       
#     } else {
#       # ratio_abs <- z_arg / q_abs
#       # r <- sqrt(1 + ratio_abs^2)
#       # # eta = r + log(ratio_abs / (1 + r))
#       # eta <- r + log(ratio_abs) - log(1 + r)
#       # 
#       # # ln_I_q ≈ |q| * eta - 0.5 * log(2 * pi * |q| * r)
#       # ln_I_q <- q_abs * eta - 0.5 * log(2 * pi * q_abs * r)
#       # 
#       # 
#       # 
#       mu_approx <- V_t + kappa * (theta - V_t) * k
# 
#       # Variance: Var[V_{t+k} | V_t] ≈ diffusion^2 * k = (sigma^2 * V_t) * k
#       variance_approx <- (sigma^2) * V_t * k
# 
#       # Ensure standard deviation is positive
#       sd_approx <- sqrt(pmax(variance_approx, 1e-8))
# 
#       # Calculate the log density using the Normal distribution
#       ln_density <- dnorm(V_t_plus_k, mean = mu_approx, sd = sd_approx, log = TRUE)
#       return(ln_density)
#     }
#   }
# 
#   log_vu_ratio <- 0
#   if(u > 1e-10 && v > 1e-10) {
#     log_vu_ratio <- log(v / u)
#   }
#   
#   ln_density <- log(C) - u - v + (q / 2) * log_vu_ratio + ln_I_q
#   
#   return(ln_density)
# }



ln_d_Heston <- function(V_t, V_t_plus_k, 
                        k = 1/252, 
                        kappa, theta, sigma,
                        min_val = 1e-300) {
  # ---- 1. Compute Constants ----
  ek <- exp(-kappa * k)
  
  c  <- (sigma^2 * (1 - ek)) / (4 * kappa)
  d  <- 4 * kappa * theta / sigma^2
  lambda <- (4 * kappa * ek / (sigma^2 * (1 - ek))) * V_t
  
  # Avoid numerical blowups
  c      <- max(c,  min_val)
  d      <- max(d,  min_val)
  lambda <- max(lambda, min_val)
  V_ratio <- V_t_plus_k / c
  
  # ---- 2. Compute noncentral chi-square density ----
  dens <- dchisq(V_ratio, df = d, ncp = lambda, log = FALSE)
  
  # Adjust for variable transformation (divide by c)
  return(log(dens / c))
}




















source("Heston_get_init.R")
source("fHMM_model.R")
#' Model fitting
#'
#' @description
#' This function fits a hidden Markov model via numerical likelihood 
#' maximization.
#'
#' @details
#' Multiple optimization runs starting from different initial values are 
#' computed in parallel if \code{ncluster > 1}.
#'
#' @param data
#' An object of class \code{\link{fHMM_data}}.
#' 
#' @param ncluster
#' Set the number of clusters for parallel optimization runs to reduce 
#' optimization time.
#' By default, \code{ncluster = 1} (no clustering).
#' 
#' @param verbose
#' Set to \code{TRUE} to print progress messages.
#' 
#' @inheritParams Heston_set_controls
#' @inheritParams get_initial_values
#' @inheritParams Heston_set_controls
#'
#' @return
#' An object of class \code{\link{fHMM_model}}.
#' 
#' @examples
#' ### 2-state HMM with normal distributions
#' 
#' # set specifications
#' controls <- Heston_set_controls(
#'   states = 2, sdds = "normal", horizon = 100, runs = 10
#' )
#' 
#' # define parameters
#' parameters <- fHMM_parameters(controls, mu = c(-1, 1), seed = 1)
#' 
#' # sample data
#' data <- prepare_data(controls, true_parameter = parameters, seed = 1)
#' 
#' # fit model
#' model <- fit_model(data, seed = 1)
#' 
#' # inspect fit
#' summary(model)
#' plot(model, "sdds")
#' 
#' # decode states
#' model <- decode_states(model)
#' plot(model, "ts")
#' 
#' # predict
#' predict(model, ahead = 5)
#'
#' @export

Heston_fit_model <- function (data, controls = data[["controls"]], fit = list(), 
                              runs = 10, origin = FALSE, accept = 1:3, gradtol = 0.01, 
                              iterlim = 1000, print.level = 0, steptol = 0.01, ncluster = 1, 
                              seed = NULL, verbose = TRUE, initial_estimate = NULL) 
{ 
  
  if (!inherits(data, "Heston_data")) {
    stop("'data' is not of class 'Heston_data'.", call. = FALSE)
  }
  if (!checkmate::test_count(ncluster, positive = TRUE)) {
    stop("'ncluster' must be a positive integer.", call. = FALSE)
  }
  if (!checkmate::test_flag(verbose)) {
    stop("'verbose' must be either TRUE or FALSE.", call. = FALSE)
  }
  
  data[["controls"]] <- Heston_set_controls(controls = controls, hierarchy = data[["controls"]][["hierarchy"]], 
                                            states = data[["controls"]][["states"]], sdds = data[["controls"]][["sdds"]], 
                                            horizon = data[["controls"]][["horizon"]], period = data[["controls"]][["period"]], 
                                            data = data[["controls"]][["data"]], fit = fit, runs = runs, 
                                            origin = origin, accept = accept, gradtol = gradtol, 
                                            iterlim = iterlim, print.level = print.level, steptol = steptol)
  
  initial_values <- Heston_get_init(data = data, ncluster = ncluster, # problem
                                    seed = seed, verbose = verbose, initial_estimate = initial_estimate)
  
  # check nstate match the true regime
  runs <- length(initial_values)
  
  target <- Heston_nLL_hmm
  
  if (verbose) {
    pb <- progress::progress_bar$new(format = "[:bar] :percent, :eta ETA", 
                                     total = runs, width = 45, clear = TRUE, complete = "=", 
                                     incomplete = "-", current = ">")
    pb$message("Maximizing likelihood...")
  }
  
  
  
  
  start_time <- Sys.time()
  if (ncluster == 1) {
    mods <- list()
    for (run in seq_len(runs)) {
      if (verbose) 
        pb$tick(0)
      suppressWarnings({
        
        parameter_history <- list(
          par_uncon = list(),
          nll_value = numeric()
        )
        
        mod <- try(
          stats::nlm(f = target, p = initial_values[[run]], 
                     observations = data[["data"]], controls = data[["controls"]], 
                     iterlim = data[["controls"]][["fit"]][["iterlim"]],  # maximum number of iterations
                     steptol = data[["controls"]][["fit"]][["steptol"]], 
                     gradtol = data[["controls"]][["fit"]][["gradtol"]], 
                     print.level = data[["controls"]][["fit"]][["print.level"]], 
                     hessian = FALSE), silent = FALSE) # silent off
        
      })
      
      accept_run <- !inherits(mod, "try-error") && mod[["code"]] %in% 
        data[["controls"]][["fit"]][["accept"]]
      
      if (accept_run) {
        mods[[run]] <- mod
      }
      else {
        mods[[run]] <- NA
      }
      if (verbose) 
        pb$tick()
    }
  }
  
  end_time <- Sys.time()
  lls <- -unlist(sapply(mods, `[`, "minimum"), use.names = FALSE)
  if (all(is.na(lls))) {
    stop("None of the estimation runs ended successfully.\n", 
         "Adapt 'accept' or increase 'runs' in 'controls'.", 
         call. = FALSE)
  }
  if (verbose) 
    message("Approximating Hessian...")
  fisher <- pracma::hessdiag(f = target, x = mods[[which.max(lls)]][["estimate"]], 
                             observations = data[["data"]], controls = data[["controls"]])
  if (all(fisher > 0)) {
    inverse_fisher <- 1/fisher
  }
  else {
    hessian <- suppressWarnings(pracma::hessian(f = target, 
                                                x0 = mods[[which.max(lls)]][["estimate"]], observations = data[["data"]], 
                                                controls = data[["controls"]]))
    inverse_fisher <- diag(MASS::ginv(hessian))
  }
  if (verbose) 
    message("Fitting completed!")
  mod <- mods[[which.max(lls)]]
  ll <- -mod[["minimum"]]
  estimate <- mod[["estimate"]]
  class(estimate) <- "parUncon"
  estimation_time <- ceiling(difftime(end_time, start_time, 
                                      units = "mins"))
  out <- fHMM_model(data = data, estimate = estimate, nlm_output = mod, 
                    estimation_time = estimation_time, ll = ll, lls = lls, 
                    gradient = mod$gradient, inverse_fisher = inverse_fisher, 
                    decoding = NULL)
  
  out <- Heston_reorder_states(out, state_order = "mean")
  return(out)
}











Heston_trace <- new.env()
Heston_trace$par_list <- list()


BW_result <- Gen_fit(Gamma, mu, kappa, theta, sigma, rho, V_ = TRUE , plot_path = TRUE)



library(reshape2)
library(ggplot2)
library(patchwork) # 用于排列图表（如果您想分别创建图表对象）
true_params <- list(kappa = c(10, 5), theta = c(0.03, 0.6), sigma = c(0.05, 0.05))












y_ranges <- list(
  kappa1 = c(0, 15),
  kappa2 = c(0, 15),
  theta1 = c(0, 1),
  theta2 = c(0, 1),
  sigma1 = c(0, 0.1),
  sigma2 = c(0, 0.1)
)
plot_single_param <- function(df, param_name, true_val = NULL) {
  
  df_sub <- df[df$Parameter == param_name, ]
  
  ggplot(df_sub, aes(x = Iteration, y = Value)) +
    geom_line(linewidth = 0.8) +
    { if(!is.null(true_val)) geom_hline(yintercept = true_val,
                                        color = "red", linetype = "dashed") } +
    ylim(y_ranges[[param_name]]) +
    labs(
      title = paste("Trace:", param_name),
      x = "Iteration",
      y = "Value"
    ) +
    theme_minimal(base_size = 14)
}
p_kappa1 <- plot_single_param(trace_long, "kappa1", true_params$kappa[1])
p_kappa2 <- plot_single_param(trace_long, "kappa2", true_params$kappa[2])

p_theta1 <- plot_single_param(trace_long, "theta1", true_params$theta[1])
p_theta2 <- plot_single_param(trace_long, "theta2", true_params$theta[2])

p_sigma1 <- plot_single_param(trace_long, "sigma1", true_params$sigma[1])
p_sigma2 <- plot_single_param(trace_long, "sigma2", true_params$sigma[2])
print(p_kappa1)
print(p_kappa2)

print(p_theta1)
print(p_theta2)

print(p_sigma1)
print(p_sigma2)



