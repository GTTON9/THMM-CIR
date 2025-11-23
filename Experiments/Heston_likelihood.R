

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

Heston_nLL_hmm <- function(parUncon, observations, controls) {

  class(parUncon) <- "parUncon"
  T <- length(observations)

  nstates <- controls[["states"]][1]

  par <- parUncon2par_heston(parUncon, controls, FALSE, numerical_safeguard = TRUE)
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

  allprobs <- matrix(NA_real_, nstates, T-1)
  for (i in 1:nstates) {
    if (sdd == "Heston") {
        
        allprobs[i, ] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
        # allprobs[i, ] <- get_transition_density_Gaussain(observations, kappa[i], theta[i], sigma[i])
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

get_transition_density_Gaussain <- function(observations, kappa, theta, sigma){
  
  T_obs <- length(observations)
  
  all_ln_densities <- numeric(T_obs - 1) 
  
  for (t in 1:(T_obs - 1)) {
    
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1] # k=1
    transition_density <- get_gaussian_approx_density(
      V_t, V_t_plus_k, 
      kappa = kappa, theta = theta, sigma = sigma, 
      dt = 1/252)
    
    all_ln_densities[t] <- transition_density
  }
  
  # print(all_ln_densities)
  return(all_ln_densities)
}


get_transition_density_heston_ln <- function(observations, kappa, theta, sigma) {
  T_obs <- length(observations)
  
  all_ln_densities <- numeric(T_obs - 1) 
  
  for (t in 1:(T_obs - 1)) {
 
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1] # k=1
    
    ln_transition_density <- ln_d_Heston(
      V_t = V_t, 
      V_t_plus_k = V_t_plus_k, 
      k = 1, kappa = kappa, theta = theta, sigma = sigma
    )
    
    all_ln_densities[t] <- exp(ln_transition_density)
  }

 #  print(all_ln_densities)
  return(all_ln_densities)
}




ln_d_Heston<- function(V_t, V_t_plus_k, k = 1, kappa, theta, sigma,
                               min_val = 1e-300, z_large_threshold = 50) {
  # numeric safety
  sigma <- max(sigma, 1e-12)
  exp_k <- exp(-kappa * k)
  den <- (1 - exp_k) * (sigma^2)
  if (abs(den) < 1e-16) return(-Inf)
  
  C <- 2 * kappa / den
  q <- (2 * kappa * theta) / (sigma^2) - 1
  
  # ensure nonnegative V inputs
  V_t <- pmax(V_t, 0)
  V_tp <- pmax(V_t_plus_k, 0)
  
  u <- C * V_t * exp_k
  v <- C * V_tp
  
  # clamp to avoid exact zeros causing log(0)
  u <- pmax(u, min_val)
  v <- pmax(v, min_val)
  
  sqrt_u <- sqrt(u)
  sqrt_v <- sqrt(v)
  z <- 2 * sqrt_u * sqrt_v  # 2*sqrt(uv)
  
  # compute log(v/u) robustly
  log_vu_ratio <- 0
  if (u > 0 && v > 0) {
    # use log difference (more stable)
    log_vu_ratio <- log(v) - log(u)
  }
  
  ln_Iq <- NA_real_
  
  # Try scaled bessel first (preferred if reliable)
  success_scaled <- FALSE
  if (is.finite(z) && z >= 0) {
    # Protect bessel call in tryCatch
    scaled_val <- tryCatch(besselI(z, nu = q, expon.scaled = TRUE),
                           error = function(e) NA_real_,
                           warning = function(w) NA_real_)
    if (!is.na(scaled_val) && is.finite(scaled_val) && scaled_val > 0) {
      ln_Iq <- z + log(scaled_val)  # log(I_q(z)) = z + log(I_q(z)*exp(-z))
      success_scaled <- TRUE
    }
  }
  
  # If scaled bessel is not good, use asymptotics depending on regime
  if (!success_scaled) {
    # decide regime: large z (argument dominates), large order, or uniform Olver
    if (z > z_large_threshold) {
      # large-argument asymptotic (z >> q case also lands here)
      # log I_q(z) ≈ z - 0.5*log(2*pi*z) + small correction
      ln_Iq <- z - 0.5 * log(2 * pi * pmax(z, 1e-300))
      
    } else if (abs(q) > 50 && (z / abs(q)) < 0.8) {
      # large order, small ratio (use large-nu expansion)
      q_abs <- abs(q)
      ln_Iq <- q_abs * (1 + log(pmax(z/(2*q_abs), 1e-300))) - 0.5 * log(2 * pi * q_abs)
      
    } else {
      # uniform Olver-type expansion (safe middle-ground)
      q_abs <- max(abs(q), 1e-12)
      t <- z / q_abs
      r <- sqrt(1 + t^2)
      eta <- r + log( pmax(t, 1e-300) ) - log(1 + r)
      ln_Iq <- q_abs * eta - 0.5 * log(2 * pi * q_abs * r)
    }
  }
  
  # Now assemble log-density but avoid catastrophic cancellation
  # Original: ln = log(C) - u - v + (q/2) * log(v/u) + ln_Iq
  # Use identity: -u - v + z = - (sqrt(u) - sqrt(v))^2
  # If ln_Iq used asymptotic with leading term z, separate that out.
  
  # Determine if ln_Iq was formed with explicit 'z' included.
  # We will compute the combination -u - v + ln_Iq stably:
  
  # If ln_Iq uses leading z term (our asymptotic choices include it),
  # compute residual = ln_Iq - z and then use identity for -u - v + z.
  # If scaled bessel was used, ln_Iq already contains z; same logic applies.
  
  residual <- ln_Iq - z  # may be small negative (e.g., -0.5*log(2*pi*z))
  # compute stable quadratic term:
  quad_term <- - (sqrt_u - sqrt_v)^2
  
  # assemble:
  ln_density <- log(abs(C)) + quad_term + residual + (q / 2) * log_vu_ratio
  
  # If any NaN/Inf appear, fall back to a safe large-negative value:
  if (!is.finite(ln_density)) ln_density <- -1e300
  
  return(ln_density)
}
