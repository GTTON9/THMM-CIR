#' Log-likelihood function of an (H)HMM
#'
#' @description
#' This function computes the log-likelihood value of a (hierarchical) hidden 
#' Markov model for given observations and parameter values.
#'
#' @inheritParams Heston_parameters
#' @inheritParams parameter_transformations
#' @param observations
#' A \code{numeric} \code{vector} of time-series data.
#' 
#' In the hierarchical case (\code{hierarchy = TRUE}), a \code{matrix} with 
#' coarse-scale data in the first column and corresponding fine-scale data in 
#' the rows.
#' @inheritParams set_controls
#' @param negative
#' Either \code{TRUE} to return the negative log-likelihood value (useful for
#' optimization) or \code{FALSE} (default), else.
#'
#' @return
#' The (negative) log-likelihood value.
#'
#' @examples
#' ### HMM log-likelihood 
#' controls <- set_controls(states = 2, sdds = "normal")
#' parameters <- Heston_parameters(controls)
#' parUncon <- par2parUncon(parameters, controls)
#' observations <- 1:10
#' ll_hmm(parUncon, observations, controls)
#' 
#' ### HHMM log-likelihood 
#' controls <- set_controls(
#'   hierarchy = TRUE, states = c(2, 2), sdds = c("normal", "normal")
#' )
#' parameters <- Heston_parameters(controls)
#' parUncon <- par2parUncon(parameters, controls)
#' observations <- matrix(dnorm(110), ncol = 11, nrow = 10)
#' ll_hmm(parUncon, observations, controls)
#'
#' @export

ll_hmm <- function(
  parUncon,
  observations,
  controls = list(),
  hierarchy = FALSE,
  states = if (!hierarchy) 2 else c(2, 2),
  sdds = if (!hierarchy) "normal" else c("normal", "normal"), 
  negative = FALSE,
  check_controls = TRUE
) {

  if (isTRUE(check_controls)) {
    controls <- set_controls(
      controls = controls, hierarchy = hierarchy, states = states, sdds = sdds
    )
  } else {
    controls <- structure(
      oeli::merge_lists(
        controls,
        list("hierarchy" = hierarchy, "states" = states, "sdds" = sdds)
      ),
      class = "Heston_controls"
    )
  }
  nll <- if (isTRUE(controls[["hierarchy"]])) {
    nLL_hhmm(
      parUncon = parUncon, observations = observations, controls = controls
    )
  } else {
    nLL_hmm(
      parUncon = parUncon, observations = observations, controls = controls
    )
  }
  ifelse(isTRUE(negative), nll, -nll)
}

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

    } 
    
    else {
      stop("Unknown state-dependent distribution", call. = FALSE)
    }
  }

  ll <- -LL_HMM_R(allprobs, Gamma, delta)

  # print(paste("log_likelihood:",ll))
  # fileConn<-file("output.txt")
  # writeLines(paste("log_likelihood:",ll), append = TRUE)

  return(ll)
  # -LL_HMM_Rcpp(allprobs, Gamma, delta, nstates, T)

}







get_transition_density_heston <- function(observations, kappa, theta, sigma) {

  T_obs <- length(observations)
  
  all_densities <- numeric(T_obs - 1) 

  for (t in 1:(T_obs - 1)) {
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1] # k=1
    
      
    if (V_t <= 0 || V_t_plus_k <= 0) {
      
      all_densities[t] <- 1e-300 
      next
    }
    likelihood <- d_Heston(
      V_t = V_t, 
      V_t_plus_k = V_t_plus_k, 
      k = 1, kappa = kappa, theta = theta, sigma = sigma
    )
    
    all_densities[t] <- likelihood
  }
  # print(all_densities)
  return(all_densities) 
}



d_Heston <- function(V_t, V_t_plus_k, k = 1, kappa, theta, sigma) {
 
  C_numerator <- 2 * kappa
  C_denominator <- (1 - exp(-kappa * k)) * (sigma^2)
  C <- C_numerator / C_denominator
 
  # print(paste("C:",C))
  q_numerator <- 2 * kappa * theta
  q_denominator <- sigma^2
  q <- (q_numerator / q_denominator) - 1

  u <- C * V_t * exp(-kappa * k)
  v <- C * V_t_plus_k


  if (u <= 0 || v <= 0) {
    print("badbadbadbadbadbadbad")
    return(1e-300) 
  }
  
  argument <- 2 * sqrt(u * v)
  
  print(paste("q:",q))
  print(paste("argument:", argument))
  
  I_q <- besselI(argument, nu = q)
  print(paste("I_q:",I_q))
  density_factor <- (v / u)^(q / 2)
  
  transition_density <- C * exp(-u - v) * density_factor * I_q


  if (is.nan(transition_density) || transition_density < 0) {
    #print("bad")
 
  }
  print(transition_density)
  return(transition_density)
}


  

ln_d_Heston <- function(V_t, V_t_plus_k, k = 1, kappa, theta, sigma) {
  
  C_numerator <- 2 * kappa
  C_denominator <- (1 - exp(-kappa * k)) * (sigma^2)
  
  
  C <- C_numerator / C_denominator
  
  q_numerator <- 2 * kappa * theta
  q_denominator <- sigma^2

  q <- (q_numerator / q_denominator) - 1
  
  u <- C * V_t * exp(-kappa * k)
  v <- C * V_t_plus_k
  
  argument <- 2 * sqrt(u * v)

  # e^(-x)I_v(x), so take log term
  scaled_I_q <- besselI(argument, nu = q, expon.scaled = TRUE)
  # log(I_q)
  ln_I_q <- -1*argument + log(scaled_I_q)
  
  # ln(d_Heston) = log(C) - u - v + (q / 2) * log(v / u) + ln_I_q
  ln_transition_density <- log(C) - u - v + (q / 2) * log(v / u) + ln_I_q
  
  return(ln_transition_density)
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
  
  
  #print(all_ln_densities)
  return(all_ln_densities)
}

