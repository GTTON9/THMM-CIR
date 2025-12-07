source("Heston_likelihood.R")
#' Initialization of numerical likelihood optimization
#' 
#' @description
#' This helper function generates a set of initial values for the numerical
#' optimization of the model likelihood function.
#' 
#' @param initial_estimate
#' Optionally defines an initial estimate for the numerical likelihood 
#' optimization. Good initial estimates can improve the optimization process.
#' Can be:
#' - \code{NULL} (the default), in this case
#'   - applies a heuristic to calculate a good initial estimate
#'   - or uses the true parameter values (if available and 
#'     \code{data$controls$origin} is \code{TRUE})
#' - or an object of class \code{parUncon} (i.e., a \code{numeric} of 
#'   unconstrained model parameters), for example the estimate of a 
#'   previously fitted model (i.e. the element \code{model$estimate}). 
#'   
#' @param seed
#' Set a seed for the generation of initial values.
#' No seed by default.
#' 
#' @inheritParams fit_model
#' 
#' @return
#' A \code{list}, where each element is an object of class \code{parUncon}.
#'
#' @keywords internal

Heston_get_init<- function(
    data, ncluster = 1, seed = NULL, verbose = TRUE, initial_estimate = NULL
) {
  
  ### input checks
  checkmate::assert_class(data, "Heston_data")
  checkmate::assert_number(ncluster)
  checkmate::assert_flag(verbose)
  controls <- data[["controls"]]
  checkmate::assert_class(controls, "Heston_controls")
  runs <- controls[["fit"]][["runs"]]
  
  ### set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ### define likelihood function
  target <- ifelse(!data[["controls"]][["hierarchy"]], Heston_nLL_hmm, Heston_nLL_hmm)
  
  ### define function to compute log-likelihood value at initial estimate
  compute_ll_at_initial_estimate <- function(initial_estimate) {
  
    a <-target(
      parUncon = initial_estimate,
      observations = data[["data"]],
      controls = controls
    )
    return(a)
  }
  
  ### define check function for initial estimate
  check_initial_estimate <- function(initial_estimate, verbose, return_value) {

    # this line may use some random paramter to get the length, ignore the output
    # expected_length <- length(par2parUncon_heston(Heston_parameters(controls,), controls))
    expected_length <- suppressMessages(suppressWarnings({
      invisible(capture.output(
        result <- length(par2parUncon_heston(Heston_parameters(controls), controls)),
        type = c("output", "message")
      ))
      
      # 2. 返回结果
      result
    }))
    
    test_initial_estimate <- oeli::test_numeric_vector(
      initial_estimate, finite = TRUE, any.missing = FALSE, len = expected_length
    )
    
    if (!test_initial_estimate) {
      ll <- NA_real_
      if (verbose) {
        
        error_msg <- oeli::check_numeric_vector(
          initial_estimate, finite = TRUE, any.missing = FALSE, len = expected_length
        )
        message("'initial_estimate' is bad: ", error_msg)
      }
      
    } else {
 
      ll <- try(compute_ll_at_initial_estimate(initial_estimate), silent = TRUE) 
     
    }
    
    
    if (!(inherits(ll, "try-error") || is.na(ll) || is.nan(ll) || abs(ll) > 1e100)) {
      if (return_value) {
        return(ll)
      } else {
        return(TRUE)
      }
    } else {
      if (verbose) {
        message("Evaluating the log-likelihood at 'initial_estimate' failed")
      }
      if (return_value) {
        return(NA_real_) 
      } else {
        return(FALSE)
      }
    }
  }
  
  ### initialize at pre-defined initial estimate
  
  initial_heuristic <- function(data, states, positive_mu) {
    
    # initial clustering
    cluster <- stats::kmeans(
      data[!is.na(data)], centers = states, iter.max = 100, nstart = 100
    )$cluster
    print(cluster)
    # Gamma <- matrix(0.1 / (states - 1), nrow = states, ncol = states)
    # diag(Gamma) <- 0.9
    
    
    
    
    Gamma <- matrix(NA, states, states)
    
    for (i in 1:states) {
      indices_i <- which(cluster[1:(length(cluster) - 1)] == i)
      
      if (length(indices_i) == 0) {
        Gamma[i, ] <- 1 / states 
        next
      }
      
      N_i <- length(indices_i)
      numerator_row <- numeric(states) 
      
      for (j in 1:states) {
        N_ij <- sum(cluster[indices_i + 1] == j)
        numerator_row[j] <- max(N_ij, 1e-2) 
      }

      S_i <- sum(numerator_row)
      
      Gamma[i, ] <- numerator_row / S_i 
    }
    
    kappa_est <- numeric(states)
    theta_est <- numeric(states)
    sigma_est <- numeric(states)
    
    
    dt <- 1/252
    
    Delta_Xt <- diff(data)
    Xt_lag <- data[-length(data)]
    
    cluster_lag <- cluster[-length(cluster)]
    
    for (s in seq_len(states)) {
      
      Xt_s <- data[cluster == s]
      
      # theta (numerical protection)
      theta_s <- mean(Xt_s, na.rm = TRUE)
      theta_est[s] <- max(theta_s, 1e-6)
      
      
      
      indices_s <- which(cluster_lag == s)
      
      if (length(indices_s) < 3) {
        warning(paste("State", s, "has too few observations for OLS. Using default kappa = 1.0"))
        kappa_est[s] <- 2
        next 
      }
      
      Delta_Xt_s <- Delta_Xt[indices_s]
      Xt_lag_s <- Xt_lag[indices_s]
      
      # sigma (numerical protection)
      sigma_s_base <- sd(Xt_s, na.rm = TRUE)
      sigma_est[s] = sqrt( var(Delta_Xt_s) / (dt * mean(Xt_lag_s)) )
      
      ols_model_s <- lm(Delta_Xt_s ~ Xt_lag_s)
      
      beta_hat <- coef(ols_model_s)["Xt_lag_s"]
      alpha_hat <- coef(ols_model_s)["(Intercept)"]
      
    
      kappa_s_raw <- -beta_hat / dt
      
      # numerical protection
      kappa_est[s] <- max(kappa_s_raw, 1e-3)
      
      # if(s == 1 ){
      #   kappa_est[s] = 2
      #   theta_est[s] = 0.1
      #   sigma_est[s] = 0.1
      # }else{
      #   kappa_est[s] = 1
      #   theta_est[s] = 0.2
      #   sigma_est[s] = 0.1
      # }
    }
    list(
      "cluster" = cluster,
      "pars" = list(
        "kappa" = kappa_est, 
        "theta" = theta_est, 
        "sigma" = sigma_est, 
        "Gamma" = Gamma
      )
    )
  }
  
  
  
  
  ### applying heuristic
  if (is.null(initial_estimate)) {
    if (controls[["hierarchy"]]) {
      
      ### heuristic for coarse-scale
      heuristic_cs <- initial_heuristic(
        data = data[["data"]][, 1], states = controls[["states"]][1],
        positive_mu = controls[["sdds"]][[1]][["name"]] %in% c("gamma", "poisson")
      )
      cluster_cs <- heuristic_cs[["cluster"]]
      initial_estimate_list_cs <- heuristic_cs[["pars"]]
      
      ### heuristic for fine-scale
      Gamma_star <- list()
      mu_star <- list()
      sigma_star <- list()
      for (s in seq_len(controls[["states"]][1])) {
        cluster_fs_data <- as.vector(data[["data"]][s == cluster_cs, -1])
        
        heuristic_fs <- initial_heuristic(
          cluster_fs_data, states = controls[["states"]][2],
          positive_mu = controls[["sdds"]][[2]][["name"]] %in% c("gamma", "poisson")
        )
        
        Gamma_star[[s]] <- heuristic_fs[["pars"]][["Gamma"]]
        mu_star[[s]] <- heuristic_fs[["pars"]][["mu"]]
        sigma_star[[s]] <- heuristic_fs[["pars"]][["sigma"]]
      }
      initial_estimate_list_fs <- list(
        "Gamma_star" = Gamma_star, "mu_star" = mu_star, "sigma_star" = sigma_star
      )
      
      ### combine coarse-scale and fine-scale
     
      initial_estimate <- par2parUncon_heston(
        do.call(
          what = Heston_parameters,
          args = c(
            initial_estimate_list_cs, initial_estimate_list_fs, list(controls)
          )
        ),
        controls
      )
    } else {
      
      initial_estimate_list <- initial_heuristic(
        data = data[["data"]], states = controls[["states"]],
        positive_mu = FALSE 
      )[["pars"]]
      
      print("Initial Estimate")
      print(initial_estimate_list)
      
      initial_estimate <- par2parUncon_heston(
        do.call(
          what = Heston_parameters,
          args = c(initial_estimate_list, list(controls))
        ),
        controls
      )

    }
  }

  ### jitter 'initial_estimate'
  jitter_initial_estimate <- function(initial_estimate, N) {
    par_names <- names(initial_estimate)
    jittered <- matrix(
      initial_estimate, nrow = N, ncol = length(initial_estimate), byrow = TRUE
    )
    parameter_class <- gsub("_.*", "", names(initial_estimate))
    for (class in unique(parameter_class)) {
      id <- which(parameter_class == class)
      jittered[, id] <- jitter(jittered[, id], factor = 10)
    }

    lapply(seq_len(N), function(i) {
      structure(jittered[i, ], names = par_names, class = c("parUncon", "numeric"))
    })
  }
  
  
  ### generate and check start values
  if (verbose) {
    message("Checking start values...")
  }
  N <- runs * 2
  
  initial_values <- jitter_initial_estimate(initial_estimate, N)

  initial_values[[1]] <- initial_estimate
  ll_at_initial_values <- rep(NA_real_, N)
  
  while (TRUE) {
    
    ind <- which(is.na(ll_at_initial_values))
    N <- length(initial_values)
    
    ### stopping criterium
    if (length(ind) == 0 && N == runs) {
      break
    }
    
    ### catch cases, where parallelization does not make sense
    if (length(ind) < ncluster) {
      ncluster <- max(1, min(ncluster, length(ind)))
      
    }
   
    ### compute log-likelihood values
    if (length(ind) > 0) {
      if (ncluster <= 1) {
        
        for (i in ind) {
        
          ll_at_initial_values[i] <- check_initial_estimate(# this will call Heston_parameters
            initial_values[[i]], verbose = FALSE, return_value = TRUE)
   

        }
        
      } else {
        
        
        cluster <- parallel::makeCluster(ncluster)
        doSNOW::registerDoSNOW(cluster)
        
        ll_at_initial_values_update <- foreach::foreach(
          N = seq_along(ind), .packages = "fHMM"
        ) %dopar% {
          check_initial_estimate(
            initial_values[[ind[N]]], verbose = FALSE, return_value = TRUE
          )
        }
        
        parallel::stopCluster(cluster)
        ll_at_initial_values[ind] <- ll_at_initial_values_update
      }
      
    }

    ### replace initial values that lead to NA
    ind <- which(is.na(ll_at_initial_values))
    
    if (length(ind) > 0) {
      initial_values[ind] <- jitter_initial_estimate(initial_estimate, length(ind))
    }

    ### drop largest negative log-likelihood value
    largest <- which.max(ll_at_initial_values)
    if (length(largest) >= 1) {
      ll_at_initial_values <- ll_at_initial_values[-largest]
      initial_values <- initial_values[-largest]
    }
    
   
  }
  
  ### return list of initial values

  return(initial_values)
}








