Heston_nLL_hmm <- function(parUncon, observations, controls) {
  
  class(parUncon) <- "parUncon"
  T <- length(observations)
  nstates <- controls[["states"]][1]
  
  par <- parUncon2par_heston(parUncon, controls, FALSE, numerical_safeguard = TRUE)
  
  
  if (exists("parameter_history", envir = .GlobalEnv)) {
    assign("parameter_history", 
           within(get("parameter_history", envir = .GlobalEnv), {
             # 记录 nlm 正在优化的向量 parUncon
             par_uncon[[length(par_uncon) + 1]] <- parUncon 
           }), 
           envir = .GlobalEnv)
  }
  
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
  
  # 记录 NLL 值 (ll)
  if (exists("parameter_history", envir = .GlobalEnv)) {
    assign("parameter_history", 
           within(get("parameter_history", envir = .GlobalEnv), {
             nll_value <- c(nll_value, ll)
           }), 
           envir = .GlobalEnv)
  }
  
  return(ll)
}






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
      
      
      parameter_history <- list(
        par_uncon = list(),
        nll_value = numeric())
        
        
      suppressWarnings({
        mod <- try(
          stats::nlm(f = target, p = initial_values[[run]], 
                     observations = data[["data"]], controls = data[["controls"]], 
                     iterlim = data[["controls"]][["fit"]][["iterlim"]],  # maximum number of iterations
                     steptol = data[["controls"]][["fit"]][["steptol"]], 
                     gradtol = data[["controls"]][["fit"]][["gradtol"]], 
                     print.level = data[["controls"]][["fit"]][["print.level"]], 
                     hessian = FALSE), silent = FALSE) # silent off
      })
      
      
      par_matrix <- do.call(rbind, parameter_history$par_uncon)
      num_iterations <- nrow(par_matrix)
      
      
      if (num_iterations > 0) {
        # 创建迭代次数向量
        iteration_index <- 1:num_iterations
        
        # 依赖库
        if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
        if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'")
        
        # 将矩阵转换为长格式数据框 (方便 ggplot2 绘图)
        plot_data <- as.data.frame(par_matrix)
        plot_data$Iteration <- iteration_index
        
        # 将列名 V1, V2...改为 Par_1, Par_2...
        param_names <- paste0("Par_", 1:ncol(par_matrix))
        colnames(plot_data)[1:ncol(par_matrix)] <- param_names
        
        plot_long <- tidyr::pivot_longer(plot_data, 
                                         cols = tidyr::starts_with("Par_"), 
                                         names_to = "Parameter", 
                                         values_to = "Value")
        
        # 将真实值转换为长格式数据框
        # 确保真实参数向量 real_par_uncon 存在且长度正确
        if (length(real_par_uncon) == ncol(par_matrix)) {
          real_values_long <- data.frame(
            Parameter = param_names,
            Real_Value = real_par_uncon
          )
          
          # 绘图代码
          p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = Iteration, y = Value, color = Parameter)) +
            ggplot2::geom_line() +
            # 添加真实值水平线
            ggplot2::geom_hline(data = real_values_long, 
                                ggplot2::aes(yintercept = Real_Value, color = Parameter), 
                                linetype = "dashed", 
                                alpha = 0.7) +
            ggplot2::facet_wrap(~ Parameter, scales = "free_y") +
            ggplot2::labs(title = paste0("HMM Parameter Convergence (Run ", run, ")"),
                          y = "Parameter Value (Unconstrained)",
                          x = "Optimization Iteration") +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none")
          
          print(p)
          
        } else {
          warning("Cannot plot convergence: Real (unconstrained) parameter vector is missing or has incorrect length.")
        }
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





model_hmm <- Heston_fit_model(data_hmm) 

final_model <- decode_states_heston(model_hmm) 
states_estimate <- final_model$decoding






