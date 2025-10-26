simulate_Reg <- function(
    series_length = 250,
    initial_state = 0,
    Reg_tran = matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2),
    plot = T,
    seed = 999
) {
  set.seed(seed)
  Reg_chain <- numeric(series_length)
  Reg_chain[1] <- initial_state
  
  for (t in 2:series_length) {
    current_state <- Reg_chain[t-1]
    matrix_row_index <- current_state + 1 
    transition_probabilities <- Reg_tran[matrix_row_index, ]
    
    next_state <- sample(
      x = c(0, 1), 
      size = 1, 
      prob = transition_probabilities
    )
    Reg_chain[t] <- next_state
  }
  
  if (plot) {
    plot(Reg_chain, type = 's', 
         ylim = c(-0.1, 1.1),
         main = "Simulated Regime Chain (0=Calm, 1=Turbulent)",
         xlab = "Time Step", 
         ylab = "Regime",
         yaxt = 'n' 
    )
    axis(2, at = c(0, 1), labels = c("Calm", "Turbulent"))
  }
  
  invisible(Reg_chain) 
}
