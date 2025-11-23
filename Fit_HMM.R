
library(fHMM)
fit_HMM <- function(series, real_reg = NULL, sdds = "t"){
  set.seed(999)
  series_length <- length(series)
  start_date <- as.Date("2024-01-01")
  date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
  
  
  
  my_data_df <- data.frame(
    Date = date_sequence,
    price = seriess
  )
  
  # followed the example
  series_control <- set_controls( 
    states      = 2,     # 2 state
    sdds        = sdds,         
    date_column = "Date",
    file        = my_data_df, 
    data_column = "price",      
    logreturns  = TRUE,         
    from        = date_sequence[1],             
    to          = date_sequence[length(date_sequence)]
  )
  
  
  data_hmm <- prepare_data(series_control)

  model_hmm <- fit_model(data_hmm) # fit HMM
  model_hmm <- decode_states(model_hmm)
  
  
  plot(model_hmm$decoding - 1,          
            col = "blue",     
            xlab = "Time Step",  
            ylab = "Regime",           
            main = "Estimated Regime Chain",
            ylim = range(c(series$decoding - 1, real_reg))
  )
  
  if(! is.null(real_reg)){
    
    lines(real_reg,
          col = "red",     
          type = "l",
          lwd = 3
          
    )
    legend("topleft",              
           legend = c("Estimated Regime", "True Regime"),
           col = c("blue", "red"),
           lty = c(1, 2),
           cex = 0.5
    )
  }
  return(model_hmm)
}

  


