

get_HMM_param <- function(model){
  parUncon <- model$estimate
  
  
  inv_logit <- function(x) {
    return(exp(x) / (1 + exp(x)))
  }
  
  inv_log <- function(x) {
    return(exp(x))
  }
  
  gamma21_uncon <- parUncon["gammasUncon_21"] 
  gamma12_uncon <- parUncon["gammasUncon_12"] 
  
  P21 <- inv_logit(gamma21_uncon) # P(State 2 -> State 1)
  P12 <- inv_logit(gamma12_uncon) # P(State 1 -> State 2)

  P11 <- 1 - P12
  P22 <- 1 - P21
  
  Gamma <- matrix(c(P11, P21, P12, P22), nrow = 2, byrow = TRUE, 
                  dimnames = list(c("State 1", "State 2"), c("State 1", "State 2")))
  
  
  mu1 <- parUncon["muUncon_1"]
  mu2 <- parUncon["muUncon_2"]
  
  sigma1 <- inv_log(parUncon["sigmaUncon_1"])
  sigma2 <- inv_log(parUncon["sigmaUncon_2"])

  df1 <- inv_log(parUncon["dfUncon_1"])
  df2 <- inv_log(parUncon["dfUncon_2"])
  
  emission_params <- data.frame(
    State = c("State 1", "State 2"),
    Mu = c(mu1, mu2),
    Sigma = c(sigma1, sigma2),
    DF = c(df1, df2)
  )
  
  
  
  df2 <- inv_log(parUncon["dfUncon_2"])
  return(list(
    Gamma = Gamma,
    Emission_Params = emission_params
  ))
}









