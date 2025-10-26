
library(fHMM)
library(depmixS4 )
source("simulate_heston.R")
set.seed(9999)



simulate_Reg <- function(
    series_length = 252,
    initial_state = 0,
    Reg_tran = matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2),
    plot = T
) {
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


# simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M, method = "E")




# fit the data using depmix
get_reg <- function(Reg_chain, series){
  # Reg_chain: simulated true regime
  # series: data data 
  
  returns <- diff(log(series)) 
  plot(returns, type="l", xlab='', ylab="Returns")
  
  hmm <- depmix(returns ~ 1, family = gaussian(), nstates = 2, data=data.frame(returns=returns))
  hmmfit <- fit(hmm, verbose = FALSE)
  
  post_probs <- posterior(hmmfit)
  layout(1:3)
  plot(Reg_chain, type='s', main='True Regimes', xlab='', ylab='Regime')
  plot(post_probs$state, type='s', main='Detected Regimes', xlab='', ylab='Regime')
  matplot(post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability')
  legend(x='topright', c('Calm','Turbulent'), fill=1:2, bty='n')
  
}


parameter_perm <- matrix(
  c(
    0.05, 2.0, 0.1, 0.2, -0.1,  0.05, 2.0, 0.1, 0.05, -0.1, # sigma
    0.05, 2.0, 1, 0.05, -0.1,  0.05, 2.0, 0.5, 0.05, -0.1,  # theta
    0.05, 2.0, 0.1, 0.05, -0.1,  0.05, 1.0, 0.1, 0.05, -0.1,  # kappa
    0.7, 2.0, 0.1, 0.05, -0.1,  0.05, 2.0, 0.1, 0.05, -0.1# mu
  ),
  nrow = 4,
  ncol = 10,
  byrow = TRUE
)

colnames(parameter_perm) <- c(
  "mu0", "kappa0", "theta0", "sigma0", "rho0",
  "mu1", "kappa1", "theta1", "sigma1", "rho1"
)

S0 <- 100.0; v0 <- 0.04; T <- 1.0; N <- 252; M <- 1;


series_mat <- matrix(NA, nrow = 5, ncol = 253)
set.seed(2999)
Reg_chain <- simulate_Reg(series_length = 250)
plot(Reg_chain)
for( i in 4:4){
  set.seed(100)
  Reg_param <- matrix(parameter_perm[i,], nrow = 2, byrow = T)
  syn_series <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M, method = "M")
  series_mat[i,] <- syn_series$S_paths[1,]
  matplot(syn_series$S_paths[1,], type = "l", 
          main = "Heston Model Asset Price Paths", 
          xlab = "Time Step", ylab = "Asset Price S_t")
  
  S_simulated <- syn_series$S_paths[1,]
  get_reg(Reg_chain, S_simulated)
}
plot(series_mat[1,],type = "s")



mu_0 <- 0.15
kappa_0 <- 2.0
theta_0 <- 0.03
sigma_0 <- 0.15
rho_0 <- -0.1

# Regime 1 (Turbulent): high volatility
mu_1 <- 0.00
kappa_1 <- 1.0
theta_1 <- 0.20
sigma_1 <- 0.6
rho_1 <- -0.8

Reg_param_optimal <- matrix(
  c(
    mu_0, kappa_0, theta_0, sigma_0, rho_0,  # Regime 0 (Calm)
    mu_1, kappa_1, theta_1, sigma_1, rho_1   # Regime 1 (Turbulent)
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)


Reg_param_theta_sig <- matrix(
  c(
    1, 2.0, 0.03, 0.15, -0.1,  0.05, 2.0, 0.03, 0.15, -0.1
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)



syn_series <- simulate_heston(S0, v0, Reg_chain, Reg_param_1, T, N, M=1, method = "M")

matplot(syn_series$S_paths[1,], type = "l", 
        main = "Heston Model Asset Price Paths", 
        xlab = "Time Step", ylab = "Asset Price S_t")

S_simulated <- syn_series$S_paths[1,]
get_reg(Reg_chain, S_simulated)







library(fHMM)
dax <- download_data(symbol = "^GDAXI")

controls <- set_controls( 
  states      = 3,
  sdds        = "t",
  file        = dax,
  date_column = "Date",
  data_column = "Close",
  logreturns  = TRUE,
  from        = "2000-01-01",
  to          = "2022-12-31"
)
?set_controls
data <- prepare_data(controls)
summary(data)
model <- fit_model(data)
model <- decode_states(model)
model <- compute_residuals(model)
summary(model)

events <- fHMM_events(
  list(dates = c("2001-09-11", "2008-09-15", "2020-01-27"),
       labels = c("9/11 terrorist attack", "Bankruptcy Lehman Brothers", "First COVID-19 case Germany"))
)
plot(model, plot_type = c("sdds","ts"), events = events)
plot(model, plot_type = "pr")
plot(as.matrix(dax[6]),type = "s")




stats::dnorm(
  x = 1,
  mean = 0,
  sd = 1)



