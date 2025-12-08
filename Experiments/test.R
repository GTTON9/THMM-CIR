source("Model_Simulation.R")

N <- 252
v0 <- 0.03
S0 <- 100
n_intraday <- 2400
nstates = 2
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(10, 5)
theta <- c(0.03, 0.6)
sigma <- c(0.2, 0.2)
rho <- c(-0.1, -0.1)


N <- 252 * n_intraday
Reg_chain_year <- simulate_Reg(series_length = N/n_intraday, Reg_tran = Gamma)
Reg_chain <- rep(Reg_chain_year, each = n_intraday)

Reg_param <- cbind(mu, kappa, theta, sigma, rho)
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E",interp = F, seed = 999)
S_simulated <- sim_series$S_paths

V_simulated <- sim_series$V_paths


par(mfrow = c(1, 2))

  plot(
    S_simulated,
    type = "l",
    main = "S path",
    col = "blue")
  
  plot(
    V_simulated,
    type = "l",
    main = "V path",
    col = "red"
  )

par(mfrow = c(1, 1))




S <- S_simulated
n_days <- 252
n_intraday <- 2400

RV_V <- numeric(n_days)

for (t in 1:n_days) {
  idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
  S_day <- S[idx]
  r_day <- diff(log(S_day))      # length 3
  RV_V[t] <- sum(r_day^2) * 250 
}


plot(V_simulated[seq(1, length(V_simulated), by = n_intraday)], type = 'l', col = "black")
lines(RV_V, col = "blue")
lines(lowess(RV_V, f = 0.02), col = "red")



RV_V <- lowess(RV_V, f = 0.02)$y




series_length <- length(RV_V)
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)

source("Gen_fit.R")
my_data_df <- data.frame(
  Date = date_sequence,
  Var = RV_V
)
colnames(my_data_df) <- c("Date", "Var")



# followed the example
series_control <- Heston_set_controls( 
  states      = 2,     # 2 state
  sdds        = "Heston",         
  date_column = "Date",
  file        = my_data_df, 
  data_column = "Var",      
  logreturns  = FALSE,         
  from        = date_sequence[1],             
  to          = date_sequence[length(date_sequence)],
  runs = 10
)


data_hmm <- prepare_data(series_control)
model_hmm <- Heston_fit_model(data_hmm) 

final_model <- decode_states_heston(model_hmm) 
states_estimate <- final_model$decoding
states_estimate

plot(states_estimate, col = "blue")
lines(Reg_chain_year+1)
param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)


p <- plot_viterbi(RV_V, nstates, param$Gamma, 
                  param$kappa, param$theta, 
                  param$sigma, Reg_chain_year)
plot(p)
param



























source("Model_Simulation.R")

source("Gen_fit.R")
N <- 252
v0 <- 1
S0 <- 100
nstates = 2
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(2, 1)
theta <- c(0.1, 0.2)
sigma <- c(0.05, 0.05)
rho <- c(-0.1, -0.1)

N <- 252 * n_intraday
Reg_chain_year <- simulate_Reg(series_length = N/n_intraday, Reg_tran = Gamma)
Reg_chain <- rep(Reg_chain_year, each = n_intraday)

Reg_param <- cbind(mu, kappa, theta, sigma, rho)
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E",interp = F, seed = 999)
S_simulated <- sim_series$S_paths

V_simulated <- sim_series$V_paths


par(mfrow = c(1, 2))

plot(
  S_simulated,
  type = "l",
  main = "S path",
  col = "blue")

plot(
  V_simulated,
  type = "l",
  main = "V path",
  col = "red"
)

par(mfrow = c(1, 1))




S <- S_simulated
n_days <- 252
n_intraday <- 2400

RV_V <- numeric(n_days)

for (t in 1:n_days) {
  idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
  S_day <- S[idx]
  r_day <- diff(log(S_day))      # length 3
  RV_V[t] <- sum(r_day^2) * 250 
}


plot(V_simulated[seq(1, length(V_simulated), by = n_intraday)], type = 'l', col = "black")
lines(RV_V, col = "blue")
lines(lowess(RV_V, f = 0.05), col = "red")



RV_V <- lowess(RV_V, f = 0.05)$y




series_length <- length(RV_V)
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)

source("Gen_fit.R")
my_data_df <- data.frame(
  Date = date_sequence,
  Var = RV_V
)
colnames(my_data_df) <- c("Date", "Var")



# followed the example
series_control <- Heston_set_controls( 
  states      = 2,     # 2 state
  sdds        = "Heston",         
  date_column = "Date",
  file        = my_data_df, 
  data_column = "Var",      
  logreturns  = FALSE,         
  from        = date_sequence[1],             
  to          = date_sequence[length(date_sequence)],
  runs = 10
)


data_hmm <- prepare_data(series_control)
model_hmm <- Heston_fit_model(data_hmm) 

final_model <- decode_states_heston(model_hmm) 
states_estimate <- final_model$decoding
states_estimate

plot(states_estimate, col = "blue")
lines(Reg_chain_year+1)
param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)


p <- plot_viterbi(RV_V, nstates, param$Gamma, 
                  param$kappa, param$theta, 
                  param$sigma, Reg_chain_year)
plot(p)
param
