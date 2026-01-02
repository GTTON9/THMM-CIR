# source("Model_Simulation.R")
# source("Gen_fit.R")

v0 <- 0.03
S0 <- 100
nstates = 2
T = 1
n_days <- 252
n_intraday <- 2400
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(10, 5)
theta <- c(0.03, 0.6)
sigma <- c(0.05, 0.05)
rho <- c(-0.1, -0.1)




v0 <- 0.1
S0 <- 100
nstates = 2
T = 1
n_days <- 252
n_intraday <- 2400
mu <- c(0.3, 0.3)
kappa <- c(10, 5)
theta <- c(0.1, 0.5)
sigma <- c(0.1, 0.1)
rho <- c(-0.1, -0.1)


N <- 252* n_intraday
Reg_chain_year <- simulate_Reg(series_length = N/n_intraday * T, Reg_tran = Gamma)
Reg_chain <- rep(Reg_chain_year, each = n_intraday)
Reg_param <- cbind(mu, kappa, theta, sigma, rho)
sim_series_1_year  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E",interp = F, seed = 999)
S_simulated <- sim_series_1_year$S_paths
plot(S_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)], type = "l")
V_simulated <- sim_series_1_year$V_paths


par(mfrow = c(1, 2))

plot(
  S_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)],
  type = "l",
  main = "S path",
  col = "blue")

plot(
  V_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)],
  type = "l",
  main = "V path",
  col = "red"
)

par(mfrow = c(1, 1))




S <- S_simulated


RV_V <- numeric(n_days)

for (t in 1:n_days) {
  idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
  S_day <- S_simulated[idx]
  r_day <- diff(log(S_day))      # length 3
  RV_V[t] <- sum(r_day^2) * 252 
}


plot(RV_V, col = "grey")
lines(RV_V, col = "blue")


# RV_V <- lowess(RV_V, f = 0.03)$y




series_length <- length(RV_V)
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)







# source("Gen_fit.R")
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












library(ggplot2)





true_vol <- V_simulated[seq(1, length(V_simulated), by = n_intraday)]

result <- plot_cir_confidence_simple(
  true_vol = V_daily,
  param = param,
  states_estimate = states_estimate,
  dt = 1/252,
  interval_step = 10
)

print(result$plot)



print(summary(result$ci_data))


param



S_daily <- S_simulated[seq(n_intraday, length(S_simulated), by = n_intraday)]
V_daily <- V_simulated[seq(n_intraday, length(V_simulated), by = n_intraday)]


result_mle <- estimate_mu_simple_mle(
  S_daily = S_daily,
  S_intraday = S_simulated,
  n_intraday = n_intraday,
  states_estimate = states_estimate,
  nstates = 2,
  dt = 1/252
)

result_mle$mu






