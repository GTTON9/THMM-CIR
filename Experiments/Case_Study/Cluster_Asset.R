
source("Model_Simulation.R")
T = 1
N = 25000
v0 = 0.03
S0 =100
states <- 2 

Reg_chain <- simulate_Reg(series_length = N)
Reg_chain_year <- simulate_Reg(series_length = N/100)
Reg_chain <- rep(Reg_chain_year, each = 100)

Reg_param <- matrix(
  
  # feller condition: 2 * kappa * theta > sigma^2
  c( # mu, kappa, theta, sigma, rho
    0.5,   10,    0.03,   0.05,  -0.1, # calm
    0.5,   5,     0.6,   0.05,  -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE 
)

Reg_param <- matrix(

  # feller condition: 2 * skappa * theta > sigma^2
  c( # mu, kappa, theta, sigma, rho
    0.5,   2,    0.1,   0.05,  -0.1, # calm
    0.5,   1,     0.2,   0.05,  -0.1 # turbulent
  ),
  nrow = 2,
  ncol = 5,
  byrow = TRUE
)


v0 = 1





# 62-144 regime change
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T = 1, N, M=1, method = "E")
S_simulated <- sim_series$S_paths
plot(S_simulated,  type = "l")

V_simulated <- sim_series$V_paths
plot(V_simulated, type = "l")


S <- S_simulated
n_days <- 250
n_intraday <- 100

RV_V <- numeric(n_days)

for (t in 1:n_days) {
  idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
  S_day <- S[idx]
  r_day <- diff(log(S_day))      # length 3
  RV_V[t] <- sum(r_day^2) * 250 
}



plot(V_simulated[seq(1, length(V_simulated), by = 100)], type = 'l', ylim = c(0,1), col = "black")
lines(RV_V, col = "blue")
lines(lowess(RV_V, f = 0.1), col = "red")
legend("topright", 
       legend = c("True Variance (V_t)", 
                  "Realized Variance Estimate", 
                  "Smoothed RV (lowess, f=0.1)"),
       col = c("black", "blue", "red"),
       lwd = c(2, 1.5, 2),
       lty = c(1, 1, 2),
       bty = "n",
       cex = 0.5)



RV_V <- lowess(RV_V, f = 0.1)$y






set.seed(999)
series_length <- length(RV_V)
start_date <- as.Date("2024-01-01")
date_sequence <- seq(from = start_date, by = "day", length.out = series_length)


my_data_df <- data.frame(
  Date = date_sequence,
  Var = RV_V
)
colnames(my_data_df) <- c("Date", "Var")


source("Gen_fit.R")

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

param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
param




p <- plot_viterbi(RV_V, nstates, param$Gamma, 
                  param$kappa, param$theta, 
                  param$sigma, Reg_chain_year[1:249])
plot(p)




Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
N <- 250
v0 <- 1
S0 <- 100
nstates = 2
mu <- c(0.5, 0.5)
kappa <- c(2, 1)
theta <- c(0.1, 0.2)
sigma <- c(0.1, 0.1)
rho <- c(-0.1, -0.1)




states_estimate <- viterbi(RV_V, nstates, Gamma, kappa, theta, sigma)
states_estimate

plot(states_estimate, col = "blue")
lines(Reg_chain_year+1)


Gen_fit_1(Gamma, mu, kappa, theta, sigma, rho,  V_ = TRUE , plot_path = TRUE)

Gen_fit_1(Gamma, mu, kappa, theta, sigma, rho, input_ = "S_",  V_ = TRUE , plot_path = TRUE)

