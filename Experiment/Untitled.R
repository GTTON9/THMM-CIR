library(gsl) # modified Bessel Function
library(fHMM)
install.packages("fHMM")
# install.packages("gsl")

kappa <- 2.0  
theta <- 0.03
sigma <- 0.15

k <- 1.0    
V_t <- 0.03   
V_t_plus_k <- 0.04
V_t_plus_k <- seq(0.001,0.07, 0.001)

C_numerator <- 2 * kappa
C_denominator <- (1 - exp(-kappa * k)) * (sigma^2)
C <- C_numerator / C_denominator

q_numerator <- 2 * kappa * theta
q_denominator <- sigma^2
q <- (q_numerator / q_denominator) - 1

u <- C * V_t * exp(-kappa * k)

v <- C * V_t_plus_k


argument <- 2 * sqrt(u * v)
I_q <- besselI(argument, nu = q)


density_factor <- (v / u)^(q / 2)
transition_density <- C * exp(-u - v) * density_factor * I_q
transition_density
plot(transition_density, type = "l")


G_alpha <- (2 * theta * kappa)/(sigma^2)
G_beta <- (2 * kappa)/(sigma^2)


pdf_values <- dgamma(
  x = seq(0.00, 0.1, length.out = 500),
  shape = G_alpha,
  rate = G_beta
)
plot(pdf_values, type = "l")












controls <- set_controls(
  states  = 2,
  sdds    = "gamma",
  horizon = 1000
)
par <- fHMM_parameters(
  controls = controls,
  Gamma    = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2), 
  mu       = c(1, 3), 
  sigma    = c(1, 3)
)
sim <- simulate_hmm(
  controls        = controls,
  true_parameters = par
)
plot(sim$data, col = sim$markov_chain, type = "b")

(parUncon <- par2parUncon(par, controls))


ll_hmm(parUncon, sim$data, controls)
ll_hmm(parUncon, sim$data, controls, negative = TRUE)
remove.packages("fHMM")
 install.packages("fHMM")
