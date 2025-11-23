install.packages("oeli")
library(oeli)
install.packages("Rpcc")
library(Rpcc)
install.packages(c("Rcpp", "RcppArmadillo", "checkmate", "pracma", "progress"))


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

sourceCpp("src/ll.cpp")
sourceCpp("src/RcppExports.cpp")
install.packages("gfortran")
library(gfortran)

devtools::load_all()

ll_hmm(parUncon, sim$data, controls)
ll_hmm(parUncon, sim$data, controls, negative = TRUE)

optimization <- nlm(
  f = ll_hmm, p = parUncon, observations = sim$data, controls = controls, negative = TRUE
)





install.packages("fHMM")
library("fHMM")


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

optimization <- nlm(
  f = ll_hmm, p = parUncon, observations = sim$data, controls = controls, negative = TRUE
)

(estimate <- optimization$estimate)

class(estimate) <- "parUncon"
estimate <- parUncon2par(estimate, controls)

par$Gamma
estimate$Gamma

par$mu
estimate$mu

par$sigma
estimate$sigma





print("a")





