install.packages("FME")
library(FME)


SEIR <- function (pars, S_0 = 49889, E_0 = 100, I_0 = 10, R_0 = 1) {

  derivs <- function (time, y, pars) {
  with (as.list(c(pars, y)), {
  dS <- (-beta * S * I)
  dE <- (beta * S * I) - (delta * E)
  dI <- (delta * E) - (gamma * I)
  dR <- (gamma * I)

  return(list(c(dS, dE, dI, dR), infected = I * N))
  })
  }

  # calc. total population
  N <- S_0 + E_0 + I_0 + R_0

  # initial conditions
  y <- c(S = S_0/N, E = E_0/N, I = I_0/N, R = R_0/N)

  times <- seq(0, 40, by=1)
  out <- ode(y = y, parms = pars, times = times, func = derivs)

  as.data.frame(out)
  }

# Parameters
contact_rate = 10                     # number of contacts per day
transmission_probability = 0.10       # transmission probability
infectious_period = 5                 # infectious period
latent_period = 2                     # latent period

# Compute values of beta (tranmission rate) and gamma (recovery rate).
beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period

# Compute Ro - Reproductive number.
Ro = beta_value / gamma_value

# Disease dynamics parameters
pars = c(beta = beta_value, gamma = gamma_value, delta = delta_value)

out <- SEIR(pars = pars)

# Plot the modeled rate of infected to exposed population
par(mfrow = c(1, 2))
plot(out$time, out$infected,
  main = "Infected Population", ylab = "infected",xlab = "time", type = "b")
plot(out$time, out$E, main = "Exposed Poopulation", ylab = "-", xlab = "time", type = "b")
par(mfrow = c(1, 1))

# The FME algorithms will be tested on simulated data. Such synthetic experiments
# are often used to study parameter identifiability or to test fitting routines.

# Simulated data for synthetic experiment
DataInfected <- cbind(
  time = out$time,
  infected = round(out$infected + rnorm(sd = 7.0, n = length(out$infected)), 0),
  sd = 7.0)

DataE <- cbind(
  time = out$time,
  E = out$E + rnorm(sd = 0.02, n = length(out$E)),
  sd = 0.02)

# Weighted residuals of the model output versus the data and calculated
# sum of squared residuals
SEIRcost <- function(pars) {
  out <- SEIR(pars = pars)
  cost <- modCost(model = out, obs = DataInfected, err = "sd")
  return(cost)
#   return(modCost(model = out, obs = DataE, err = "sd", cost = cost))
  }

# Residuals of model and pseudodata
SEIRcost(pars)$model
plot(SEIRcost(pars), xlab="time")

# Sensitivity functions of model output to parameters
Sfun <- sensFun(SEIRcost, pars)
summary(Sfun)
plot(Sfun, which = c("infected"), xlab="time", lwd = 2)

# Pairwise plot of sensitivity functions
pairs(Sfun, which = c("infected"), col = c("blue"))

# Multivariate parameter identifiability
# Collinearity for all parameter combinations (wrt 2 parameters: infected and E)
ident <- collin(Sfun)
head(ident, n = 20)

# Collinearity plot
plot(ident, log = "y")

# Collinearity for all parameter combinations (wrt 1 parameter: infected or E)
collin(Sfun, parset = c("beta", "gamma"), which = "infected")
# collin(Sfun, parset = c("beta", "gamma"), which = "E")

# Fitting the model to data

# The log transformation:
# (1) ensures that the parameters remain positive during the fitting, and
# (2) deals with the fact that the parameter values are spread over six orders of magnitude
SEIRcost2 <- function(lpars) {
  SEIRcost(c(exp(lpars)))
}

# Model is fitted to the data, and best-fit parameters and residual sum of squares shown.
Fit <- modFit(f = SEIRcost2, p = log(pars))
exp(coef(Fit))
deviance(Fit)
summary(Fit)

# For comparison, the initial model and the best-fit model are plotted against the data
ini <- SEIR(pars = c(pars))
final <- SEIR(pars = c(exp(coef(Fit))))

par(mfrow = c(1,1))
plot(DataInfected, xlab = "time", ylab = "infected")
lines(ini$time, ini$infected, lty = 2)
lines(final$time, final$infected)
legend("topright", c("data", "initial", "fitted"), lty = c(NA,2,1), pch = c(1, NA, NA))

# plot(DataE, xlab = "time", ylab = "E")
# lines(ini$time, ini$E, lty = 2)
# lines(final$time, final$E)
# par(mfrow = c(1, 1))

# MCMC
# Bayesian methods can be used to derive the data-dependent probability distribution
# of the parameters.

# The algorithm is started with the optimal parameter set (Fit$par)

# The prior error variance var0 is chosen to be the mean of the unweighted squared
# residuals from the model.
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled

# Fit MCMC
MCMC <- modMCMC(f = SEIRcost2, p = Fit$par, niter = 5000, jump = cov0,
  var0 = var0, wvar0 = 0.1, updatecov = 50)

MCMC$pars <- exp(MCMC$pars)
summary(MCMC)

# Results of the MCMC application.
plot(MCMC, Full = TRUE)

# Pairs plot of the MCMC application.
pairs(MCMC, nsample = 1000)

# The effect of the parameter uncertainty on the model output
sR <- sensRange(func = SEIR, parms = pars, parInput = MCMC$par)

# Sensitivity range based on parameter distribution as generated with the MCMC application.
plot(summary(sR), xlab = "time")

parRange <- cbind(min = 0.75 * pars, max = 1.25 * pars)
crlfun <- function (pars) return(meanS = mean(SEIR(pars)$S))
CRL <- modCRL(fun = crlfun, parRange = parRange, num = 500)
cor(CRL)

plot(CRL, ylab = "susceptibles", trace = TRUE)
