install.packages("devtools")
install.packages("tidyverse")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggmap)

devtools::install_github("RamiKrispin/coronavirus", force = TRUE)
library(coronavirus)
data("coronavirus")

# total cases by country and date
summary_df <- coronavirus %>%
  select(country = Country.Region, date, type, cases) %>%
  group_by(country, date, type) %>%
  summarise(total_cases = sum(cases)) %>%
  pivot_wider(names_from = type,
              values_from = total_cases)

# cumulative summary of total cases by by country, state, and date
summary_df <- summary_df %>%
  group_by(country) %>%
  mutate(confirmed = cumsum(confirmed)) %>%
  mutate(death = cumsum(death)) %>%
  mutate(recovered = cumsum(recovered)) %>%
  filter(country == "US" | country == "China" | country == "Italy")

# Viz 1 - Flatten the curve graph

# step 1 -
# china, italy, us - aligned by first diagnosis date (in lieu of date)
aligned_df <- summary_df %>%
  group_by(country) %>%
  mutate(t = 0 + cumsum(confirmed > 450)) %>%
  filter(t >= 1)


# step 2 -
# plot confirmed cases over time (t0 - tn)
ggplot(data=aligned_df, aes(x=t, y=confirmed, group=country)) +
  geom_line(aes(color=country))+
  ggtitle("COVID-19 \n Rate of Confirmed Cases")+
  xlab("t (Days Since 450 Confirmed)")+
  ylab("Confirmed Cases")+
  labs(caption="t0 corresponds to the first date having 450 total confirmed cases: \n China: 01/22/2020, Italy: 02/26/2020, US: 03/08/2020")

# step 3 -
# model epidemic using SEIR (susceptible, exposed, infectious, recovered) model

install.packages("deSolve")
library(deSolve)              # differential equation solver

# https://www.idmod.org/docs/hiv/model-seir.html

seir_model = function(current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  E = state_values [2]        # exposed
  I = state_values [3]        # infectious
  R = state_values [4]        # recovered

  with (
    as.list (parameters),     # variable names within parameters can be used
         {
           # compute derivatives
           dS = (-beta * S * I)
           dE = (beta * S * I) - (delta * E)
           dI = (delta * E) - (gamma * I)
           dR = (gamma * I)

           # combine results
           results = c(dS, dE, dI, dR)
           list (results)
         }
    )
}

# Parameters
contact_rate = 5                      # number of contacts per day
transmission_probability = 0.06       # transmission probability
infectious_period = 5                 # infectious period
latent_period = 2                     # latent period

# Compute values of beta (tranmission rate) and gamma (recovery rate).
beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period

# Compute Ro - Reproductive number.
Ro = beta_value / gamma_value

# Disease dynamics parameters
parameter_list = c(beta = beta_value, gamma = gamma_value, delta = delta_value)

# Initial values for sub-populations.
US_POP = 327200000
W = US_POP - 2727   # susceptible hosts
X = 1000            # infectious hosts
Y = 727             # recovered hosts
Z = 1000            # exposed hosts

# Compute total population.
N = W + X + Y + Z

# Initial state values for the differential equations.
initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)

# Output timepoints.
timepoints = seq(0, 40, by=1)

# Simulate the COVID-19 epidemic.
# The R function `lsoda` from deSolve provides an interface to the FORTRAN ...
# ... ODE (first-order ordinary differential equations) solver of the same name
output = lsoda(y=initial_values, times=timepoints, func=seir_model, parms=parameter_list)

susceptible_pop <- round(output[, 2] * US_POP, 0)
exposed_pop <- round(output[, 3] * US_POP, 0)
infected_pop <- round(output[, 4] * US_POP, 0)
recovered_pop <- round(output[, 5] * US_POP, 0)

total_cases <- exposed_pop + infected_pop + recovered_pop
total_cases

# US Forecast
country <- rep("US", 40)
confirmed <- total_cases[2:length(total_cases)]
date <- seq(as.Date("2020-03-15"), length.out=40, by="days")

forecast <- data.frame(country, date, confirmed)

aligned_df <- bind_rows(aligned_df, forecast)
aligned_df <- aligned_df %>%
  group_by(country) %>%
  mutate(t = 0 + cumsum(confirmed > 450)) %>%
  filter(t >= 1)

ggplot(data=aligned_df, aes(x=t, y=confirmed, group=country)) +
  geom_line(aes(color=country))+
  ggtitle("COVID-19 \n Rate of Confirmed Cases")+
  xlab("t (Days Since 450 Confirmed)")+
  ylab("Confirmed Cases")+
  labs(caption="t0 corresponds to the first date having 450 total confirmed cases: \n China: 01/22/2020, Italy: 02/26/2020, US: 03/08/2020")

# # susceptible hosts over time
# plot (S ~ time, data = output, type='b', ylim = c(0,1), col = 'blue', ylab = 'S, E, I, R', main = 'COVID-19 Epidemic')
# par (new = TRUE)      # remain on same frame

# # exposed hosts over time
# plot (E ~ time, data = output, type='b', ylim = c(0,1), col = 'pink', ylab = '', axes = FALSE)
# par (new = TRUE)      # remain on same frame

# # infectious hosts over time
# plot (I ~ time, data = output, type='b', ylim = c(0,1), col = 'red', ylab = '', axes = FALSE)
# par (new = TRUE)      # remain on same frame

# # recovered hosts over time
# plot (R ~ time, data = output, type='b', ylim = c(0,1), col = 'green', ylab = '', axes = FALSE)


# step 4 - work backwards using FME
install.packages("FME")
library(FME)

SEIR <- function (pars, S_0 = 49889, E_0 = 100, I_0 = 10, R_0 = 1) {

  derivs <- function (time, y, pars) {
  with (as.list(c(pars, y)), {
  # Compute values of beta (tranmission rate) and gamma (recovery rate).
  beta <- contact_rate * transmission_probability
  gamma <- 1 / infectious_period
  delta <- 1 / latent_period

  # Compute Ro - Reproductive number.
  Ro <- beta / gamma

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

# Disease dynamics parameters
pars = c(contact_rate = 10,
  transmission_probability = 0.10,
  infectious_period = 5,
  latent_period = 2)

out <- SEIR(pars = pars)

# Plot the modeled rate of infected to exposed population
par(mfrow = c(1, 2))
plot(out$time, out$infected,
  main = "Infected Population", ylab = "infected",xlab = "time", type = "b")
plot(out$time, out$E, main = "Exposed Poopulation", ylab = "-", xlab = "time", type = "b")
par(mfrow = c(1, 1))

# Simulated data for synthetic experiment
DataInfected <- cbind(
  time = out$time,
  infected = round(out$infected + rnorm(sd = 4.5, n = length(out$infected)), 0),
  sd = 4.5)

DataE <- cbind(
  time = out$time,
  E = out$E + rnorm(sd = 0.01, n = length(out$E)),
  sd = 0.01)

# Weighted residuals of the model output versus the data and calculated
# sum of squared residuals
SEIRcost <- function(pars) {
  out <- SEIR(pars = pars)
  cost <- modCost(model = out, obs = DataInfected, err = "sd")
  return(modCost(model = out, obs = DataE, err = "sd", cost = cost))
  }

# Residuals of model and pseudodata
SEIRcost(pars)$model
plot(SEIRcost(pars), xlab="time")

# Sensitivity functions of model output to parameters
Sfun <- sensFun(SEIRcost, pars)
summary(Sfun)
plot(Sfun, which = c("infected", "E"), xlab="time", lwd = 2)

# Pairwise plot of sensitivity functions
pairs(Sfun, which = c("infected", "E"), col = c("blue", "green"))

# Multivariate parameter identifiability
# Collinearity for all parameter combinations (wrt 2 parameters: infected and E)
ident <- collin(Sfun)
head(ident, n = 20)

# Collinearity plot
plot(ident, log = "y")

# Collinearity for all parameter combinations (wrt 1 parameter: infected or E)
collin(Sfun, parset = c("contact_rate", "infectious_period"), which = "infected")
collin(Sfun, parset = c("contact_rate", "infectious_period"), which = "E")

# Fitting the model to data

# The log transformation:
# (1) ensures that the parameters remain positive during the fitting, and
# (2) deals with the fact that the parameter values are spread over six orders of magnitude
# SEIRcost2 <- function(lpars) {
#   SEIRcost(c(exp(lpars)))
# }

# Model is fitted to the data, and best-fit parameters and residual sum of squares shown.
# Fit <- modFit(f = SEIRcost2, p = log(pars))
Fit <- modFit(f = SEIRcost, p = pars)
deviance(Fit)
summary(Fit)

# For comparison, the initial model and the best-fit model are plotted against the data
ini <- SEIR(pars = c(pars))
final <- SEIR(pars = c(coef(Fit)))

par(mfrow = c(1,2))
plot(DataInfected, xlab = "time", ylab = "infected")
lines(ini$time, ini$infected, lty = 2)
lines(final$time, final$infected)
legend("topright", c("data", "initial", "fitted"), lty = c(NA,2,1), pch = c(1, NA, NA))

plot(DataE, xlab = "time", ylab = "E")
lines(ini$time, ini$E, lty = 2)
lines(final$time, final$E)
par(mfrow = c(1, 1))

# MCMC
# Bayesian methods can be used to derive the data-dependent probability distribution
# of the parameters.

# The algorithm is started with the optimal parameter set (Fit$par)

# The prior error variance var0 is chosen to be the mean of the unweighted squared
# residuals from the model.
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5

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

# TODO: refactor parameters to take contact_rate, number of contacts per day, transmission_probability, infectious_period, infectious period, and latent_period.
# TODO: in lieu of using synthetic data, use COVID-19 from a region in China.
