# devtools::install_github("RamiKrispin/coronavirus", force = TRUE)
# library(coronavirus)
# data("coronavirus")

install.packages("devtools")
install.packages("tidyverse")

library(dplyr)
library(tidyr)

setwd("~/github/coronavirus/")
coronavirus <- read.csv(file = 'covid19.csv')
# head(coronavirus)

# # total cases by country and date
# coronavirus <- coronavirus %>%
#   select(country = Country.Region, date, type, cases) %>%
#   group_by(country, date, type) %>%
#   summarise(total_cases = sum(cases)) %>%
#   pivot_wider(names_from = type,
#               values_from = total_cases)

summary_df <- coronavirus %>%
  select(country = Country.Region, date, confirmed, recovered, deaths) %>%
  mutate(date = mdy(date))  %>%
  group_by(country, date) %>%
  summarise(confirmed = sum(confirmed), recovered = sum(recovered), deaths = sum(deaths))  %>%
  arrange(country, date) %>%
  filter(country == "US" | country == "China" | country == "Italy") %>%
  mutate(infected = (confirmed - recovered - deaths)) %>%
  filter(date < as.Date("2020-03-25"))

ggplot(summary_df, aes(x=date, y=infected)) +
  geom_line(aes(color=country)) +
  theme_linedraw() +
  ggtitle("COVID-19 Infected Individuals") +
  xlab('t') + ylab('Infected') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")

# # cumulative summary of total cases by by country, state, and date
# summary_df <- coronavirus %>%
#   select(country = Country.Region, date, confirmed, recovered, deaths) %>%
#   group_by(country, date) %>%
#   summarise(confirmed = sum(confirmed), recovered = sum(recovered), deaths = sum(deaths)) %>%
#   group_by(country) %>%
#   mutate(confirmed = cumsum(confirmed)) %>%
#   mutate(deaths = cumsum(deaths)) %>%
#   mutate(recovered = cumsum(recovered)) %>%
#   mutate(infected = (confirmed - deaths)) %>%
#   filter(country == "US" | country == "China" | country == "Italy")

china_df <- summary_df %>%
  group_by(country) %>%
  mutate(time = -1 + cumsum(confirmed >= 0)) %>%
  filter(country == "China")

# Simulated data for synthetic experiment
ChinaInfected <- cbind(
  time = china_df$time,
  infected = round(china_df$infected),
  sd = sd(china_df$infected))


install.packages("FME")
library(FME)

SEIR <- function (pars, S_0 = 120000, E_0 = 3000, I_0 = 500, R_0 = 0) {

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

  times <- seq(0, 62, by=1)
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

ggplot(as.data.frame(ChinaInfected), aes(x=time, y=infected)) +
  geom_point(shape=1, color="black", size=2) +
  theme_linedraw() +
  ggtitle("COVID-19 Infected Individuals in China") +
  xlab('t') + ylab('Infected (China)') +
  theme(plot.title = element_text(hjust = 0.5))

# # Plot the modeled rate of infected to exposed population
# par(mfrow = c(1, 2))
# plot(out$time, out$infected,
#   main = "Infected Population", ylab = "infected",xlab = "time", type = "b")
# plot(out$time, out$E, main = "Exposed Poopulation", ylab = "-", xlab = "time", type = "b")
# par(mfrow = c(1, 1))


# Weighted residuals of the model output versus the data and calculated
# sum of squared residuals
SEIRcost <- function(pars) {
  out <- SEIR(pars = pars)
  cost <- modCost(model = out, obs = ChinaInfected, err = "sd")
  return(cost)
#   return(modCost(model = out, obs = DataE, err = "sd", cost = cost))
  }

# Residuals of model and pseudodata
SEIRcost(pars)$model
# plot(SEIRcost(pars), xlab="time")

# Sensitivity functions of model output to parameters
Sfun <- sensFun(SEIRcost, pars)
summary(Sfun)

# plot(Sfun, which = c("infected"), xlab="time", lwd = 2)
ggplot(Sfun, aes(x=x)) +
  geom_line(aes(y = beta, color = "#00BA38")) +
  geom_line(aes(y = gamma, color= "#F8766D")) +
  geom_line(aes(y = delta, color= "#619CFF")) +
  scale_color_discrete(name = "Parameters", labels = c("Beta", "Delta", "Gamma")) +
  ggtitle("Parameter Sensitivity wrt Infected") +
  xlab('t') + ylab('Sensitivity') +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")

# Pairwise plot of sensitivity functions
# pairs(Sfun, which = c("infected"), col = c("blue"))
ggpairs(Sfun, columns = c("beta", "gamma", "delta"), title = "Parameter Sensitivity Correlation") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")

# Multivariate parameter identifiability
# Collinearity for all parameter combinations (wrt 2 parameters: infected and E)
ident <- collin(Sfun)
head(ident, n = 20)

# Collinearity plot
# plot(ident, log = "y")
ggplot(ident, aes(x=N, y=collinearity)) + geom_point() +  theme_linedraw() +
  ggtitle("Collinearity for all parameter combinations") +
  xlab('N (number of parameters)') + ylab('Collinearity') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")

# Collinearity for all parameter combinations (wrt 1 parameter: infected or E)
# collin(Sfun, parset = c("beta", "gamma"), which = "infected")
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
round(exp(coef(Fit)), 5)
deviance(Fit)
# summary(Fit)

exp(coef(Fit))["beta"] / exp(coef(Fit))["gamma"]

# For comparison, the initial model and the best-fit model are plotted against the data
ini <- SEIR(pars = c(pars))
final <- SEIR(pars = c(exp(coef(Fit))))

# par(mfrow = c(1,1))
# plot(ChinaInfected, xlab = "time", ylab = "infected")
# lines(ini$time, ini$infected, lty = 2)
# lines(final$time, final$infected)
# legend("topleft", c("data", "initial", "fitted"), lty = c(NA,2,1), pch = c(1, NA, NA))

ggplot(NULL, aes(x=time, y=infected)) +
  geom_point(data = china_df, shape=1, color="black", size=2) +
  geom_line(data=ini, aes(color = "#619CFF")) +
  geom_line(data=final, aes(color= "#F8766D")) +
  theme_linedraw() +
  ggtitle("SEIR Model: Initial vs. Fitted") +
  xlab('t') + ylab('Infected (China)') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  scale_color_discrete(name = "Legend", labels = c("Initial", "Fitted"))

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
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5

# Fit MCMC
MCMC <- modMCMC(f = SEIRcost2, p = Fit$par, niter = 5000, jump = cov0,
  var0 = var0, wvar0 = 0.1, updatecov = 50)

MCMC$pars <- exp(MCMC$pars)
summary(MCMC)

# Results of the MCMC application.
# plot(MCMC, Full = TRUE)

# Pairs plot of the MCMC application.
# pairs(MCMC, nsample = 1000)
ggpairs(as.data.frame(MCMC$pars),
  columns = c("beta", "gamma", "delta"), title = "MCMC Parameter Estimates") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5))

# The effect of the parameter uncertainty on the model output
sR <- sensRange(func = SEIR, parms = pars, parInput = MCMC$par)

# Sensitivity range based on parameter distribution as generated with the MCMC application.
# plot(summary(sR), xlab = "time")

sR_summary <- cbind(Category = rownames(summary(sR)), summary(sR))
sR_summary$Category <- gsub('[[:digit:]]+', '', sR_summary$Category)
sR_summary$Category <- factor(sR_summary$Category, levels = c("S", "E", "I", "R", "infected"))

ggplot(sR_summary[(sR_summary$Category != "infected"), ], aes(x=x, color=Category)) +
  geom_line(aes(y = q05), linetype="longdash") +
  geom_line(aes(y = q95), linetype="longdash") +
  geom_line(aes(y = Mean)) +
  theme_linedraw() +
  ggtitle("China SEIR Fitted Model") +
  xlab('t') + ylab('% Susceptibles') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")

# parRange <- cbind(min = 0.75 * pars, max = 1.25 * pars)
# crlfun <- function (pars) return(meanS = mean(SEIR(pars)$S))
# CRL <- modCRL(fun = crlfun, parRange = parRange, num = 500)
# cor(CRL)

# plot(CRL, ylab = "susceptibles", trace = TRUE)

beta = 0.50510
gamma = 0.07098
delta = 0.28953

# calculate R0
R0 = beta / gamma

# back into SEIR model assumptions from optimal model params discovered
# beta
transmission_proba <- seq(0.05, 1, .05)
contact_rate <- c()
for (x in  transmission_proba) {
contact_rate <- c(contact_rate, 0.50510 / x)}
infectious_force <- cbind(transmission_proba, contact_rate)

ggplot(as.data.frame(infectious_force), aes(x=transmission_proba, y=contact_rate)) +
  geom_point(shape=1, color="black", size=2) +
  theme_linedraw() +
  ggtitle("Infectious Force") +
  xlab('Transmission Probability') + ylab('Contact Rate') +
  theme(plot.title = element_text(hjust = 0.5))

# gamma
infectious_period <- 1 / gamma

# delta
latent_period <- 1 / delta