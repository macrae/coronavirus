On Estimating the Reproductive Rate of the COVID-19 Virus

By Sean MacRae
Date: 03/30/2020

# Introduction

I've recently come across many compelling data visualizations and statistical models about the COVID-19 virus. The ones that really get my attention demonstrate the effectiveness of social distancing <insert link>, like this nifty little simulation by the Washington Post <insert link>. Many of the COVID-19 analyses explore epidemiological models, like this one <insert link> that uses an SEIR model to estimate the effect of social distancing.

I love learning about new models that describe some natural phenomenon. I have experienced many times when an approach to solving a problem in one field translates well into the problem-space in another domain. I've seen actuarial models used to forecast conversion funnels of paid marketing campaigns, and survival analysis used to estimate the price elasticity of student loan interest rates.

I wanted to learn more, so I decided that while I'm quarantined with my wife and two kids watching Frozen (no less than two times a day), I would learn by doing. I set out to model the John Hopkins University data <insert link> with an SEIR model because since COVID-19 hit the global scene, SEIR has been getting a bit of press.

# The SEIR Model

SEIR stands for Susceptible, Exposed, Infected, Recovered. The model is described by the rates at which a virus causes susceptible people within a population to be exposed, infected, and recover from the contagion.

<insert equations here>

... describe equations ....

Where beta is the _transmission rate_ of the virus, gamma is the _recovery rate_ of an individual, and delta is the _incubation rate_ or the rate at which latent individuals become infectious. These model parameters are defined as:

```
beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period
```

An interesting SEIR metric is `R0` the _reproductive rate_, which is calculated by:

```
R0 = beta_value / gamma_value
```

`R0` is the rate of transmission to recovery, which proxies the number of people that a single individual infects during the epidemic. If `R0` is 2, for example, then on average, each person infects 2 other people before recovering. The reproductive rate is an interesting metric because it's a measure of virality - the higher the `R0`, the more viral the epidemic.

Since the SEIR model is defined using the _rate of change_ between these four conditions, the model is entirely dependent on the initial conditions of the population (e.g., susceptibles, exposed, infected, and recovered). Some initial conditions may be entirely unknown, for example, the data we have reports "confirmed cases" (synonymous with "infected") but it does not report on the "susceptible" or "exposed" population because we don't (or can't) measure that. We can only try and make some informed assumptions. We must also make assumptions about model parameters beta, gamma, and delta, or the following values:

```
contact_rate = 10                     # number of contacts per day
transmission_probability = 0.10       # transmission probability
infectious_period = 5                 # infectious period
latent_period = 2                     # latent period
```

# Modeling COVID-19

To model COVID-19 with an SEIR model, we need to make some assumptions about the initial conditions of the population as well as the model parameters used to calculate beta, gamma, and delta. One way to do this is through research and analysis, which is better than using our own judgment and intuition, but what if some of these values are actually unknowable? We may have some research that supports a strong belief about the infectious period of the virus because that measurement is better suited to a controlled scientific experiment. Still, something like the transmission probability may be a factor of many social or environmental causes, and quite likely, can vary from region to region. Similarly, the contact rate between individuals of a population isn't something we can measure at all. In essence, it's only something we can have a belief, or make an assumption, about.

These unknowable, yet key, assumptions is where the inspiration for my own COVID-19 analysis came from. I didn't want to assume anything, nor did I want to experiment with arbitrarily controlling for and adjusting this sizable set of assumptions in an attempt to model the data. What if instead of "working forward" by setting assumptions and evaluating the output against the data, we did the inverse? What if we tried "working backward" by taking some data and finding the optimal set of assumptions that describe it?

The hypothesis I wish to test is that the number of infected individuals can be described using an SEIR model. As the virus spreads, at any given time, different countries are in different parts of the infection curve. If my goal is to inversely model the parameters given the data, then I need a population that has a more seasoned infection curve.

# TODO: normalize by 1,000 people (per capita)
<insert confirmed cases by country graph>
<insert infected by country graph>

Since patient zero was in Wuhan, China, it makes sense for the most complete curve to be from China, where the spread of the virus has plateaued. That being the case, I will rephrase the hypothesis to the number of infected individuals _in China_ can be described using an SEIR model.

# Model Research

Having defined the goal of this research through the formulation of a hypothesis, I started researching how people have solved similar problems. I definitely do not want to reinvent the wheel, so I started Googling things like "inverse modeling differential equations." Luckily for me, I stumbled upon a library in R called FME <https://cran.r-project.org/web/packages/FME/vignettes/FME.pdf> that does precisely that! From one of their R vignettes:

```
Mathematical simulation models are commonly applied to analyze experimental or environmental data and eventually to acquire predictive capabilities. Typically these models depend on poorly defined, unmeasurable parameters that need to be given a value. Fitting a model to data, so-called inverse modeling, is often the sole way of finding
reasonable values for these parameters. There are many challenges involved in inverse model applications, e.g., the existence of non-identifiable parameters, the estimation of parameter uncertainties, and the quantification of the implications of these uncertainties on model predictions.
```

Excellent! That's what I want to do. Having now identified a hypothesis to test and a library that supports the task, it's time to blow the dust of my R Studio install and start modeling!

All code for this analysis can be found: <link to GitHub page here>.

# Model Definition

```
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

  times <- seq(0, 58, by=1)
  out <- ode(y = y, parms = pars, times = times, func = derivs)

  as.data.frame(out)
  }
```

# Parameter Optimization with FME

To test the hypothesis that the number of infected individuals _in China_ can be described using an SEIR model, we need to first make some (maybe arbitrary) assumptions about the eight initial conditions and model parameters in question. We need to set the values for susceptibles, exposed, infected, recovered, contact rate, transmission probability, infectious period, and latent period.

```
# population & virus assumptions
contact_rate = 10                     # number of contacts per day
transmission_probability = 0.10       # transmission probability
infectious_period = 5                 # infectious period
latent_period = 2                     # latent period

# beta (tranmission rate), gamma (recovery rate), and gamma (incubation rate)
beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period

# Ro (reproductive rate)
Ro = beta_value / gamma_value

# model hyperparameters
pars = c(beta = beta_value, gamma = gamma_value, delta = delta_value)
```

We then define a cost function that relates these model assumptions to the model output, infected individuals (at time `t`). The FME library has a function `modcost` that takes the model output and updates the estimated cost.

```
SEIRcost <- function(pars) {
  out <- SEIR(pars = pars)
  return(modCost(model = out, obs = ChinaInfected, err = "sd"))
  }
```

## Sensitivity Analysis
FME contains a handy `sensFun` method that estimates the sensitivity of the model output relative to the parameter values. When applied to observed data, these sensitivity functions help identify important parameters given the absolute magnitudes of the sensitivity summary values. We can even use this method to rank the importance of parameters on the number of infected individuals in China.

<insert sensitivity graph here>

# TODO:
Describe the sensitivity graph here...

# Fitting the Model to Data

The FME function `modFit` does a lot of sophisticated optimization that I don't have the mental bandwidth to dive into. Still, lucky for me their API is as simple as:

```
modFit(f = SEIRcost, p = pars)
```

The `modFit` function will optimize the parameters, `pars,` to the data, `ChinaInfected,` using the `SEIRCost` function, which traces back to the `SEIR` model definition defined above.

<insert sensitivity graph here>

Here's a graph that shows the actual data for infected individuals in China over time, the initial curve based on some arbitrary assumptions, and the fit curve _after parameter optimization_. The optimal parameter values can be retrieved with `coef`:

```
coef(Fit)
```

```
beta     gamma     delta
0.499    0.071     0.294
```

# TODO:
Interpret and describe optimal model parameters here...

# TODO: explore contact rate and transmission probability relationship, assuming beta == 0.499
<insert contact rate v transmission probability graph here>

# Bayesian Methods

What's really cool about the `FME` library are the Bayesian methods that can be used to derive data-dependent probability distributions of parameters.

## MCMC

The function `modMCMC` implements a Markov chain Monte Carlo (MCMC) method to sample from probability distributions by constructing a Markov chain with the desired distribution.... need to wordsmith... The API to perform this task is very user friendly:

```
MCMC <- modMCMC(f = SEIRcost, p = Fit$par, niter = 2000, jump = cov0,
                var0 = var0, wvar0 = 0.1, updatecov = 50)
```

Calling `pairs(MCMC, nsample = 1000)` on the MCMC object yields a pair plot of parameter-dependent distributions which shows the pairwise relationship in the upper panel, the correlation coefficients in the lower
panel, and the marginal distribution for each parameter, represented by a histogram, on the
diagonal.

<insert pairwise comparison of hyperparameters>

The strong relation between parameters beta (transmission rate) and delta (incubation rate) makes sense. Consider how the SEIR derivatives are defined, specifically, the rate at which susceptibles are exposed to the virus.

# TODO:
State benefit of Bayesian estimation of parameters - marginal distributions for possible values - and the correlation between values...

## Parameter Uncertainty (w.r.t Infected Individuals)

```
sensRange(func = SEIR, parms = pars, parInput = MCMC$par)
```

# Infected individuals sensitivity based on parameter distribution
plot(summary(sR), xlab = "time")

<insert sensitivity range graph based on  MCMC>

# TODO: interpret the sensitivity graph

# Estimated Reproductive Rate of COVID-19

As defined in the SEIR Model section above, the reproductive rate `R0` is defined as the rate of transmission over the rate of recovery. Equivalently it is the ratio of beta to gamma.

# TODO: R0 analysis and commentary

# Conclusion

# TODO: write a conclusion - yes, the infected population in China can be modeled relatively will with an SEIR model, and we can use an inverse modeling approach to identify a range of optimal, correlated parameters. Note interesting finds and recommendations for the next steps...
