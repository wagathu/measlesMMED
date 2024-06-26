
# Importing the packages --------------------------------------------------

pacman::p_load(data.table, deSolve, tidyverse)

# Collecting the parameter values -----------------------------------------
durationInfection <- 5
latentPeriod <- 8
N0 <- 5e5
parms <- c(
  b = 1.10/(60 * 365),                         # The birth rate
  mu = 1/(60 * 365),                           # The death rate
  beta0 = 3.6,                                 # The transmission parameter
  alpha = .011,                                  # The amplitude 
  sigma = 1 / latentPeriod,                    # Progression rate from pre-infectious to Infectious
  delta = 0.01,                                # Disease induced deaths
  gamma = 1 / durationInfection                # Recovery rate
)

# The state variables -----------------------------------------------------

y <- c(S = .99 * N0,
       E = 0,
       I = .1 * N0,
       R = 0)

# Some more parameters of interest
deltat <- 1
maxout <- 365 * 50

# The SEIR model ----------------------------------------------------------

seir <- \(t, y, parms) {
  with(c(as.list(y), parms), {
    N = sum(y)
    births = b * N
    beta <- beta0 * (1 + alpha * cos(2 * pi * (t / 365)))
    beta <- beta
    dSdt = births - mu * S - beta * S * I / N
    dEdt = beta * S * I / N - mu * E - sigma * E
    dIdt = sigma * E - I * (mu + delta) - gamma * I
    dRdt = gamma * I - mu * R
    return(list(c(dSdt, dEdt, dIdt, dRdt)))
  })
}

# Solving the differential equations --------------------------------------

modelResults <- lsoda(
  y = y,
  parms = parms,
  func = seir,
  times = seq(deltat, maxout, by = deltat)
)

modelDf <- modelResults |> 
  as.data.frame() |> 
  mutate(N = S + E + I + R)

modellong <- modelDf |> pivot_longer(cols = -time)

# Plotting ----------------------------------------------------------------

modelResults <- subset(modellong, time > 100 & name == "I") |>
  ggplot() +
  geom_line(aes(x = time, y = value, color = name)) +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black')
  ) +
  labs(x = 'days', y = 'Measles cases') +
  scale_x_continuous()
modelResults

