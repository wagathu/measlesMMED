
# Importing the packages --------------------------------------------------

pacman::p_load(data.table, deSolve, tidyverse)
source('measles.R')

# Collecting the parameter values -----------------------------------------
durationInfection <- 5/7
latentPeriod <- 8/7
N0 <- 6374980
parms <- c(
  b = 1.10/(60 * 52),                         # The birth rate
  mu = 1/(60 * 52),                           # The death rate
  beta0 = 3.6,                                 # The transmission parameter
  alpha = .11,                                 # The amplitude 
  sigma = 1 / latentPeriod,                    # Progression rate from pre-infectious to Infectious
  delta = 0.001/52,                              # Disease induced deaths
  gamma = 1 / durationInfection                # Recovery rate
)

# The state variables -----------------------------------------------------

y <- c(S = .0645 * 6374898,
       E = 0,
       I = 82,
       R = 0)

# Some more parameters of interest
deltat <- 1
#maxout <- 365 * 50

# The SEIR model ----------------------------------------------------------

seir <- \(t, y, parms) {
  with(c(as.list(y), parms), {
    N = sum(y)
    births = b * N
    beta <- beta0 * (1 + alpha * cos(2 * pi * (t / 52)))
    dSdt = births - mu * S - beta * S * I / N
    dEdt = beta * S * I / N - mu * E - sigma * E
    dIdt = sigma * E - I * (mu + delta) - gamma * I
    dRdt = gamma * I - mu * R
    return(list(c(dSdt, dEdt, dIdt, dRdt)))
  })
}

# Solving the differential equations --------------------------------------

modelResultsDf <- lsoda(
  y = y,
  parms = parms,
  func = seir,
  times = seq(1, nrow(demoDf), by = 1)
)
modelDf <- modelResultsDf |> 
  as.data.frame() |> 
  mutate(N = S + E + I + R,
         time = demoDf$date
         )

# Plotting ----------------------------------------------------------------

# modelResults <- subset(modellong, time > 100) |>
#   ggplot() +
#   geom_line(aes(x = time, y = value, color = name)) +
#   theme_classic() +
#   theme(
#     axis.line = element_line(color = 'black'),
#     axis.text = element_text(color = 'black'),
#     axis.title = element_text(color = 'black')
#   ) +
#   labs(x = 'days', y = 'Measles cases') +
#   scale_x_log10()
# modelResults

modelResults <- modelDf |>
  filter(time > as.Date('1950-01-01')) |> 
  ggplot() +
  geom_line(aes(x = time, y = I), col = 'red') +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black'),
    plot.title = element_text(color = 'black', hjust = .5)
  ) +
  labs(x = 'days', y = 'Measles cases') 
modelResults
# ggsave(
#   'images/birth0.01375.png',
#   width = 10,
#   height = 6,
#   dpi = 1e3,
#   bg = NULL
# )

# Comparing to our own data
dfCompare <- modelDf |>
  data.frame() |> 
  mutate(date = time) |> 
  dplyr::select(-time) |> 
  merge(demoDf |> select(date,cases), by = 'date')

pltCompare <- dfCompare |>
  filter(date > as.Date('1950-01-01')) |> 
  ggplot(aes(x = date)) +
  geom_point(aes(y = cases), col = 'red') +
  geom_line(aes(y = I)) +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black'),
    plot.title = element_text(color = 'black', hjust = .5)
  ) 
pltCompare
# ggsave(
#   'images/birth0.01375.png',
#   width = 10,
#   height = 6,
#   dpi = 1e3,
#   bg = NULL
# )

# Calculating the negative likelihood -------------------------------------

# nlikelihood <- \(par, obsDat) {
#   simDat <- lsoda(
#     y = y,
#     parms = par,
#     func = seir,
#     times = seq(1, nrow(demoDf), by = 1)
#   )
#   simDat <- simDat |>
#     as.data.frame() |> 
#     mutate(P = I / (S + E + I + R))
#   
#   nlls <- -dbinom(x = obsDat$cases,
#                   size = obsDat$population,
#                   prob = simDat$P,
#                   log = T)
#   return(sum(nlls))
# }
# 
# # Solving using the optim function
# optim.vals <- optim(par = parms, #isolated
#                     fn = nlikelihood,
#                     obsDat = data.frame(demoDf),
#                     control = list(trace = 3, maxit = 150),
#                     method = "Nelder-Mead")
# 
# simDat <- lsoda(
#   y = y,
#   parms = optim.vals$par,
#   func = seir,
#   times = seq(1, nrow(demoDf), by = 1)
# )