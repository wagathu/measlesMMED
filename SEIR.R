
# Importing the packages --------------------------------------------------

pacman::p_load(data.table, deSolve, tidyverse, patchwork)
source('measles.R')

# Collecting the parameter values -----------------------------------------
durationInfection <- 5
latentPeriod <- 8
N0 <- 6374980
# parms <- c(
#   b = 1.10/(60 * 52),                            # The birth rate
#   mu = 1/(60 * 52),                              # The death rate
#   beta0 = 3.6,                                   # The transmission parameter
#   alpha = .11,                                   # The amplitude 
#   sigma = 1 / latentPeriod,                      # Progression rate from pre-infectious to Infectious
#   delta = 0.001/52,                              # Disease induced deaths
#   gamma = 1 / durationInfection                  # Recovery rate
# )

parms <- c(
  beta0 = 2.6,  # The transmission parameter
  alpha = 0.08  
)


# The state variables -----------------------------------------------------

y <- c(
  S = .1 * 6374898,
  E = .001* 6374898,
  I = .001* 6374898,
  R = 0,
  Ci = 0
)

# Some more parameters of interest
deltat <- 1
times <- seq(from = as.Date('1944-01-03'), to = as.Date('1994-12-31'), by = 1)

# The SEIR model ----------------------------------------------------------

seir <- \(t, y, parms) {
  with(c(as.list(y), parms), {
    N = sum(y)
    
  # Fixed parameters
    delta = 0
    mu = .02/365                                     # The death rate
    b = .02/365                                      # The birth rate
    sigma = 1 / latentPeriod                       # Progression rate from pre-infectious to Infectious
    gamma = 1 / durationInfection                  # Recovery rate
    
    births = b * N
    beta <- beta0 * (1 + alpha * cos(2 * pi * (t)))
    dSdt = births - mu * S - beta * S * I / N
    dEdt = beta * S * I / N - mu * E - sigma * E
    dIdt = sigma * E - I * (mu + delta) - gamma * I
    dRdt = gamma * I - mu * R
    dCidt = sigma * E
    return(list(c(dSdt, dEdt, dIdt, dRdt, dCidt)))
  })
}

# Solving the differential equations --------------------------------------

modelResultsDf <- lsoda(
  y = y,
  parms = parms,
  func = seir,
  times = seq(1, length(times), by = 1),
)
modelDf <- modelResultsDf |> 
  as.data.frame() |> 
  mutate(N = S + E + I + R,
         P = I/N,
         time = times
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

# Several birth rates -----------------------------------------------------

diffBirths <- \(parms, births, ...) {
  map2(births, as.character(births), \(b, name) {
    # The parameters for the model
    parms <- c(b, parms)
    # The initial states
    y <- c(
      S = .1 * 6374898,
      E = .001 * 6374898,
      I = .001 * 6374898,
      R = 0,
      Ci = 0
    )
    
    # The function for the model
    seir <- \(t, y, parms) {
      with(c(as.list(y), parms), {
        N = sum(y)
        
        # Fixed parameters
        delta = 0
        mu = .02 / 365                                     # The death rate
        b = .02 / 365                                      # The birth rate
        sigma = 1 / latentPeriod                         # Progression rate from pre-infectious to Infectious
        gamma = 1 / durationInfection                    # Recovery rate
        
        births = b * N
        beta <- beta0 * (1 + alpha * cos(2 * pi * (t)))
        dSdt = births - mu * S - beta * S * I / N
        dEdt = beta * S * I / N - mu * E - sigma * E
        dIdt = sigma * E - I * (mu + delta) - gamma * I
        dRdt = gamma * I - mu * R
        dCidt = sigma * E
        return(list(c(dSdt, dEdt, dIdt, dRdt, dCidt)))
      })
    }
    cat(crayon::bold(crayon::cyan('Runing the model with birth rate:')), name, '\n')
    
    modelResultsDf <- lsoda(
      y = y,
      parms = parms,
      func = seir,
      times = seq(1, length(times), by = 1),
    )
    modelDf <- modelResultsDf |>
      as.data.frame() |>
      mutate(time = times) |>
      mutate(week_end = ceiling_date(time, "week", change_on_boundary = FALSE) - 2) |>
      group_by(date = week_end) %>%
      summarize(
        weekly_incidence = last(Ci) - first(Ci),
        .groups = 'drop',
        across(c(S, E, I, R), ~ sum(.)),
        N = S + E + I + R,
        P = I / N
      )
    return(modelDf)
  })
}

births <- c(seq(.01, .05, by = .01))
dfs <- diffBirths(parms = parms, births = births) %>%
  setNames(c(paste0('birth_', as.character(births)))) %>%
  map2(., names(.), \(x, y) {
    x |>
      data.table() %>%
      .[, birth := y] %>%
      .[date > as.Date('1950-12-31'),]
  }) |> rbindlist()

# Plotting the weekly measles cases combined
pltWeeklyIncidence <- copy(dfs) |> 
  _[ ,birth := case_when(
    birth == 'birth_0.01' ~ '0.01',
    birth == 'birth_0.02' ~ '0.02',
    birth == 'birth_0.03' ~ '0.03',
    birth == 'birth_0.04' ~ '0.04',
    TRUE ~ '0.05'
  )] |> 
  ggplot(aes(x = date)) +
  geom_line(aes(y = weekly_incidence, colour = birth)) +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black'),
    plot.title = element_text(color = 'black', hjust = .5)
  ) +
  labs(x = 'Year', y = 'Weekly measles cases', col = 'birth rate')
pltWeeklyIncidence
ggsave(
  'images/pltWeeklyIncidence.png',
  width = 10,
  height = 6,
  dpi = 1e3,
  bg = NULL
)

# Plotting the weekly measles cases but individually
individualWk <- copy(dfs) |>
  _[ ,birth := case_when(
    birth == 'birth_0.01' ~ '0.01',
    birth == 'birth_0.02' ~ '0.02',
    birth == 'birth_0.03' ~ '0.03',
    birth == 'birth_0.04' ~ '0.04',
    TRUE ~ '0.05'
  )] %>%
  split(.$birth) %>%
  map2(., as.character(names(.)),
    \(x,y) { 
      x |> 
      ggplot(aes(x = date)) +
      geom_line(aes(y = weekly_incidence)) +
      theme_classic() +
      theme(
        axis.line = element_line(color = 'black'),
        axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black'),
        plot.title = element_text(color = 'black', hjust = .5)
      ) + 
        labs(x = 'Year', y = 'Weekly measles cases', title = paste('Birth rate:', y)) 
    }
  )
individual <- wrap_plots(individual, nrow = 2)
individual
ggsave(
  'images/combinedWeeklyIncidence.png',
  width = 20,
  height = 18,
  dpi = 1e3,
  bg = NULL
)

# Calculating the negative likelihood -------------------------------------

nlikelihood <- \(par, obsDat) {
  simDat <- lsoda(
    y = y,
    parms = par,
    func = seir,
    times = seq(1, nrow(demoDf), by = 1)
  )
  simDat <- simDat |>
    as.data.frame() |>
    mutate(P = I / (S + E + I + R))

  nlls <- -dbinom(x = obsDat$cases,
                  size = obsDat$population,
                  prob = simDat$P,
                  log = T)
  return(sum(nlls))
}


# Solving using the optim function
optim.vals <- optim(par = parms, #isolated
                    fn = nlikelihood,
                    obsDat = data.frame(demoDf),
                    control = list(trace = 3, maxit = 150),
                    method = "Nelder-Mead")

# simDat <- lsoda(
#   y = y,
#   parms = optim.vals$par,
#   func = seir,
#   times = seq(1, nrow(demoDf), by = 1)
# )