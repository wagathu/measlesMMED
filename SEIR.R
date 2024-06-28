
#' This script simulated data from the SEIR model. We also did some model fitting
#' using the Maximum Likelihood Estimate which also present in the script
#' 
#' All the parameters used in the script (apart from the sinsidal function) were
#' adapted from the paper: A simple Model for Complex Dynamical Transitions in 
#' Epidemics' which can be found from the following link:
#' 
#' https://pubmed.ncbi.nlm.nih.gov/10650003/
#' 
#' The sinosidual function was obtain from the paper: 'A century of transitions in New 
#' York City's measles dynamics from the following link: 
#' 
#' https://royalsocietypublishing.org/doi/suppl/10.1098/rsif.2015.0024
#' 
#' The file was prepared by Brian Njuguna and Group 4 project members on MMED-2024
#' 

# Importing the packages --------------------------------------------------

pacman::p_load(data.table, deSolve, tidyverse, patchwork)

# Importing some data -----------------------------------------------------

demoDf <- fread('data/demoDF.csv')
myDat <- fread('data/myDat.csv')

# Collecting the parameter values -----------------------------------------

durationInfection <- 5 # The duration of infection
latentPeriod <- 8 # The latent period
N0 <- 6374980 # Initial population
#' The initial population to be used. The population was obtained from the 1944 total population

parms <- c(
  beta0 = 2.06,  # The transmission parameter
  alpha = 0.08  # The seasonal amplitude
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

#' Simulate the SEIR Model
#'
#' This function simulates the Susceptible-Exposed-Infectious-Recovered (SEIR) model over a specified 
#' time period. It uses a set of initial conditions and parameters to generate the model dynamics. 
#' Some parameters are internally defined within the function to facilitate parameter isolation 
#' during likelihood fitting.
#'
#' @param t A numeric vector of time points at which to evaluate the model.
#' @param y A numeric vector of initial conditions for the compartments: 
#' Susceptible (S), Exposed (E), Infectious (I), Recovered (R), and Cumulative Incidence (Ci).
#' @param parms A list of parameters for the SEIR model including:
#' \describe{
#'   \item{beta0}{Baseline transmission rate.}
#'   \item{alpha}{Amplitude of seasonal variation in transmission rate.}
#'   \item{latentPeriod}{Average duration of the latent period.}
#'   \item{durationInfection}{Average duration of the infectious period.}
#' }
#' @return A list of the derivatives of the SEIR compartments at each time point. The list contains:
#' \describe{
#'   \item{dSdt}{Rate of change of the Susceptible population.}
#'   \item{dEdt}{Rate of change of the Exposed population.}
#'   \item{dIdt}{Rate of change of the Infectious population.}
#'   \item{dRdt}{Rate of change of the Recovered population.}
#'   \item{dCidt}{Rate of change of the Cumulative Incidence.}
#' }
#' @examples
#' \dontrun{
#' times <- seq(0, 100, by = 1)
#' init <- c(S = 999, E = 1, I = 0, R = 0, Ci = 0)
#' parameters <- list(beta0 = 0.3, alpha = 0.1, latentPeriod = 5, durationInfection = 7)
#' seir(times, init, parameters)
#' }
#' @export
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
# Plotting ----------------------------------------------------------------

modelResults <- modelDf |>
  filter(date > as.Date('1950-01-01')) |> 
  ggplot() +
  geom_line(aes(x = date, y = weekly_incidence), col = 'red') +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black'),
    plot.title = element_text(color = 'black', hjust = .5)
  ) +
  labs(x = 'days', y = 'Measles cases') 
modelResults
ggsave(
  'images/birth0.01375.png',
  width = 10,
  height = 6,
  dpi = 1e3,
  bg = NULL
)

# Comparing to our own data
dfCompare <- modelDf |>
  data.frame() |> 
  mutate(cases = demoDf$cases)

pltCompareFirst <- dfCompare |>
  filter(date > as.Date('1950-01-01')) |> 
  ggplot(aes(x = date)) +
  geom_point(aes(y = cases), col = 'red') +
  geom_line(aes(y = weekly_incidence)) +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black'),
    plot.title = element_text(color = 'black', hjust = .5)
  ) +
  labs(x = 'Year', y = 'Weekly incidence')
pltCompareFirst
ggsave(
  'images/birth0.02.png',
  width = 10,
  height = 6,
  dpi = 1e3,
  bg = NULL
)

# Several birth rates -----------------------------------------------------

#' Simulate the SEIR Model for Different Birth Rates
#'
#' This function simulates the SEIR model for different birth rates. It takes in a 
#' set of parameters and a vector of birth rates, evaluates the SEIR model at
#' each birth rate, and returns the results as a list of data frames.
#'
#' @param parms A list of parameters for the SEIR model, including:
#' \describe{
#'   \item{beta0}{Baseline transmission rate.}
#'   \item{alpha}{Amplitude of seasonal variation in transmission rate.}
#'   \item{latentPeriod}{Average duration of the latent period.}
#'   \item{durationInfection}{Average duration of the infectious period.}
#' }
#' @param births A numeric vector of birth rates to evaluate.
#' @param ... Additional arguments passed to the model.
#' @return A list of data frames, each containing the results of the SEIR model 
#' simulation for a given birth rate. Each data frame includes:
#' \describe{
#'   \item{date}{Date of the observation (week ending).}
#'   \item{weekly_incidence}{Weekly incidence of the disease.}
#'   \item{S}{Number of susceptible individuals.}
#'   \item{E}{Number of exposed individuals.}
#'   \item{I}{Number of infectious individuals.}
#'   \item{R}{Number of recovered individuals.}
#'   \item{N}{Total population (S + E + I + R).}
#'   \item{P}{Prevalence of the disease (I / N).}
#' }
#' @examples
#' \dontrun{
#' parms <- list(beta0 = 0.3, alpha = 0.1, latentPeriod = 5, durationInfection = 7)
#' births <- c(0.01, 0.02, 0.03, 0.04, 0.05)
#' results <- diffBirths(parms, births)
#' }
#' @export
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
        mu = b / 365                                     # The death rate
        b = b / 365                                      # The birth rate
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
      .[, cases := demoDf$cases] %>%
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
individual <- wrap_plots(individualWk, nrow = 2)
individual
ggsave(
  'images/combinedWeeklyIncidence.png',
  width = 20,
  height = 18,
  dpi = 1e3,
  bg = NULL
)

# The SEIR Plot
wrap_plots(copy(dfs) |>
             _[ ,birth := case_when(
               birth == 'birth_0.01' ~ '0.01',
               birth == 'birth_0.02' ~ '0.02',
               birth == 'birth_0.03' ~ '0.03',
               birth == 'birth_0.04' ~ '0.04',
               TRUE ~ '0.05'
             )] %>%
             .[, .(date, S, E, I, R, birth)] %>% 
             .[, melt(
               .SD,
               id.vars = c('date', 'birth')
             )] %>%
             split(.$birth) %>%
             map2(., names(.),
               ~ .x |> 
                 ggplot(aes(x = date)) +
                 geom_line(aes(y = value, colour = variable)) +
                 theme_light() +
                 theme(
                   #axis.line = element_line(color = 'black'),
                   axis.text = element_text(color = 'black'),
                   axis.title = element_text(color = 'black'),
                   plot.title = element_text(color = 'black', hjust = .5)
                 ) + 
                 labs(x = 'Year', y = '', title = paste('Birth rate:', .y)) +
                 scale_y_continuous(labels = scales::number_format())
  ), guides = 'collect')

ggsave(
  'images/SEIR.png',
  width = 20,
  height = 18,
  dpi = 1e3,
  bg = NULL
)

# Comparing with the actual data ------------------------------------------


compares <-  copy(dfs) |>
  _[, birth := case_when(
    birth == 'birth_0.01' ~ '0.01',
    birth == 'birth_0.02' ~ '0.02',
    birth == 'birth_0.03' ~ '0.03',
    birth == 'birth_0.04' ~ '0.04',
    TRUE ~ '0.05'
  )] %>%
  split(.$birth) %>%
  map2(., as.character(names(.)), \(x, y) {
    x |>
      ggplot(aes(x = date)) +
      geom_point(aes(y = cases, col = '1', ), show.legend = TRUE) +
      geom_line(aes(y = weekly_incidence, col = '2', ), show.legend = TRUE) +
      theme_classic() +
      theme(
        axis.line = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        plot.title = element_text(color = 'black', hjust = .5, size = 14),
        legend.text = element_text(color = 'black', size = 14)
      ) +
      labs(x = 'Year',
           y = 'Weekly measles cases',
           title = paste('Birth rate:', y)) +
      scale_color_manual(
        values = c("2", "1"),
        labels = c("Actual measles cases", "Cases from model"),
        guide = guide_legend(override.aes = list(linetype = c(NA, 2))) 
        
      ) 
  })

leg <- as_ggplot(get_legend(compares))

pltCompare <- wrap_plots(compares, leg, nrow = 2, guides = 'collect')  &
  theme(legend.position = "bottom")
ggsave(
  'images/pltCompare.png',
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
