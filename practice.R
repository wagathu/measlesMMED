# SEIR Model with birth and death rate
require(ggplot2);require(deSolve)
D <- 5
IN <- 8
N0 <- 500000
latent <- 8
parms <- c(
  gamma = 1 / D, # duration of infection
  beta0 = 3.6,
  sigma = 1 / latent,
  rho = 0.6,
  alpha = 0.11,
  alpha2 = .02,
  delta = 0.01,
  sigma_se = 0.1, # standard deviation for gamma white noise
  cc = 2e-3, # cohort parameter for seasonal birth rate variation
  Ve = 0.9, # vaccine efficacy
  iota = .31,
  B = 1 / (60 * 52), # annual birth rate
  mu = 1 / (60 * 52) # death rate
)
# We need to set the inital states
y <- c(
  S = 0.0311 * N0,
  E = 5e-4 * N0,
  I = 2.66e-4 * N0,
  R = 0.968 * N0
)
deltat <- 7
times <- seq.Date(from = as.Date(min(demoDf$date)), to = as.Date(max(demoDf$date)), by = deltat)
sier <- \(t, y, parms) {
  with(c(as.list(y), parms), {
    
    date2 <- times[ceiling(t / deltat) + 1]
    year <- as.numeric(format(date2, "%Y"))
    day_of_year <- as.numeric(format(date2, "%j"))
    week_of_year <- as.numeric(format(date2, "%U"))
    
    # Extract B_t and vacc_rate from demoDf for the current date
    current_data <- demoDf[date == date2,]
    B_t <- current_data$B_t
    vacc_rate <- .7
    
    # Adjust birth rate for the first week of September
    if (week_of_year == 35) {
      birth_rate <- cc * B / deltat + (1 - cc) * B_t
    } else {
      birth_rate <- (1 - cc) * B_t
    }
    beta <- beta0 * (1 + alpha2 * cos(2 * pi * (t / 365)))
    
    # Incorporate stochasticity
    births <- rpois(1, birth_rate * deltat)
    vaccination <- rpois(1, vacc_rate * Ve * deltat)
    N = sum(y)
    lambda <- beta * S /N * rgamma(1, shape = 1, scale = sigma_se)
    dSdt <- births - lambda * (I + iota)^alpha - mu * S - vaccination
    dEdt <- lambda * (I + iota)^alpha - sigma * E - mu * E
    dIdt <- sigma * E - gamma * I - mu * I
    dRdt <- gamma * I - mu * R
    return(list(c(dSdt, dEdt, dIdt, dRdt)))
  })
}

tt <- lsoda(
  y = y,
  times = seq(0, 52 * 50, by = deltat),
  func = sier,
  parms = parms,
  rtol = 1e-15,
  atol = 1e-15
) |> 
  as.data.frame()


modelResults <- tt |>
  ggplot() +
  geom_line(aes(x = time, y = I)) +
  geom_vline(aes(xintercept = 8401),
             linetype = 2,
             col = 'red') +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black')
  ) +
  labs(x = 'Year', y = 'Measles cases')
modelResults
# ggsave(
#   'images/modelResults.png',
#   width = 10,
#   height = 6,
#   dpi = 1e3,
#   bg = NULL
# )
# 
# 

