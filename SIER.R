# SEIR Model with birth and death rate
require(ggplot2)
D <- 5
IN <- 10
N0 <- 7780000
parms <- c(
  b = 1/(60 * 365), # birth rate
  mu = 1/(60 * 365), # death rate
  gamma = 1/D, # The duration of infection
  beta0 = 3.6,
  sigma = 1/IN,
  N = N0,
  rho = .66,
  alpha = .11,
  delta = 0.01
)

# We need to set the inital states
y <- c(
  S = 0.065*N0,
  E = 20,
  I = 20,
  R = 0) 
deltat <- 1
maxout = 365*50

sier <- \(t, y, parms) {
  with(c(as.list(y), parms), {
    if (t >= 8401) {
      rho <- rho
    } else {
      rho <- 0
    }
    beta <- beta0 * (1 + alpha * cos(2 * pi * (t / 365)))
    dSdt <- b * N * (1 - rho) - mu * S - beta * S * I / N
    dEdt <- beta * S * I / N - sigma * E - mu * E
    dIdt <- sigma * E - gamma * I - mu * I - delta * I
    dRdt <- b * rho + gamma * I - mu * R
    return(list(c(dSdt, dEdt, dIdt, dRdt)))
  })
}

tt <- lsoda(
  y = y,
  times = seq(deltat, maxout, by = deltat),
  func = sier,
  parms = parms
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

