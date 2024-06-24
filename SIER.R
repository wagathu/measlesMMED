# SEIR Model with birth and death rate
require(ggplot2)
D <- 5
IN <- 10
N0 <- 7780000
parms <- c(
  b = 1/(10 * 365), # birth rate
  mu = 1/(60 * 365), # death rate
  gamma = 1/D, # The duration of infection
  beta0 = 3.6,
  sigma = 1/IN,
  N = N0,
  rho = .2,
  vacc.yr = 365 * 2,
  alpha = .11
)

# We need to set the inital states
y <- c(
  S = 0.065*N0,
  E = 20,
  I = 20,
  R = 0) 
deltat <- 0.1
maxout = 365*5

sier <-\(t, parms,y) {
  with(c(as.list(y), parms), {
   rho <-  ifelse(t >= vacc.yr, rho, 0)
    beta  = beta0 * (1 + alpha * cos(2 * pi * t))
    dSdt <- b * N * (1 - rho)  - mu * S - beta * S * I/N
    dEdt <- beta * S * I/N - sigma * E -  mu * E
    dIdt <- sigma * E  - gamma * I -  mu * I
    dRdt <- b * rho + gamma * I - mu * I
    return(list(c(dSdt,dEdt, dIdt, dRdt)))
  })
}

tt <- lsoda(
  y = y,
  times = seq(deltat, maxout, by = deltat),
  func = sier,
  parms = parms
) |> 
  as.data.frame()


tt |> 
  ggplot() +
  geom_line(aes(x = time, y = I)) +
  theme_light() +
  annotate(
    geom = "segment",
    x = 365 * 2,
    y = 15000,
    xend = 365 * 2,
    yend = 50000,
    arrow = arrow(length = unit(0.15, "cm"), ends = "first"),
    color = "#2b8cbe"
  ) +
  annotate(
    geom = "text",
    x = 365 * 2,
    y = 20000,
    label = "Start of vaccination",
    vjust = -30,
    hjust = 0.5,
    color = "#2b8cbe",
    size = 3.5
  )

