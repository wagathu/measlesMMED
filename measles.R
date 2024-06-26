
# Importing the packages --------------------------------------------------

pacman::p_load(data.table, tidyverse, purrr, deSolve, ggalt, zoo, lubridate, patchwork)

# Importing data ----------------------------------------------------------

df <- fread('data/londonMeaslesWeekly.csv')
demography <- fread('data/londonDemography.csv')

# Creating the weekly measles cases ---------------------------------------

df2 <- df |> 
  mutate(date = ymd(paste(year, month, day, sep = '-'))) |> 
  dplyr::select(date, cases)

dfMonth <- df %>%
  mutate(date = make_date(year, month, day),
         yearMonth = as.yearmon(date)
  ) |> 
  dplyr::select(yearMonth, cases) |> 
  group_by(yearMonth) |> 
  summarise(cases = sum(cases)) 


# Combining for demography ------------------------------------------------

demoDf <- df2 |> 
  mutate(year = year(date)) |> 
  merge(demography |> dplyr::select(year = Year, births = Births, pop = Pop, vax = Vax),
        by = 'year'
        ) |> 
  mutate(B_t = (births/pop)/52,
         vacc_rate = (vax/pop)/52
         )

modelDf <- demoDf |> 
  dplyr::select(date, B_t, vacc_rate) |> 
  mutate(date = ymd(date))

# Plotting for the measles inicidence -------------------------------------

pltMeaslesCases <- df2 |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = cases)) +
  geom_vline(aes(xintercept = 1967),
             linetype = 2,
             col = 'red') +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black')
  ) +
  labs(x = 'Year', y = 'Measles cases')
pltMeaslesCases
ggsave(
  'images/pltMeaslesCases.png',
  width = 10,
  height = 6,
  dpi = 1e3,
  bg = NULL
)


# Monthly cases
pltMeaslesCasesMonth <- dfMonth |>
  ggplot(aes(x = yearMonth)) +
  geom_line(aes(y = cases)) +
  geom_vline(aes(xintercept = 1967),
             linetype = 2,
             col = 'red') +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black')
  ) +
  labs(x = 'Year', y = 'Measles cases')
pltMeaslesCasesMonth 
ggsave(
  'images/pltMeaslesCasesMonth.png',
  width = 10,
  height = 6,
  dpi = 1e3,
  bg = NULL
)


difftime(as.Date('1967-01-01'), min(dfMonth$yearMonth), units = 'weeks')

