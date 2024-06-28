
#' This script loads the data for the weekly measles cases (London 1944 - 1994), yearly 
#' births in London (1944 - 1994) and the year population in London (1944 - 1994)
#' 
#' The script also plots the weekly and monthly measles incidence and then finally 
#' saves it to the /images folder in this repository

# Importing the packages --------------------------------------------------

pacman::p_load(data.table,
               tidyverse,
               purrr,
               ggalt,
               zoo,
               lubridate,
               patchwork)

# Importing data ----------------------------------------------------------

df <- fread('data/londonMeaslesWeekly.csv')
b <- fread('data/births.csv')
pop <- fread('data/population.csv')

# Creating the weekly measles cases ---------------------------------------

demography <- merge(b |>
                      mutate(year = as.integer(year))
                    , pop, by = 'year')

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
  merge(demography ,
        by = 'year'
        ) |> 
  mutate(date = ymd(date),
         P = cases/population,
         B_t = as.numeric(births)/population)
fwrite(demoDf, 'data/demoDF.csv', row.names = F)

myDat <- demoDf |>
  mutate(
    time = date,
    numPos = cases,
    numSamp = population,
    P = cases / population
  ) |>
  dplyr::select(time, numPos, numSamp, P)
fwrite(myDat, 'data/myDat.csv', row.names = F)

# Plotting for the measles inicidence -------------------------------------

pltMeaslesCases <- df2 |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = cases)) +
  geom_vline(aes(xintercept = as.Date('1967-01-01')),
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
  geom_vline(aes(xintercept = zoo::as.yearmon('1967-01-01')),
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
