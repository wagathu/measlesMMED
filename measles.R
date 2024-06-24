
# Importing the packages --------------------------------------------------

pacman::p_load(data.table, tidyverse, purrr, deSolve, ggalt, zoo, lubridate, patchwork)

# Importing data ----------------------------------------------------------

df <- fread('data/londonMeaslesWeekly.csv')


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

# Plotting for the measles inicidence -------------------------------------

pltMeaslesCases <- df2 |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = cases)) +
  geom_vline(aes(xintercept = 1960),
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
  geom_vline(aes(xintercept = 1960),
             linetype = 2,
             col = 'red') +
  theme_classic() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.text = element_text(color = 'black'),
    axis.title = element_text(color = 'black')
  ) +
  labs(x = 'Year', y = 'Measles cases')
pltMeaslesCases + pltMeaslesCasesMonth 
ggsave(
  'images/pltMeaslesCasesMonth.png',
  width = 10,
  height = 6,
  dpi = 1e3,
  bg = NULL
)
