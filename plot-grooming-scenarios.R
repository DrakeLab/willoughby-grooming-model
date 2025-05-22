# plot different grooming scenarios 

# load packages 
library(tidyr)
library(ggplot2)
library(dplyr)

# Define parasite density range
x <- seq(0.01, 100, length.out = 1000)

# Define grooming scenarios 

## no grooming 
death_none <- function(x) rep(0, length(x))

## programmed to a certain number of fleas per day (G)
death_programmed <-  function(x, G = 10) {
  pmin(1, G / x)  # cap at 1 to avoid unrealistic death rates > 1
} 

## saturate to a % of fleas (e.g., stimulus driven) perG
death_stimulus <- function(x, d0 = 0, dmax = 1, perG  = 0.75) {
  k <- x * (1 - perG) / perG
  d0 + (dmax * x) / (k + x)
}

## costly, stimulus driven to a certain point, but then there is a cost 
death_costly <- function(x, max = 1, mu = 10, sigma = 10) {
  max * exp(-((x - mu)^2) / (2 * sigma^2))
}

# Evaluate each function
df <- data.frame(
  x = x,
  None = death_none(x),
  Programmed = death_programmed(x),
  Stimulus = death_stimulus(x),
  Costly = death_costly(x)
)


df_long <- pivot_longer(df, -x, names_to = "grooming scenario", values_to = "DeathRate")

# Plot
ggplot(df_long, aes(x = x, y = DeathRate, color = `grooming scenario`)) +
  geom_line(size = 1.2) +
  labs(
    title = "Per Capita Grooming-Induced Death Rate (Max = 1)",
    x = "ectoparasite population size",
    y = "death rate (`expression(mu)`)"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")

# 

# 1. Define parameter ranges
x_vals <- seq(1, 100, length.out = 100)  # parasite population size 
G_vals <- seq(1, 100, length.out = 100)   # grooming capacity (%) / per capita mortality rate

# 2. Create a grid of all combinations
grid <- expand.grid(
  pop_size = x_vals,
  death_rate = G_vals
)

# 3. Define total groomed for each scenario

# Programmed: can remove up to G fleas per day
grid$programmed <- pmin(grid$pop_size, grid$death_rate)

# Stimulus-driven: can remove up to %perG per day 
k_vals <-
grid$stimulus <- grid$pop_size*grid$death_rate/100 # thinking if saturating at 0.75% 

# Costly: Gaussian centered at mu = 50, sigma = 15
mu <- 50
sigma <- 15
efficiency <- exp(-((grid$pop_size - mu)^2) / (2 * sigma^2))
grid$costly <- grid$death_rate * efficiency

# 4. Reshape for plotting
plot_data <- pivot_longer(grid, cols = c("programmed", "stimulus", "costly"),
                          names_to = "scenario", values_to = "total_groomed")

plot_data$scenario <- factor(plot_data$scenario,
                             levels = c("programmed", "stimulus", "costly"))

# 5. Plot faceted heatmap
ggplot(plot_data, aes(x = pop_size, y = death_rate, fill = total_groomed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  facet_wrap(~scenario) +
  labs(
    x = "ectoparasite population size",
    y = "grooming efficiency",
    fill = "fleas groomed per day"
  ) +
  theme_minimal()
