# Code to recreate Drake and Lodge 2006 
# Authors: Anna Willoughby

rm(list = ls()) # clear the workspace
graphics.off() # close graphics windows

# load packages 
library(tidyverse)
library(deSolve)

# set the parameters
parameters <- c(
  r = 0.1, # growth rate 
  K = 600, # carrying capacity 
  lambda = 0.6, # birth rate 
  theta = 20, # Allee effect 
  iota = 0, # immigration rate (first let's have at 0)
  n_inv = 100, # nuisance population size 
  n_ext = 1,  # extinct population size
  alpha = 1 
)

# other equations
# n_star = (-b +/- sqrt((b^2-4ac))/2a)
# a = -r/k
# b = r*(k-theta)/k
# c = theta*(r-lambda)
# v(n) = alpha*n

# set up the system of differential equations-----------------------------------
# Eq 1: logistic growth
Eq1 <- function(t, state, parameters) {{
  with(as.list(c(state, parameters)), {{
    dndt <- r*n -(r / K)*(n^2)
    list(dndt)
  }})
}}

# Eq 2: logistic growth with Allee 
Eq2 <- function(t, state, parameters) {{
  with(as.list(c(state, parameters)), {{
    dndt <- r*n -(r / K)*(n^2) -  ((lambda*theta)/(theta + n))*n
    list(dndt)
  }})
}}

# Eq 3: logistic growth with immigration
Eq3 <- function(t, state, parameters) {{
  with(as.list(c(state, parameters)), {{
    dndt <- r*n -(r / K)*(n^2) + iota
    list(dndt)
  }})
}}


# Eq 4: logistic growth with Allee and immigration 
Eq4 <- function(t, state, parameters) {{
  with(as.list(c(state, parameters)), {{
    dndt <- r*n -(r / K)*(n^2) - ((lambda*theta)/(theta + n))*n + iota
    list(dndt)
  }})
}}

# first passage probability of attaining n_inv before n_ext

fpp = 1 - (integrate(exp(-2*integrate(Eq/alpha*n)), lower = n, upper = n_inv)) / (integrate(exp(-2*(Eq/alpha*n)), lower = n_ext , upper = n_inv ))

# Will the population invade when smaller than the nuisance threshold? (Figure 1)

# create vector of possible initial population sizes 
n_values <- seq(1,100, by = 1)

# solve for fpp with Eq3 
Eq = Eq3 
for(n in n_values){
  state = n 
  sol = as.data.frame(
    ode(y = state, func = fpp, parms = parameters)
  )
}

## SCRAPS BELOW 

### Figure 2: Stationary distribution of population sizes 

# set up a data frame to be populated with initial and mean population 
df <- data.frame(initial_population_n = integer(), mean_n = numeric())
# run a loop to calculate mean pop size 
for(n in n_values){
  state <- n
  times <- seq(0, 1000, by = 0.2)
  sol = as.data.frame(
    ode(y = state, times = times, func = Eq3, parms = parameters)
  )
  sol_df <- data.frame(n, mean(sol[[2]]))
  df <- rbind(df, sol_df)
}

df$m <- mean(logistic_allee_imm(df$n))
phi = 2*integrate() 

plot(sol)

sol = as.data.frame(
  ode(y = 100, times = times, func = Eq4, parms = parameters)
)


# nonlinear birth-death-immigration process