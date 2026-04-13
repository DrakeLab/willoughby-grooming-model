# Resource analysis 

# Sensitivity analysis 
set.seed(1993)
# libraries 
library(parallel) # for running in parallel 
library(purrr) 
library(scales) # heatmap plot

# DEFINE FUNCTIONS
# function for birth rate with Allee effect (theta)
birthrate <- function(x, b0, theta, K){
  b <- ifelse(x == 0, 0, b0*x*(x/(theta+x))) # baseline births with Allee
  b_dd <- b*x*(1 - x/ K) # birth decreases at high densities
  b_wA_dd <- pmax(0, b_dd)
  b_wA_dd # birth rate with Allee and density dependence
}

# function for immigration rate (iota)
# baseline iota, and resource linked increase
immigrationrate <- function (x, iota_baseline, iota_r, scenario){ 
  iota = iota_baseline + scenario * iota_r   
  iota
}

# function for programmed grooming 
# mid-grooming at a baseline level when resources are at 30 
groomrate_prg <- function(pop.size, groom_prg_bl, resource_scenario, groom_resource_link){
  death_rate <- pmin(pop.size, groom_prg_bl * resource_scenario * groom_resource_link)
  death_rate
}

# function for stimulus grooming deathrate based on pop size 
# mid-grooming at a baseline level when resources are at 30 
groomrate_stm <- function(pop.size, groom_stim_bl, resource_scenario, groom_resource_link){
  death_rate <- pop.size*( (groom_stim_bl * resource_scenario * groom_resource_link) /100)
  death_rate
}

# Function to get cumulative groom rate from programmed and stimulus grooming 
get_groom_rate <- function(pop.size, groom_prg_bl, groom_stim_bl, 
                           resource_scenario, groom_resource_link) {
  prg_rate = groomrate_prg(pop.size, groom_prg_bl, resource_scenario, groom_resource_link)
  stm_rate = groomrate_stm(pop.size, groom_stim_bl, resource_scenario, groom_resource_link) 
  groom_rate =  pmin(pop.size, prg_rate + stm_rate)
  groom_rate
}


# Function for density dependent deaths # old version
deathrate_dd <- function(pop.size, dmax = 1, x_mid = 50, s = 20) {
  dmax / (1 + exp(-(pop.size - x_mid) / s))
}

# set initial conditions parameters for RESOURCE analysis 
## birth rate is standard across all scenarios, right now as 37.5 flea births per year (0. per day)
## we are going to vary parameters for:  
# theta_vals <- seq(0, 20, by = 2) # will do just 8 

# iota_bl_vals <- seq(0, 10, by = 0.5) # will do just 5

# will set iota_r for 0 and 0.1
iota_r_vals <- c(0,0.1) 
scenario_vals <- seq(1,100, by = 1) # all resource changes 
groom_prg_bl_vals <- c(0,20) # grooming baseline [two scenarios: no programmed or baseline of 20]

groom_stim_bl_vals <- c(0, 20) # grooming baseline  [two scenarios: no stimulus or baseline of 20]

# construct full matrix of theta and iota options 
ti_axes <- expand.grid(
  #theta = theta_vals, # for now doing one theta value
  #iota_baseline = iota_bl_vals,  # fluctuating the baseline value of iota
  iota_r = iota_r_vals,  # for now doing one iota resource link value
  scenario = scenario_vals, # for now, doing one scenario of mid resource (r = 50)
  groom_prg_baseline = groom_prg_bl_vals, 
  groom_stm_baseline = groom_stim_bl_vals 
)

# Construct parameter set
param_grid <- ti_axes
param_grid$b0 <- log(78.75)/89 # constant [log(female fecundity per lifetime)/generation time]
param_grid$initial_pop <- 1
param_grid$theta <- 8
param_grid$iota_baseline <- 5
param_grid$groom_resource_link <- 1/50 

# remove combination of strategy for now
param_grid <- param_grid  %>%
  dplyr::filter(!(groom_prg_baseline==20 & groom_stm_baseline ==20))
# remove ir = 0, g = 0 as already have 100 sims of that 
param_grid <-  param_grid  %>%
  dplyr::filter(!(iota_r==0.0 & groom_stm_baseline ==0 & groom_prg_baseline == 0))
# remove values of r all odd numbers from resources 
param_grid <- param_grid %>% 
  dplyr::filter(scenario %% 2 == 0)

# set duration 
# Individual Experiment for Sensitivity 
num_simulations <- 90  # Number of times to run the simulation
max.time <- 730 # 2 year max for a squirrel  

# Define the task grid: 
task_grid <- expand.grid(
  param_id = 1:nrow(param_grid),
  sim_id  = 1:num_simulations)

# Define running the whole model as a function 
run_one_parameter_set <- function(task_i){
  
  param_id <- task_grid[task_i, 1]
  sim_id   <- task_grid[task_i, 2]
  
  message("Parameter set ", param_id, ", Simulation", sim_id)
  
  parms <- param_grid[param_id, ]
  
  events <- 0
  time <- 0
  x <- parms$initial_pop
  
  out_list <- list()
  
  while (time[length(time)] <= max.time) {
    
    events <- events + 1
    current.size <- x[length(x)]
    current.time <- time[length(time)]
    
    b <- birthrate(current.size, parms$b0, parms$theta, 50) # 50 as K
    g <- get_groom_rate(current.size, 
                        parms$groom_prg_baseline,
                        parms$groom_stm_baseline, 
                        parms$scenario, parms$groom_resource_link)
    #message("Current size", current.size, ", groom rate", g) # mid function checks
    d <- deathrate_dd(current.size) 
    iota <- immigrationrate(current.size, 
                            parms$iota_baseline, parms$iota_r, 
                            parms$scenario)
    
    total.rate <- b + g + d + iota
    
    if (is.na(total.rate) || total.rate <= 0) break
    
    increment.time <- -log(1 - runif(1)) / total.rate
    
    event_probs <- c(d, g, b, iota) / total.rate
    event_outcome <- sample(
      c("death","groomed","birth", "immigration"),
      size = 1,
      prob = event_probs
    )
    
    change <- ifelse(event_outcome %in% c("death", "groomed"), -1, 1)
    
    new_size <- max(0, current.size + change)
    new_time <- current.time + increment.time
    
    out_list[[events]] <- data.frame(
      event = events,
      time = new_time,
      population = new_size,
      event_type = event_outcome,
      sim_id = sim_id,
      param_id = param_id
    )
    
    time <- c(time, new_time)
    x <- c(x, new_size)
  }
  
  if (length(out_list) == 0) {
    data.frame(
      event = numeric(),
      time = numeric(),
      population = numeric(),
      event_type = character(),
      sim_id = numeric(),
      param_id = numeric()
    )
  } else {
    do.call(rbind, out_list)
  }}

# Now run for the full simulation across the 10 clusters 


# Now run for the full simulation across the 10 clusters 
system.time({
  
  cl <- makeCluster(10)
  
  # Move libraries to the clusters 
  clusterEvalQ(cl, library(tidyverse))
  
  # Move main environment items to each of the clusters
  clusterExport(cl, 
                c("param_grid","max.time", "task_grid", 
                  "birthrate", "deathrate_dd", "immigrationrate",
                  "get_groom_rate", "groomrate_prg", "groomrate_stm"
                ))
  
  simulation <- parLapply(cl, 1:nrow(task_grid), run_one_parameter_set)
  
  on.exit(stopCluster(cl), add = TRUE)
})
saveRDS(simulation, "data/sim_4gs_T8_BI5_allR_412_90-ss.Rdata") # these are full model outputs

# Save all the parameters used 

output_params <- param_grid
# create id col 
output_params$param_id = 1:nrow(output_params)
# merge task grid and params 
output_params <- left_join(task_grid, output_params, by = "param_id")

# Sample to just get end of day from each one 
# function to calculate daily snapshot 

end_of_day_flea_pop <- function(sim){
  
  init_row <- data.frame(  # pre load initial conditions simulations 
    event = 0,
    time = 0,
    population = 1,
    event_type = "initial",
    sim_id = sim$sim_id[1],
    param_id = sim$param_id[1]
  )
  
  sim %>% # Apply operations to the dataset
    bind_rows(init_row, .) %>%
    mutate(time_unit = floor(time)) %>% # Convert time to integer units
    group_by(time_unit) %>%
    slice_tail(n = 1) %>% # Get the last observation per unit
    ungroup() %>%
    select(time_unit, population) %>% # Keep relevant columns
    
    # create full sequence
    tidyr::complete(time_unit = 0:780) %>%
    
    # carry previous value forward
    tidyr::fill(population, .direction = "down") %>% 
    dplyr::filter(time_unit != 0) %>% # get rid of the initial seeding
    pull(population) # only get a vector of the population sizes 
}


# run for each simulation 
system.time({
  # add eods to output
  output_params <- output_params %>%
    mutate(
      eods = map(simulation, end_of_day_flea_pop),
      day_180 = map_dbl(eods, ~ if (length(.x) >= 180) .x[180] else NA_real_),
      day_780 = map_dbl(eods, ~ if (length(.x) >= 780) .x[780] else NA_real_)
    )
}) 

saveRDS(output_params, file = "data/resource_params_4gs_T8_BI5_allR_412_90-ss.rds")
