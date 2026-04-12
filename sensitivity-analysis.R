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

# set initial conditions parameters for SENSITIVITY analysis 
## birth rate is standard across all scenarios, right now as 37.5 flea births per year (0. per day)
## we are going to vary parameters for:  
theta_vals <- seq(0, 20, by = 2) # 11 different values

iota_bl_vals <- seq(0, 10, by = 0.5) # 21 different values Keep for all the same

# iota_r_vals will be set to 0, no resource changes
# scenario_vals <- seq(1,100, by = 1) # for now, we will explore at a mid resource level (r = 30)
groom_prg_bl_vals <- c(0,20) # grooming baseline [two scenarios: no programmed or baseline of 20]

groom_stim_bl_vals <- c(0, 20) # grooming baseline  [two scenarios: no stimulus or baseline of 20]

# construct full matrix of theta and iota options 
ti_axes <- expand.grid(
  theta = theta_vals, # for now doing one theta value
  iota_baseline = iota_bl_vals,  # fluctuating the baseline value of iota
  #iota_r = iota_r_vals,  # for now doing one iota resource link value
  #scenario = scenario_vals, # for now, doing one scenario of mid resource (r = 50)
  groom_prg_baseline = groom_prg_bl_vals, 
  groom_stm_baseline = groom_stim_bl_vals 
)

# Construct parameter set
param_grid <- ti_axes
param_grid$b0 <- log(78.75)/89 # constant [log(female fecundity per lifetime)/generation time]
param_grid$initial_pop <- 1
param_grid$iota_r <- 0
param_grid$scenario <- 50
param_grid$groom_resource_link <- 1/50 

# set duration 
# Individual Experiment for Sensitivity 
num_simulations <- 60  # Number of times to run the simulation
# max.events <- 1000 # maximum events
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
saveRDS(simulation, "data/sim_4gs_allT_allBI_mg50_411_100.Rdata") # these are full model outputs 

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

saveRDS(sens_df, file = "data/sensitivity_params_4gs_allT_allBI_mg50_412_100.rds")

# load all sims 
sens_df <- readRDS(file = "data/sensitivity_params_4gs_allT_allBI_mg50_412_100.rds")

# calculate k and var:m 
sens_df_long <- sens_df %>%
  pivot_longer(cols = c(day_180, day_780),
               names_to = "age_class",
               values_to = "ecto_burden")

# merge the age class burdens into a list per parameter set 
sens_metrics <- sens_df_long %>% 
  group_by(param_id, theta, iota_baseline, b0, initial_pop,
           iota_r, groom_prg_baseline, groom_stm_baseline, scenario,
           groom_resource_link, age_class) %>% 
  summarise(n_sims = n_distinct(unique_id), 
            burden_vals = list(ecto_burden))

# get the age class lists into the same row 
sens_metrics_wide <- sens_metrics %>%
  group_by(param_id, theta, iota_baseline, b0, initial_pop,
           iota_r, groom_prg_baseline, groom_stm_baseline, scenario,
           groom_resource_link) %>%
  summarise(
    juvenile_vals = list(unlist(burden_vals[age_class == "day_180"])),
    adult_vals    = list(unlist(burden_vals[age_class == "day_780"])),
    .groups = "drop")

# sample into different age class mixes 
param_age_summary <- sens_metrics_wide %>%
  mutate(
    adult_only = map(adult_vals, ~ sample(.x, 100)),
    juvenile_only = map(juvenile_vals, ~ sample(.x, 100)),
    
    equal_mix = map2(juvenile_vals, adult_vals,
                      ~ c(sample(.x, 50),
                          sample(.y, 50))),
    
    adult30_juv70 = map2(juvenile_vals, adult_vals,
                       ~ c(sample(.x, 70),
                           sample(.y, 30))),
    
    adult70_juv30 = map2(juvenile_vals, adult_vals,
                          ~ c(sample(.x, 30),
                              sample(.y, 70)))
  )
# calculate aggregation measures 
param_age_summary <- param_age_summary  %>%
  mutate(
    mean_adult100 = map_dbl(adult_only, mean),
    mean_juv100 = map_dbl(juvenile_only, mean),
    mean_adult50juv50 = map_dbl(equal_mix, mean),
    mean_adult30juv70 = map_dbl(adult30_juv70, mean),
    mean_adult70juv30 = map_dbl(adult70_juv30, mean),
    var_adult100 = map_dbl(adult_only, var),
    var_juv100 = map_dbl(juvenile_only, var),
    var_adult50juv50 = map_dbl(equal_mix, var),
    var_adult30juv70 = map_dbl(adult30_juv70, var),
    var_adult70juv30 = map_dbl(adult70_juv30, var)
  )

# functions to calculate aggregation measures 
calculate_k <- function(mu, v){
  mu^2 / (v - mu)
}

calculate_vmr <- function(mu, v){
  v/mu
}

# calculate the aggregation meaures 
param_age_summary <- param_age_summary %>%
  mutate(
    vmr_adult100      = calculate_vmr(mean_adult100, var_adult100),
    vmr_juv100        = calculate_vmr(mean_juv100, var_juv100),
    vmr_adult50juv50  = calculate_vmr(mean_adult50juv50, var_adult50juv50),
    vmr_adult30juv70  = calculate_vmr(mean_adult30juv70, var_adult30juv70),
    vmr_adult70juv30  = calculate_vmr(mean_adult70juv30, var_adult70juv30),
    
    k_adult100       = calculate_k(mean_adult100, var_adult100),
    k_juv100         = calculate_k(mean_juv100, var_juv100),
    k_adult50juv50   = calculate_k(mean_adult50juv50, var_adult50juv50),
    k_adult30juv70   = calculate_k(mean_adult30juv70, var_adult30juv70),
    k_adult70juv30   = calculate_k(mean_adult70juv30, var_adult70juv30)
  )



# Transform into long: 

param_age_long_num <- param_age_summary %>%
  pivot_longer(
    cols = matches("^(mean|var|k|vmr)_"),
    names_to = c(".value", "age_class"),
    names_pattern = "^(mean|var|k|vmr|fill|class)_(.*)$"
  )

# classify into dispersion type buckets
param_age_long_num <- param_age_long_num %>%
  mutate(
    fill = case_when(
      is.na(k) ~ -2,
      k < 0    ~ -1,
      k > 10   ~ 11,
      TRUE     ~ k
    ),
    class = case_when(
      is.na(k) ~ "Extinct",
      k < 0    ~ "Underdispersed",
      k > 10   ~ "Random",
      TRUE     ~ "Aggregated"
    )
  )

# define levels of adult  
param_age_long_num$age_class <- factor(param_age_long_num$age_class, 
                                       levels = c("adult100",
                                                  "adult70juv30",
                                                  "adult50juv50",
                                                  "adult30juv70",
                                                  "juv100"))


# Plot heatmap of adult samples

k_legend_df <- data.frame(
  class = c("Extinct", "Underdispersed", "Aggregated", "Random"),
  x = NA_real_,
  y = NA_real_
)

k_it <- ggplot(param_age_long_num, aes(x = theta, y = iota_baseline)) +
  geom_tile(aes(fill = fill)) +
  
  scale_fill_gradientn(
    colours = c(
      "yellow",  # extinction
      "blue",    # <0
      "black",   # 0
      "gray85",   # 10
      "red"      # >10
    ),
    values = rescale(c(
      -2,   # extinct
      -1,   # underdispersed
      0,    # start of gradient
      10,   # end of gradient
      11    # random
    )),
    limits = c(-2, 11),
    breaks = c(0, 2, 5, 10),
    labels = c("0", "2", "5", "10"),
    name = "k (0–10)"
  ) +
  
  # dummy legend for class
  geom_point(
    data = k_legend_df,
    aes(x = x, y = y, color = class),
    inherit.aes = FALSE,
    alpha = 0
  ) +
  
  scale_color_manual(
    values = c(
      "Extinct" = "yellow",
      "Underdispersed" = "blue",
      "Aggregated" = "gray50",
      "Random" = "red"
    ),
    name = "Class"
  ) +
  
  facet_grid(age_class ~ interaction(groom_prg_baseline, groom_stm_baseline),
             labeller = labeller(age_class = c(
               adult100 = "Adult (100%)",
               juv100 = "Juvenile (100%)",
               adult50juv50 = "Adult 50% / Juvenile 50%",
               adult70juv30 = "Adult 70% / Juvenile 30%",
               adult30juv70 = "Adult 30% / Juvenile 70%"),
               `interaction(groom_prg_baseline, groom_stm_baseline)` = c(
                 "0.0"   = "No grooming",
                 "20.0"  = "Programmed",
                 "0.20"  = "Stimulus-driven",
                 "20.20" = "Combined"
               ))) +
  
  labs(x = expression(paste("Allee effect (", theta, ")")), 
       y = expression(paste("Baseline immigration (", iota[b0], ")"))) +
  theme_minimal() +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1, size = 5)),
    fill  = guide_colorbar(order = 2)
  )
k_it # confirm all good


# Save Heatmap
png(filename="figures/Figure4_k_sens_412.png", width = 6, height =8, unit = 'in', res = 300)
k_it
dev.off()

# Plot VMR for supplement
vmr_it <- ggplot(param_age_long_num, 
                 aes(x = theta, y = iota_baseline)) +
 geom_tile(aes(fill = log(vmr))) +
             facet_grid(age_class ~ interaction(groom_prg_baseline, groom_stm_baseline),
                        labeller = labeller(age_class = c(
                          adult100 = "Adult (100%)",
                          juv100 = "Juvenile (100%)",
                          adult50juv50 = "Adult 50% / Juvenile 50%",
                          adult70juv30 = "Adult 70% / Juvenile 30%",
                          adult30juv70 = "Adult 30% / Juvenile 70%"),
                          `interaction(groom_prg_baseline, groom_stm_baseline)` = c(
                            "0.0"   = "No grooming",
                            "20.0"  = "Programmed",
                            "0.20"  = "Stimulus-driven",
                            "20.20" = "Combined"
                          ))) +

  labs(x = expression(paste("Allee effect (", theta, ")")), 
       y = expression(paste("Baseline immigration (", iota[b0], ")"))) +
  theme_minimal() +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1, size = 5)),
    fill  = guide_colorbar(order = 2)
  )
vmr_it # confirm all good


# Save Heatmap
png(filename="figures/FigureS4_vmr_sens_412.png", width = 6, height =8, unit = 'in', res = 300)
vmr_it
dev.off()



k_gs_pct_plot <- ggplot(k_gs_pct, aes(x = 1, y = prop, fill = class)) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_grid(~ groom_strategy) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  scale_fill_manual(
    values = c(
      "Extinct" = "yellow",
      "Underdispersed" = "blue",
      "Aggregated" = "gray50",
      "Random" = "red"
    ),
    drop = FALSE,
    name = "Dispersion"
  ) +
  labs(x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )
k_gs_pct_plot # check all good

# Save proportion of k categories
png(filename="figures/FigureX_kpop_cats.png", width = 8, height =8, unit = 'in', res = 300)
k_gs_pct_plot
dev.off()



