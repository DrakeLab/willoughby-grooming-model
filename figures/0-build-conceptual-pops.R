# Plot dispersion scenarios 
library(ggplot2) # plotting 
set.seed(1993)

# define the host and parasite pops 
n_hosts <- 100
total_parasites <- 500

# Random
rand <- rpois(n_hosts, lambda = 4)
rand <- round(rand * total_parasites / sum(rand))

# Moderate
und <- rep(5, n_hosts)   # perfectly even

# Overdispersed 
ove <- rnbinom(n_hosts, size = 0.2, mu = 5)
ove <- round(ove * total_parasites / sum(ove))

# pop comparison data frame 
hp_df <- data.frame(
  parasites = c(rand, und, ove),
  dispersion = rep(c("b. Randomly dispersed",
                     "a. Underdispersed",
                     "c. Overdispersed"), each = n_hosts)
  )

# define levels for random on the left 
hp_df$dispersion <- factor(hp_df$dispersion, levels = c("a. Underdispersed", 
                                                        "b. Randomly dispersed", 
                                                        "c. Overdispersed"))

calc_metrics <- function(parasites) {
  
  prev <- round(sum(parasites > 0) / length(parasites)*100,2)
  
  h_infected <- sum(parasites > 0)
  inte <- ifelse(h_infected > 0, 
                       sum(parasites) / h_infected ,NA_real_)  # mean intensity
  abun <- sum(parasites) / length(parasites)
  
  mu <- mean(parasites)
  v  <- var(parasites)
  
  # negative binomial aggregation parameter
  k <- if (v > mu) mu^2 / (v - mu) else NA_real_
  
  # Lloyd's mean crowding
  crowd <- mu + (v / mu) - 1 
  
  # patchiness
  patch <- crowd / mu 
  
  data.frame(
    prev = prev,
    abun = abun,
    inte = inte,
    mu = mu,
    var = v,
    k = k,
    crowd = crowd,
    patch = patch
  )
}

# Calculate metrics 
metrics <- do.call(
  rbind,
  lapply(split(hp_df$parasites, hp_df$dispersion), calc_metrics)
)

metrics$dispersion <- c("a. Underdispersed", "b. Randomly dispersed", "c. Overdispersed")
metrics$dispersion <- factor(metrics$dispersion, levels = c("a. Underdispersed", "b. Randomly dispersed", "c. Overdispersed"))


metrics$label <- sprintf(
  "Prevalance = %.2f%%\nMean abundance = %.2f\nMean intensity = %.2f\nAggregation (k) = %.2f",
  metrics$prev,
  metrics$abun,
  metrics$inte,
  metrics$k
)


# plot 

concept_pops_plot <- ggplot(hp_df, aes(x = parasites, fill=parasites==0, )) +
    geom_histogram(binwidth = 1, color="black", boundary = -0.5) +
    facet_wrap(~dispersion, scales="free_y") +
    coord_cartesian(xlim = c(0, max(hp_df$parasites, na.rm = TRUE))) + 
    scale_fill_manual(values = c("TRUE" = "gray70", "FALSE" = "red"),guide = "none") + 
    geom_label(
      data = metrics,
      aes(label = label),
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.2,
      size = 3,
      inherit.aes = FALSE) + 
    labs(x=expression(bold("Parasites per host")), 
         y=expression(bold("Frequency"))) + 
  theme_classic() +
  theme(strip.text = element_text(face = "bold"))  
 
concept_pops_plot # check all good 

png(file = "data/figures/Fig1_ConceptPops.png")
dev.off()

# Loop over all possible k values 
## Target k values: 
k_vals <- c(20, 15, 10, 5, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01) 

# Define the Host and Parasite Population Sizes 
n_hosts <- 100
total_parasites <- 500

sim_list <- lapply(k_vals, function(k_target) {
  
  parasites <- rnbinom(n_hosts, size = k_target, mu = 5)
  parasites <- round(parasites * total_parasites / sum(parasites))
  
  metrics <- calc_metrics(parasites)
  metrics$k_target <- k_target
  metrics$sum_parasites <- sum(parasites)
  
  list(
    k_target = k_target,
    parasites = parasites,
    metrics = metrics
  )
})

# build the data frame 
sim_metric_df <- do.call(rbind, lapply(sim_list, function(x) x$metrics))
sim_metric_df

# define the labels 
sim_metric_df$label <- sprintf(
  "Prevalance = %.2f%%\nMean abundance = %.2f\nMean intensity = %.2f\nAggregation (k) = %.2f",
  sim_metric_df$prev,
  sim_metric_df$abun,
  sim_metric_df$inte,
  sim_metric_df$k
)


# Plotting data frame 
plot_one_k_sim <- function(sim, metric_df, k_vals = NULL, file_name) {
  
  plot_df <- data.frame(
    parasites = sim$parasites,
    k_target = factor(sim$k_target, levels = k_vals)
  )
  
  label_df <- metric_df %>%
    dplyr::filter(k_target == sim$k_target)
  
  p <- ggplot(plot_df, aes(x = parasites, fill = parasites == 0)) +
    geom_histogram(binwidth = 1, boundary = -0.5, color = "black") +
    coord_cartesian(xlim = c(0, max(plot_df$parasites, na.rm = TRUE))) +
    scale_fill_manual(
      values = c("TRUE" = "gray70", "FALSE" = "red"),
      guide = "none"
    ) +
    geom_label(
      data = label_df,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.05,
      vjust = 1.2,
      size = 3,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Parasites per host",
      y = "Frequency"
    ) +
    theme_classic()
  
  png(file = file_name, width = 4, height = 4, unit = 'in', res = 300)
  print(p)
  dev.off()
  
  return(p)
}
# Loop through all target ks 
# invisible(lapply(sim_list, plot_one_k_sim, metric_df = sim_metric_df, k_vals = k_vals))

for (j in seq_along(sim_list)) {
  sim <- sim_list[[j]]
  plot_one_k_sim(
    sim = sim,
    metric_df = sim_metric_df,
    k_vals = k_vals,
    file_name = paste("figures/k/frame_", sprintf("%02d", j), ".png", sep = "")
  )
}

# Full facetted plot 
sim_k_plot <- ggplot(hp_df, aes(x = parasites, fill=parasites==0, )) +
  geom_histogram(binwidth = 1, color="black", boundary = -0.5) +
  facet_wrap(~dispersion, scales="free_y") +
  coord_cartesian(xlim = c(0, max(hp_df$parasites, na.rm = TRUE))) + 
  scale_fill_manual(values = c("TRUE" = "gray70", "FALSE" = "red"),guide = "none") + 
  geom_label(
    data = metrics,
    aes(label = label),
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.2,
    size = 3,
    inherit.aes = FALSE) + 
  labs(x=expression(bold("Parasites per host")), 
       y=expression(bold("Frequency"))) + 
  theme_classic() +
  theme(strip.text = element_text(face = "bold"))  

sim_k_plot # check all good 

# Build video of all frames --> Go To Terminal for using ffmeg
# ffmpeg -framerate 1 -i figures/k/frame_%02d.png -c:v libx264 -pix_fmt yuv420p figures/k/k_animation.mp4

#
# install.packages("magick")
library(magick)

# Load in image list
imgs <- image_read(list.files("figures/k", pattern = "^frame_\\d+\\.png$", full.names = TRUE))
gif <- image_animate(imgs, fps = 1) # animate
image_write(gif, "figures/k/k_animation.gif") # save 


# Plot out relationships of metrics with k 
## PREVALENCE
plot(prev ~ k , data = sim_metric_df) 

### look at a log 
plot(prev ~ log(k) , data = sim_metric_df) 

##  Ratio of Mean ABUN to MEAN INTENSITY
plot(abun/inte ~ k , data = sim_metric_df) 

##  Ratio of Mean ABUN to MEAN INTENSITY
plot(abun/inte ~ log(k) , data = sim_metric_df) 




