library(data.table)
library(cmdstanr)
library(ggplot2)
library(loo)
library(tidybayes)

# Source the helper functions
source("R/utils.R")

# Loading in data
data_files <- list(
  "index" = "data/si_raw_data/index.csv",
  "serial" = "data/si_raw_data/serial.csv"
)

data_tables <- lapply(data_files, fread)

# Combining data.tables for panels A and B, as per new interpretation
dt_onsets <- merge(
  data_tables$index[, .(days, onsets_index = nonzoonotic)],
  data_tables$serial[, .(days, onsets_serial = nonzoonotic)])[
  , onsets := onsets_index + onsets_serial][, .(days, onsets)]

# Calculate the intervals for fitting the serial interval distributions where 
# Interpret the “onsets” as frequency counts. I.e., the raw data provides counts
# of the number of cases that showed symptoms on specific days. 
# We convert this into a format representing each interval explicitly
dt_onsets_intervals <- data.table(
  time =  0:length(dt_onsets$days),
  onsets = rep(dt_onsets$days, dt_onsets$onsets))

# Compile Stan models
mod_lognormal <- cmdstan_model("stan/distributions/lognormal.stan")
mod_gamma <- cmdstan_model("stan/distributions/gamma.stan")

# Get data into format for Stan
stan_data <- list(
  N = dt_onsets_intervals[, .N],
  y = dt_onsets_intervals[, onsets],
  N_rep = 100)

# Fit models to data
fit_lognormal <- mod_lognormal$sample(stan_data, chains = 4, parallel_chains = 4)
fit_gamma <- mod_gamma$sample(stan_data, chains = 4, parallel_chains = 4)

# Extract posterior samples
dt_lognormal <- melt(
  data.table(fit_lognormal$draws("y_rep", format = "data.frame")),
  measure.vars = patterns("y_rep"))

dt_gamma <- melt(
  data.table(fit_gamma$draws("y_rep", format = "data.frame")),
  measure.vars = patterns("y_rep"))

dt_lognormal[, variable := NULL]
dt_gamma[, variable := NULL]

# Combine distribution type
dt_si_posteriors <- rbind(
  dt_lognormal[, type := "Lognormal"],
  dt_gamma[, type := "Gamma"])

# Summarise distributions
dt_summary_median <- summarise_draws(
  dt_si_posteriors, by = "type", 
  column_name = "value")[, average := "Median"]

# Summarise distributions
dt_summary_mean <- dt_si_posteriors[, .(
  me = mean(value), 
  lo = mean(value) - IQR(value),
  hi = mean(value) + IQR(value)), 
  by = "type"][, average := "Mean"]

# Creating PMF of underlying data for plotting behind fitted distributions
dt_onsets_plot <- rbind(
  data.table(time = dt_onsets_intervals$onsets, type = "Lognormal"),
  data.table(time = dt_onsets_intervals$onsets, type = "Gamma")
)

# Create the plot with observed data
p_si <- ggplot() +
  geom_histogram(
    data = dt_onsets_plot,
    aes(x = time, y = after_stat(density)),
    binwidth = 1,
    fill = "gray",
    alpha = 1,
    boundary = 0,
    closed = "left"
  ) +
  geom_density(
    data = dt_si_posteriors,
    aes(x = value, fill = type),
    alpha = 0.8
  ) +
  geom_vline(
    data = dt_summary_median,
    aes(xintercept = me),
    linetype = "dashed"
  ) +
  facet_grid(~type) +
  theme_linedraw() +
  theme(legend.position = "none") +
  labs(x = "Time (days since exposure)", y = "Density") +
  lims(x = c(0, 30))

ggsave("plots/serial_interval.png", p_si, width = 8, height = 4)

knitr::kable(
  rbind(dt_summary_median, dt_summary_mean)[order(type, average)])


