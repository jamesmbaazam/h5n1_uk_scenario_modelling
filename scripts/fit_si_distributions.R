library(data.table)
library(cmdstanr)
library(tidybayes)

#--- Source the helper functions
source("R/utils.R")

#--- Loading in data
data_files <- list(
  "index" = "data/si_raw_data/index.csv",
  "serial" = "data/si_raw_data/serial.csv"
)

data_tables <- lapply(data_files, fread)

#-- Creating vectors of intervals
intervals <- lapply(data_tables, create_intervals)

#--- Putting data into format for Stan
stan_data <- lapply(intervals, create_stan_data)

#--- Compiling Stan models
models <- list(
  gamma = cmdstan_model(stan_file = "stan/gamma.stan"),
  lognormal = cmdstan_model(stan_file = "stan/lognormal.stan")
)

#--- Fitting models
fit_results <- list(
  gamma = lapply(stan_data, fit_models, model = models$gamma),
  lognormal = lapply(stan_data, fit_models, model = models$lognormal)
)

#--- Extracting posterior predictive draws using tidybayes
posterior_draws <- list(
  gamma = lapply(fit_results$gamma, function(panel_fits) lapply(panel_fits, extract_draws)),
  lognormal = lapply(fit_results$lognormal, function(panel_fits) lapply(panel_fits, extract_draws))
)

#--- Combining posterior predictive draws
dt_draws <- rbind(
  combine_draws("index", posterior_draws, "gamma"),
  combine_draws("serial", posterior_draws, "gamma"),
  combine_draws("index", posterior_draws, "lognormal"),
  combine_draws("serial", posterior_draws, "lognormal")
)

#--- Saving combined draws
saveRDS(dt_draws, "posterior_predictive/dt_draws.rds")
