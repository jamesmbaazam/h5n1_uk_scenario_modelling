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

#---  LOO-CV analysis for picking between Gamma and Lognormal
# Create a named list to store the LOO estimates, extracting them from the nested structure
loo_results <- list(
  gamma_index_non_zoonotic = fit_results$gamma$index$nonzoonotic$loo(),
  gamma_index_zoonotic = fit_results$gamma$index$zoonotic$loo(),
  gamma_serial_non_zoonotic = fit_results$gamma$serial$nonzoonotic$loo(),
  gamma_serial_zoonotic = fit_results$gamma$serial$zoonotic$loo(),
  lognormal_index_non_zoonotic = fit_results$lognormal$index$nonzoonotic$loo(),
  lognormal_index_zoonotic = fit_results$lognormal$index$zoonotic$loo(),
  lognormal_serial_non_zoonotic = fit_results$lognormal$serial$nonzoonotic$loo(),
  lognormal_serial_zoonotic = fit_results$lognormal$serial$zoonotic$loo()
)

print(loo_compare(
  loo_results$gamma_index_non_zoonotic,
  loo_results$lognormal_index_non_zoonotic,
  loo_results$gamma_serial_non_zoonotic,
  loo_results$lognormal_serial_non_zoonotic
), simplify = FALSE)

print(loo_compare(
  loo_results$gamma_index_zoonotic,
  loo_results$lognormal_index_zoonotic,
  loo_results$gamma_serial_zoonotic,
  loo_results$lognormal_serial_zoonotic
), simplify = FALSE)
