library(data.table)
library(cmdstanr)

#--- Loading in data
dt_panel_a <- fread("data/si_raw_data/panel_a.csv")
dt_panel_b <- fread("data/si_raw_data/panel_b.csv")

#-- Creating vectors of intervals for both groups
panel_a_intervals_white <- rep(dt_panel_a$days, dt_panel_a$cases_white)
panel_a_intervals_black <- rep(dt_panel_a$days, dt_panel_a$cases_black)
panel_b_intervals_white <- rep(dt_panel_b$days, dt_panel_b$cases_white)
panel_b_intervals_black <- rep(dt_panel_b$days, dt_panel_b$cases_black)

#--- Putting data into format for Stan
stan_data_panel_a_white <- list(
  N = length(panel_a_intervals_white),
  N_rep = 1000,
  y = panel_a_intervals_white
)

stan_data_panel_a_black <- list(
  N = length(panel_a_intervals_black),
  N_rep = 1000,
  y = panel_a_intervals_black
)

stan_data_panel_b_white <- list(
  N = length(panel_b_intervals_white),
  N_rep = 1000,
  y = panel_b_intervals_white
)

stan_data_panel_b_black <- list(
  N = length(panel_b_intervals_black),
  N_rep = 1000,
  y = panel_b_intervals_black
)

#--- Compiling Stan models
mod_gamma <- cmdstan_model(stan_file = "stan/gamma.stan")
mod_lognormal <- cmdstan_model(stan_file = "stan/lognormal.stan")

#--- Fitting Gamma distributions
fit_panel_a_white_gamma <- mod_gamma$sample(
  data = stan_data_panel_a_white, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

fit_panel_a_black_gamma <- mod_gamma$sample(
  data = stan_data_panel_a_black, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

fit_panel_b_white_gamma <- mod_gamma$sample(
  data = stan_data_panel_b_white, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

fit_panel_b_black_gamma <- mod_gamma$sample(
  data = stan_data_panel_b_black, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

#--- Fitting Lognormal distributions
fit_panel_a_white_lognormal <- mod_lognormal$sample(
  data = stan_data_panel_a_white, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

fit_panel_a_black_lognormal <- mod_lognormal$sample(
  data = stan_data_panel_a_black, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

fit_panel_b_white_lognormal <- mod_lognormal$sample(
  data = stan_data_panel_b_white, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

fit_panel_b_black_lognormal <- mod_lognormal$sample(
  data = stan_data_panel_b_black, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

#--- Saving fits
saveRDS(fit_panel_a_white_gamma, "fits/fit_panel_a_white_gamma.rds")
saveRDS(fit_panel_a_black_gamma, "fits/fit_panel_a_black_gamma.rds")
saveRDS(fit_panel_b_white_gamma, "fits/fit_panel_b_white_gamma.rds")
saveRDS(fit_panel_b_black_gamma, "fits/fit_panel_b_black_gamma.rds")
saveRDS(fit_panel_a_white_lognormal, "fits/fit_panel_a_white_lognormal.rds")
saveRDS(fit_panel_a_black_lognormal, "fits/fit_panel_a_black_lognormal.rds")
saveRDS(fit_panel_b_white_lognormal, "fits/fit_panel_b_white_lognormal.rds")
saveRDS(fit_panel_b_black_lognormal, "fits/fit_panel_b_black_lognormal.rds")



