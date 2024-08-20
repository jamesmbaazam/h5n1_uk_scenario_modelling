library(data.table)
library(cmdstanr)
library(ggplot2)
library(tidybayes)

dt_panel_a <- fread("data/si_raw_data/panel_a.csv")
dt_panel_b <- fread("data/si_raw_data/panel_b.csv")

# Create vectors of intervals for both groups
panel_a_intervals_white <- rep(dt_panel_a$days, dt_panel_a$cases_white)
panel_a_intervals_black <- rep(dt_panel_a$days, dt_panel_a$cases_black)
panel_b_intervals_white <- rep(dt_panel_b$days, dt_panel_b$cases_white)
panel_b_intervals_black <- rep(dt_panel_b$days, dt_panel_b$cases_black)

stan_data_panel_a_white <- list(
  N = length(panel_a_intervals_white),
  y = panel_a_intervals_white
)

stan_data_panel_a_black <- list(
  N = length(panel_a_intervals_black),
  y = panel_a_intervals_black
)

stan_data_panel_b_white <- list(
  N = length(panel_b_intervals_white),
  y = panel_b_intervals_white
)

stan_data_panel_b_black <- list(
  N = length(panel_b_intervals_black),
  y = panel_b_intervals_black
)

# For the Gamma distribution
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

fit_panel_b_black_lognormal<- mod_lognormal$sample(
  data = stan_data_panel_b_black, chains = 4, 
  iter_warmup = 1000, iter_sampling = 2000)

#--- Extracting posterior predictive draws
res_panel_a_white_gamma <- tidybayes::spread_draws(
  fit_panel_a_white_gamma, y_rep[i]) |> data.table()

res_panel_a_black_gamma <- tidybayes::spread_draws(
  fit_panel_a_black_gamma, y_rep[i]) |> data.table()

res_panel_b_white_gamma <- tidybayes::spread_draws(
  fit_panel_b_white_gamma, y_rep[i]) |> data.table()

res_panel_b_black_gamma <- tidybayes::spread_draws(
  fit_panel_b_black_gamma, y_rep[i]) |> data.table()

res_panel_a_white_lognormal <- tidybayes::spread_draws(
  fit_panel_a_white_lognormal, y_rep[i]) |> data.table()

res_panel_a_black_lognormal <- tidybayes::spread_draws(
  fit_panel_a_black_lognormal, y_rep[i]) |> data.table()

res_panel_b_white_lognormal <- tidybayes::spread_draws(
  fit_panel_b_white_lognormal, y_rep[i]) |> data.table()

res_panel_b_black_lognormal <- tidybayes::spread_draws(
  fit_panel_b_black_lognormal, y_rep[i]) |> data.table()

dt_draws_gamma <- 
  rbind(res_panel_a_white_gamma[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Index to serial (Panel A)"],
        res_panel_a_black_gamma[, `Exposure type` := "Zoonotic"][, `Case type` := "Index to serial (Panel A)"],
        res_panel_b_white_gamma[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Serial (Panel B)"],
        res_panel_b_black_gamma[, `Exposure type` := "Zoonotic"][, `Case type` := "Serial (Panel B)"])[
          , panel := "A"][, Distribution := "Gamma"]

dt_draws_lognormal <- 
  rbind(res_panel_a_white_lognormal[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Index to serial (Panel A)"],
        res_panel_a_black_lognormal[, `Exposure type` := "Zoonotic"][, `Case type` := "Index to serial (Panel A)"],
        res_panel_b_white_lognormal[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Serial (Panel B)"],
        res_panel_b_black_lognormal[, `Exposure type` := "Zoonotic"][, `Case type` := "Serial (Panel B)"])[
          , panel := "A"][, Distribution := "Lognormal"]

dt_draws <- rbind(dt_draws_gamma, dt_draws_lognormal)

dt_draws |> 
  ggplot() + 
  geom_density(aes(x = y_rep, fill = `Distribution`), alpha = 0.2) +
  facet_grid(`Case type` ~ `Exposure type`) + 
  lims(x = c(0, 30)) + 
  labs(x = "Time (days)", y = "Probability density") + 
  theme_minimal()

