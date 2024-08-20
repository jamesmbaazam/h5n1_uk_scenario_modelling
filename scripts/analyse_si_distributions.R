library(ggplot2)
library(tidybayes)
library(loo)

#--- Loading fits
fit_panel_a_white_gamma <- readRDS("fits/fit_panel_a_white_gamma.rds")
fit_panel_a_black_gamma <- readRDS("fits/fit_panel_a_black_gamma.rds")
fit_panel_b_white_gamma <- readRDS("fits/fit_panel_b_white_gamma.rds")
fit_panel_b_black_gamma <- readRDS("fits/fit_panel_b_black_gamma.rds")
fit_panel_a_white_lognormal <- readRDS("fits/fit_panel_a_white_lognormal.rds")
fit_panel_a_black_lognormal <- readRDS("fits/fit_panel_a_black_lognormal.rds")
fit_panel_b_white_lognormal <- readRDS("fits/fit_panel_b_white_lognormal.rds")
fit_panel_b_black_lognormal <- readRDS("fits/fit_panel_b_black_lognormal.rds")

#--- Extracting posterior predictive draws using tidybayes
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

#--- Combining posterior predictive draws
dt_draws_gamma <- 
  rbind(
    res_panel_a_white_gamma[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Index to serial (Panel A)"],
    res_panel_a_black_gamma[, `Exposure type` := "Zoonotic"][, `Case type` := "Index to serial (Panel A)"],
    res_panel_b_white_gamma[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Serial (Panel B)"],
    res_panel_b_black_gamma[, `Exposure type` := "Zoonotic"][, `Case type` := "Serial (Panel B)"])[, Distribution := "Gamma"]

dt_draws_lognormal <- 
  rbind(
    res_panel_a_white_lognormal[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Index to serial (Panel A)"],
    res_panel_a_black_lognormal[, `Exposure type` := "Zoonotic"][, `Case type` := "Index to serial (Panel A)"],
    res_panel_b_white_lognormal[, `Exposure type` := "Non-zoonotic"][, `Case type` := "Serial (Panel B)"],
    res_panel_b_black_lognormal[, `Exposure type` := "Zoonotic"][, `Case type` := "Serial (Panel B)"])[, Distribution := "Lognormal"]

dt_draws <- rbind(dt_draws_gamma, dt_draws_lognormal)

#--- Plotting
dt_draws |> 
  ggplot() + 
  geom_density(aes(x = y_rep, fill = `Distribution`), alpha = 0.2) +
  facet_grid(`Case type` ~ `Exposure type`) + 
  lims(x = c(0, 30)) + 
  labs(x = "Time (days)", y = "Probability density") + 
  theme_minimal()

dt_draws[order(`Exposure type`, `Case type`)][Distribution == "Lognormal"][
  , .(mean = mean(y_rep), median = median(y_rep)), 
  by = .(`Exposure type`, `Case type`)]

#--- LOO-CV analysis for picking between Gamma and Lognormal
# Calculating LOO estimates
loo_panel_a_white_gamma <- fit_panel_a_white_gamma$loo()
loo_panel_a_black_gamma <- fit_panel_a_black_gamma$loo()
loo_panel_a_black_gamma <- fit_panel_b_white_gamma$loo()
loo_panel_b_black_gamma <- fit_panel_b_black_gamma$loo()
loo_panel_a_white_lognormal <- fit_panel_a_white_lognormal$loo()
loo_panel_a_black_lognormal <- fit_panel_a_black_lognormal$loo()
loo_panel_b_white_lognormal <- fit_panel_b_white_lognormal$loo()
loo_panel_b_black_lognormal <- fit_panel_b_black_lognormal$loo()

# Comparing LOO estimates
print(loo_compare(loo_panel_a_white_gamma, loo_panel_a_white_lognormal), simplify = FALSE)
print(loo_compare(loo_panel_a_black_gamma, loo_panel_a_black_lognormal), simplify = FALSE)
print(loo_compare(loo_panel_b_white_gamma, loo_panel_b_white_lognormal), simplify = FALSE)
print(loo_compare(loo_panel_b_black_gamma, loo_panel_b_black_lognormal), simplify = FALSE)

