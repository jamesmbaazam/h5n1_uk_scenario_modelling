library(data.table)
library(ggplot2)
library(loo)
library(tidybayes)
library(forcats)

dt_draws <- readRDS("posterior_predictive/dt_draws.rds")

#--- Adding mean, median and 95% credible intervals to data.table
dt_draws[, `:=` (
  Mean = mean(y_rep),
  Median = quantile(y_rep, 0.5),
  `Lower CrI` = quantile(y_rep, 0.025),
  `Upper CrI` = quantile(y_rep, 0.975)),
  by = .(`Exposure Type`, `Case Source`, `Distribution Type`)][
  , `Case Source` := fct_relevel(`Case Source`, c("Serial", "Index"))]

#--- Printing summary table
dt_draws[order(`Exposure Type`, `Case Source`, `Distribution Type`)
  , .(Mean, Median, `Lower CrI`, `Upper CrI`),
  by = .(`Exposure Type`, `Case Source`, `Distribution Type`)] |>
  unique()

#--- Plotting
dt_draws |> 
  ggplot() + 
  geom_density(
    aes(x = y_rep,
        fill = interaction(`Exposure Type`, `Distribution Type`, `Case Source`)),
    alpha = 0.5) +
  geom_vline(aes(
    xintercept = Median, 
    colour = interaction(`Exposure Type`, `Distribution Type`, `Case Source`)),
    linetype = "dashed") +
  facet_grid(`Distribution Type` ~ `Exposure Type` + `Case Source`) +
  lims(x = c(0, 30)) + 
  labs(
    x = "Time (days)", y = "Probability density") + 
  theme_minimal() +
  theme(legend.position = "none")
