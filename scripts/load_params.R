# Load parameters from Excel file
# params <- readxl::read_xlsx(here("data","Flu_params.xlsx"))
# 
# # Filter serial interval
# serial_interval_params <- params %>%
#   filter(epi == "serial_interval") 

# Define time-varying sensitivity function
sensitivity_function <- function(t) {
  if (t < 0) return(0)
  if (t <= 2) return(0.475 * t)
  if (t <= 7) return(0.95 - 0.19 * (t - 2))
  return(0)
}

#I. Ogi-Gittins, W.S. Hart, J. Song, R.K. Nash, J. Polonsky, A. Cori, E.M. Hill, R.N. Thompson. A simulation-based approach for estimating the time-dependent reproduction number from temporally aggregated disease incidence time series data, Epidemics,Volume 47, 2024
si_mean <- 2.6
si_sd <- 1.3

#convert to shape and scale
si_shape <- si_mean^2 / si_sd^2
si_scale <- si_sd^2 / si_mean

#Plot the serial interval distribution
#si_dist <- rgamma(10000, shape = si_shape, scale = si_scale)
#hist(si_dist, breaks = 50, main = "Serial interval distribution", xlab = "Days", freq = FALSE)

#Aditama et al. 2012 Non-zoonotic serial interval distribution (lognormal)
si_draws <- readRDS(here("posterior_predictive","dt_draws.rds")) %>% 
  filter(`Distribution Type` == "Lognormal",
         `Exposure Type` == "Non-Zoonotic",
         `Case Source` == "Serial")

#Plot with ggplot2
si_draws %>%
  ggplot(aes(x = y_rep, fill = `Distribution Type`)) +
  geom_density(alpha = 0.5) +
  labs(title = "Non-zoonotic serial interval distribution",
       x = "Days",
       y = "Density") +
  theme_minimal()
