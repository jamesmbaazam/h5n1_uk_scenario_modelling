# Load parameters from Excel file
# params <- readxl::read_xlsx(here("data","Flu_params.xlsx"))
# 
# # Filter serial interval
# serial_interval_params <- params %>%
#   filter(epi == "serial_interval") 

# Allow selection of flu type
FLU_TYPE <- "H5N1 HPAI"  # Options: "H5N1 HPAI" or "H1N1"

# Get SI parameters based on flu type
if(FLU_TYPE == "H5N1 HPAI") {
  # Load H5N1 HPAI posterior draws
  si_posteriors <- fread("data/si_posteriors.csv")
  si_draws_data <- si_posteriors[type == "Gamma", value]
  
  # Check for NAs
  if(any(is.na(si_draws_data))) {
    stop("Found NA values in H5N1 HPAI SI data")
  }
  
  # Print summary for debugging
  message("H5N1 HPAI SI summary:")
  print(summary(si_draws_data))
  
  # Function to generate draws
  si_draws <- function(n) {
    sample(si_draws_data, size = n, replace = TRUE)
  }
} else if(FLU_TYPE == "H1N1") {
  # Get H1N1 SI from epiparameter database
  flu_si <- epiparameter::epiparameter_db(
    disease = "influenza",
    epi_name = "serial interval",
    pathogen = "H1N1",
    single_epiparameter = TRUE
  )
  
  # Get parameters and check for NAs
  si_params <- get_parameters(flu_si)
  if(any(is.na(si_params))) {
    stop("Retrieved NA parameters for H1N1 SI: ", 
         paste(names(si_params[is.na(si_params)]), collapse=", "))
  }
  
  # Print parameters for debugging
  message("H1N1 SI parameters:")
  print(si_params)
  
  # Function to generate draws
  si_draws <- function(n) {
    draws <- rlnorm(n, 
                    meanlog = si_params["meanlog"], 
                    sdlog = si_params["sdlog"])
    if(any(is.na(draws))) {
      stop("Generated NA values in SI draws")
    }
    draws
  }
} else {
  stop("Invalid FLU_TYPE: must be 'H5N1 HPAI' or 'H1N1'")
}

# Generate test draws and check for errors
tryCatch({
  test_draws <- si_draws(10)
  message("Test draws successful:")
  print(summary(test_draws))
}, error = function(e) {
  stop("Error generating test SI draws: ", e$message)
})

# Load PCR sensitivity curve data
pcr_data <- fread("data/PCR_curve_summary.csv")

#COVID SI for scaling
# Nishiura H, Linton N, Akhmetzhanov A (2020). "Serial interval of novel coronavirus (COVID-19) infections." _International Journal of Infectious Diseases_. doi:10.1016/j.ijid.2020.02.060https://doi.org/10.1016/j.ijid.2020.02.060
covid_si <- epiparameter::epiparameter_db(
  disease = "COVID-19",
  epi_name = "serial interval",
  single_epiparameter = TRUE
)

covid_si_params <- get_parameters(covid_si)

# Get median serial intervals
flu_si_median <- median(si_draws(1000))
covid_si_median <- exp(covid_si_params["meanlog"]) # Convert from log scale

# Scale factor for time axis (flu SI / covid SI)
time_scale_factor <- flu_si_median / covid_si_median

# Scale the PCR data time axis
pcr_data[, days_since_infection := days_since_infection * time_scale_factor]

# Define time-varying sensitivity function using linear interpolation
# with optional uncertainty sampling
sensitivity_function <- function(t, sample_uncertainty = FALSE) {
  # Vectorized bounds check
  out_of_bounds <- t < min(pcr_data$days_since_infection) | t > max(pcr_data$days_since_infection)
  
  if (!sample_uncertainty) {
    # Just return median interpolation
    result <- approx(pcr_data$days_since_infection, pcr_data$median, xout = t)$y
    result[out_of_bounds] <- 0
    return(result)
  } else {
    # Get interpolated values for median and bounds
    median_t <- approx(pcr_data$days_since_infection, pcr_data$median, xout = t)$y
    lower_t <- approx(pcr_data$days_since_infection, pcr_data$lower_95, xout = t)$y
    upper_t <- approx(pcr_data$days_since_infection, pcr_data$upper_95, xout = t)$y
    
    # Calculate standard deviation assuming normal distribution
    sd_t <- (upper_t - lower_t)/3.92
    
    # Sample from truncated normal distribution
    samples <- rnorm(length(t), median_t, sd_t)
    samples <- pmax(0, pmin(1, samples))  # Truncate to [0,1]
    samples[out_of_bounds] <- 0
    return(samples)
  }
}

# Plot both serial interval distributions for comparison
si_plot <- ggplot() +
  geom_density(data = data.frame(x = si_draws(1000)), 
              aes(x = x, color = paste(tools::toTitleCase(FLU_TYPE), "Flu"))) +
  stat_function(aes(color = "COVID-19"), 
               fun = function(x) dlnorm(x, covid_si_params["meanlog"], covid_si_params["sdlog"])) +
  labs(title = paste(tools::toTitleCase(FLU_TYPE), "Flu vs COVID-19 Serial Interval Distributions"),
       x = "Days",
       y = "Density",
       color = "Disease") +
  theme_minimal()

# Save SI plot
ggsave("plots/serial_intervals_comparison.png", si_plot, 
       width = 8, height = 6, bg = "white", dpi = 600)
ggsave("plots/serial_intervals_comparison.pdf", si_plot, 
       width = 8, height = 6, bg = "white")

# Plot examples of scaled PCR sensitivity curves
t_values <- seq(0, max(pcr_data$days_since_infection), length.out = 100)
n_samples <- 100

# Generate samples
samples <- data.frame(
  t = rep(t_values, n_samples),
  sample_id = rep(1:n_samples, each = length(t_values)),
  sensitivity = sapply(rep(t_values, n_samples), function(t) sensitivity_function(t, TRUE))
)

# Create plot with samples and original curve
sens_plot <- ggplot() +
  geom_line(data = samples, 
            aes(x = t, y = sensitivity, group = sample_id),
            alpha = 0.1, color = "blue") +
  geom_ribbon(data = pcr_data,
             aes(x = days_since_infection, 
                 ymin = lower_95, 
                 ymax = upper_95),
             alpha = 0.2) +
  geom_line(data = pcr_data,
            aes(x = days_since_infection, y = median),
            color = "red") +
  labs(title = paste(tools::toTitleCase(FLU_TYPE), "Flu PCR Test Sensitivity Over Time"),
       subtitle = sprintf("Time axis scaled by %.2fx (flu SI / COVID SI)", time_scale_factor),
       x = "Days Since Infection",
       y = "Sensitivity") +
  theme_minimal()

# Save sensitivity plot
ggsave("plots/pcr_sensitivity_scaled.png", sens_plot, 
       width = 8, height = 6, bg = "white", dpi = 600)
ggsave("plots/pcr_sensitivity_scaled.pdf", sens_plot, 
       width = 8, height = 6, bg = "white")


