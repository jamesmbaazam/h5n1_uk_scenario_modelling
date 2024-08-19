# Install and load required packages
if (!require("pak")) install.packages("pak")
pak::pak(c("epiverse-trace/epichains", "epiverse-trace/epiparameter", "tidyverse", "truncdist", "MASS", "fitdistrplus", "ggplot2", "gridExtra"))

# Load required libraries
library(tidyverse)
library(epichains)
library(truncdist)
library(epiparameter)
library(MASS)
library(fitdistrplus)
library(ggplot2)
library(gridExtra)

# Set seed for reproducibility
set.seed(123)

# Define time-varying sensitivity function
sensitivity_function <- function(t) {
  if (t < 0) return(0)
  if (t <= 2) return(0.475 * t)
  if (t <= 7) return(0.95 - 0.19 * (t - 2))
  return(0)
}

# Define the flight_test_fun function
flight_test_fun <- function(chain_data, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay, flight_time_range, flight_duration) {
  
  results <- chain_data %>%
    mutate(
      flight_time = runif(n(), min = flight_time_range[1], max = flight_time_range[2]),
      flight_end = flight_time + flight_duration,
      pre_flight_test_t = Inf,
      pre_flight_test_res = 0,
      post_flight_test_t = Inf,
      post_flight_test_res = 0,
      isolated_pre = FALSE,
      isolated_post = FALSE,
      isolation_time = Inf
    )
  
  # Test before flight
  if (test_before_flight == 1) {
    results <- results %>%
      mutate(
        pre_flight_test_t = flight_time - pre_flight_test_delay,
        pre_flight_test_sensitivity = sapply(pre_flight_test_t, sensitivity_function),
        pre_flight_test_res = rbinom(n(), 1, pre_flight_test_sensitivity),
        isolated_pre = ifelse(pre_flight_test_res == 1, TRUE, isolated_pre)
      )
  }
  
  # Test after flight
  if (test_after_flight == 1) {
    results <- results %>%
      mutate(
        post_flight_test_t = flight_end + post_flight_test_delay,
        post_flight_test_sensitivity = sapply(post_flight_test_t, sensitivity_function),
        post_flight_test_res = rbinom(n(), 1, post_flight_test_sensitivity),
        isolated_post = ifelse(post_flight_test_res == 1, TRUE, isolated_post)
      )
  }
  
  # Earliest isolation time
  results <- results %>%
    mutate(
      isolation_time = pmin(ifelse(isolated_pre, pre_flight_test_t, Inf), ifelse(isolated_post, post_flight_test_t, Inf))
    )
  
  # Indicate whether 2nd gen infections were averted if the infector was isolated
  gen_2_averted <- results %>% 
    filter(generation > 1) %>% 
    dplyr::select(chain, infector, generation, flight_time, time) %>%
    right_join(results %>% filter(generation == 1) %>% 
                 dplyr::select(chain, isolation_time), 
               by = "chain") %>%
    mutate(averted = time >= isolation_time)
  
  return(gen_2_averted)
}

# Function to calculate R and k using negative binomial fitting
calculate_r_and_k <- function(chain_data) {
  # Count number of offspring per chain using dplyr count (don't drop empty)
  
  offspring_counts <- chain_data %>%
    mutate(chain = as.factor(chain),
           generation = as.factor(generation),
           infector = as.factor(infector)) %>%
    filter(!averted, time >= flight_time) %>% 
    group_by(chain, .drop=F) %>%
    summarise(N = n()) %>%
    pull(N)
  
  # Check if all offspring counts are zero
  if (all(offspring_counts == 0)) {
    warning("All offspring counts are zero. Unable to estimate R and k.")
    return(list(R = NA, k = NA))
  }
  
  # Fit negative binomial distribution
  tryCatch({
    nb_fit <- fitdistrplus::fitdist(as.vector(offspring_counts), "nbinom")
    r <- nb_fit$estimate["mu"]
    k <- nb_fit$estimate["size"]
    return(list(R = r, k = k))
  }, error = function(e) {
    warning("Error in fitting negative binomial distribution: ", e$message)
    return(list(R = NA, k = NA))
  })
}

# Function to run a single scenario
run_scenario <- function(sim_chains_df, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay, flight_time_range, flight_duration, scenario_name) {
  filtered_results <- flight_test_fun(
    chain_data = sim_chains_df,
    test_before_flight = test_before_flight,
    test_after_flight = test_after_flight,
    pre_flight_test_delay = pre_flight_test_delay,
    post_flight_test_delay = post_flight_test_delay,
    flight_time_range = flight_time_range,
    flight_duration = flight_duration
  )
  
  # Calculate R and k
  r_k <- calculate_r_and_k(filtered_results)
  
  cat(scenario_name, ":\n")
  cat("Number of exported infections:", nrow(filtered_results), "\n")
  if (!is.na(r_k$R) && !is.na(r_k$k)) {
    cat("Estimated R:", r_k$R, "\n")
    cat("Estimated k:", r_k$k, "\n")
  } else {
    cat("Unable to estimate R and k. There might not be enough data.\n")
  }
  cat("\n")
  
  return(list(filtered_results = filtered_results, 
              R = r_k$R,
              k = r_k$k))
}

# Simulate chains
sim_chains <- simulate_chains(
  n_chains = 1000,  # Increased for more robust R and k estimation
  statistic = "length",
  offspring_dist = rnbinom,
  stat_threshold = 2,
  generation_time = function(n) rdunif(n, 1, 4),
  mu = 2.5,
  size = 0.1
)

# Convert sim_chains to dataframe
sim_chains_df <- as.data.frame(sim_chains)

cat("\nNumber of rows in original sim_chains:", nrow(sim_chains_df), "\n\n")

# Define scenarios
scenarios <- list(
  list(name = "No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_time_range = c(1, 4), flight_duration = 0.2),
  list(name = "Pre-flight (1 day before)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_time_range = c(1, 4), flight_duration = 0.2),
  list(name = "Post-flight (1 day after)", pre = 0, post = 1, pre_delay = NA, post_delay = 1, flight_time_range = c(1, 4), flight_duration = 0.2),
  list(name = "Both (1 day before and after)", pre = 1, post = 1, pre_delay = 1, post_delay = 1, flight_time_range = c(1, 4), flight_duration = 0.2)
)

# Run all scenarios
results <- lapply(scenarios, function(s) {
  run_scenario(sim_chains_df, s$pre, s$post, s$pre_delay, s$post_delay, s$flight_time, s$flight_duration, s$name)
})

# Summarise results
results_summary <- data.frame(
  Scenario = sapply(scenarios, function(s) s$name),
  R = sapply(results, function(x) x$R),
  k = sapply(results, function(x) x$k)
)

# Calculate relative risk
no_intervention_R <- results_summary$R[results_summary$Scenario == "No testing"]
results_summary$RelativeRisk <- results_summary$R / no_intervention_R

cat("Summary of results:\n")
print(results_summary)

# Visualise the results
# Plot R values
p1 <- ggplot(results_summary, aes(x = Scenario, y = R)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = sprintf("%.2f", R)), vjust = -0.5) +
  ylim(0, max(results_summary$R, na.rm = TRUE) * 1.1) +
  ggtitle("Estimated R for Exported Infections by Testing Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot k values
p2 <- ggplot(results_summary, aes(x = Scenario, y = k)) +
  geom_bar(stat = "identity", fill = "darkred") +
  geom_text(aes(label = sprintf("%.2f", k)), vjust = -0.5) +
  ylim(0, max(results_summary$k, na.rm = TRUE) * 1.1) +
  ggtitle("Estimated k for Exported Infections by Testing Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add a plot to visualise the sensitivity function
t_values <- seq(0, 7, by = 0.1)
sensitivity_values <- sapply(t_values, sensitivity_function)

p3 <- ggplot(data.frame(t = t_values, sensitivity = sensitivity_values), aes(x = t, y = sensitivity)) +
  geom_line() +
  labs(title = "Time-varying Test Sensitivity", x = "Days since infection", y = "Test Sensitivity") +
  theme_minimal()

# Plot Relative Risk
p4 <- ggplot(results_summary, aes(x = Scenario, y = RelativeRisk)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f", RelativeRisk)), vjust = -0.5) +
  ylim(0, max(results_summary$RelativeRisk, na.rm = TRUE) * 1.1) +
  ggtitle("Relative Risk of Exported Infections by Testing Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")

# Arrange plots
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
