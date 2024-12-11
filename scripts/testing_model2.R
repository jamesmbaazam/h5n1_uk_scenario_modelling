# Load required libraries
library(tidyverse)
library(purrr)
library(gridExtra)
library(fitdistrplus)
library(tictoc)
library(patchwork)

# Set seed for reproducibility
set.seed(123)

# Define simulation parameters
n_chains <- 10000 
tmax <- 1000
daily_flight_probability <- 0.025 # 

# Define all scenarios in one list
scenarios <- list(
  list(name = "A. No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_duration = 0.2, quarantine_start_day = Inf),
  list(name = "B. Pre-flight (1 day before)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 50),
  list(name = "C. Post-flight (1 day after)", pre = 0, post = 1, pre_delay = NA, post_delay = 1, flight_duration = 0.2, quarantine_start_day = 50),
  list(name = "D. Both (1 day before and after)", pre = 1, post = 1, pre_delay = 1, post_delay = 1, flight_duration = 0.2, quarantine_start_day = 50),
  list(name = "E. Pre-flight testing (Day 50)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 50),
  list(name = "F. Pre-flight testing (Day 100)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 100)
)

# Define all scenarios in one list
scenarios <- list(
  list(name = "A. No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_duration = 0.2, quarantine_start_day = Inf),
  list(name = "B. Pre-flight Day 0", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 0),
  list(name = "C. Pre-flight Day 25", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 25),
  list(name = "D. Pre-flight Day 50", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 50),
  list(name = "E. Pre-flight Day 75", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 75),
  list(name = "F. Pre-flight Day 100", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 100)
)

# Create a data frame with all combinations of R values and size parameters
simulation_params <- expand_grid(
  R = c(2),
  size = c(0.1)
) %>% 
  mutate(R_k_id = row_number())

# Generate initial chains
initial_chains <- generate_initial_chains(simulation_params, n_chains = n_chains, stat_threshold = 15)


# First, initialize flight chains once for all scenarios
tic()
initial_flight_chains <- initialize_flight_chains(initial_chains$initial_chain[[1]], daily_flight_probability, si_draws)
toc()

tic()
# Process scenarios one at a time and combine results
all_results <- list()

for (scenario_idx in seq_along(scenarios)) {
  current_scenario <- list(scenarios[[scenario_idx]])
  
  # Process single scenario for each chain
  scenario_results <- list()
  for (chain_idx in seq_len(nrow(initial_chains))) {
    chain_result <- run_scenarios_on_chains(
      initial_flight_chains,  # Use the pre-initialized flight chains
      current_scenario, 
      daily_flight_probability, 
      si_draws
    )
    
    # Add the R and size values from initial_chains
    chain_result <- chain_result %>%
      mutate(
        R = initial_chains$R[chain_idx],
        size = initial_chains$size[chain_idx]
      )
    
    scenario_results[[chain_idx]] <- chain_result
  }
  
  # Combine results for this scenario
  scenario_results <- bind_rows(scenario_results)
  
  # Store results
  all_results[[scenario_idx]] <- scenario_results
  
  # Print progress
  print(paste("Completed scenario", LETTERS[scenario_idx], "-", scenarios[[scenario_idx]]$name))
  
  # Clear intermediate objects
  rm(scenario_results, chain_result)
  gc()
}

# Combine all results
all_results <- bind_rows(all_results)
toc()

# Create plots using combined results
p_daily <- plot_daily_cases(all_results, "Set1")
p_cumulative <- plot_cumulative_cases(all_results, "Set1")
p_avg_cumulative <- plot_avg_cumulative_cases(all_results, "Set1")

# Save individual plots
ggsave("results/daily_cases.png", p_daily, width = 12, height = 8, bg = "white", dpi = 600)
ggsave("results/daily_cases.pdf", p_daily, width = 12, height = 8, bg = "white", dpi = 600)

ggsave("results/cumulative_cases.png", p_cumulative, width = 12, height = 8, bg = "white", dpi = 600)
ggsave("results/cumulative_cases.pdf", p_cumulative, width = 12, height = 8, bg = "white", dpi = 600)

ggsave("results/avg_cumulative_cases.png", p_avg_cumulative, width = 12, height = 8, bg = "white", dpi = 600)
ggsave("results/avg_cumulative_cases.pdf", p_avg_cumulative, width = 12, height = 8, bg = "white", dpi = 600)

# Calculate and plot additional metrics
time_to_100_cases <- calculate_time_to_100_cases(all_results)
outbreak_prob <- calculate_extinction_probability(all_results)

p_time_to_100 <- plot_time_to_100_cases(time_to_100_cases, "Set1")
p_outbreak_prob <- plot_outbreak_probability(outbreak_prob, "Set1")

# Save additional plots
ggsave("results/time_to_100_cases.png", p_time_to_100, width = 12, height = 8, bg = "white", dpi = 600)
ggsave("results/time_to_100_cases.pdf", p_time_to_100, width = 12, height = 8, bg = "white", dpi = 600)

ggsave("results/outbreak_probability.png", p_outbreak_prob, width = 12, height = 8, bg = "white", dpi = 600)
ggsave("results/outbreak_probability.pdf", p_outbreak_prob, width = 12, height = 8, bg = "white", dpi = 600)

print("Simulation and visualization complete. Check the results directory for output files.")

