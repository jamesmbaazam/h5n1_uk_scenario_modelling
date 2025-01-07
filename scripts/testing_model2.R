
#############################
# Initialize Simulation
#############################

# Set random seed for reproducibility
set.seed(123)

# Define core simulation parameters
sim_params <- list(
  n_chains = 10000,
  tmax = 1000,
  stat_threshold = 12
)

# Define constant sensitivity function
sensitivity_function <- function(x) {
  return(0.4)
}

#############################
# Calculate Flight Probability
#############################

# Source data from official statistics
flight_stats <- list(
  # IPS ONS survey: https://www.visitbritain.org/research-insights/inbound-visits-and-spend-annual-uk
  us_uk_visitors_23 = 5122000,
  
  # US Census Bureau: https://www.census.gov/library/stories/2022/12/happy-new-year-2023.html
  us_population = 334233854
)

# Calculate daily probability of US->UK flight
daily_flight_probability <- (flight_stats$us_uk_visitors_23 / flight_stats$us_population) / 365

#############################
# Setup Model Parameters
#############################

# Define transmission parameters
simulation_params <- expand_grid(
  R = c(1.5,2),
  k = c(0.1,1000)
) %>% 
  mutate(R_k_id = row_number())

# Generate initial branching process chains
initial_chains <- generate_initial_chains(
  simulation_params, 
  n_chains = sim_params$n_chains, 
  stat_threshold = sim_params$stat_threshold,
  tmax = sim_params$tmax
)

# Initialize flight chains
tic()
# Create processing directory if it doesn't exist
if (!dir.exists("processing")) {
  dir.create("processing")
}

# Set overwrite flag for cache control
overwrite <- TRUE  # Change to TRUE to force regeneration of flight chains

# Process each R and k combination
all_flight_chains <- list()
n_flying_chains <- numeric(nrow(simulation_params))

for(param_idx in 1:nrow(simulation_params)) {
  # Define the cache file path for this parameter combination
  flight_chains_cache <- sprintf("processing/initial_flight_chains_R%.1f_k%.1f.qs", 
                               simulation_params$R[param_idx],
                               simulation_params$k[param_idx])
  
  # Load or generate flight chains for this parameter set
  if (file.exists(flight_chains_cache) && !overwrite) {
    message(sprintf("Loading cached flight chains for R=%.1f, k=%.1f...",
                   simulation_params$R[param_idx],
                   simulation_params$k[param_idx]))
    initial_flight_chains <- qs::qread(flight_chains_cache)
  } else {
    message(sprintf("Generating new flight chains for R=%.1f, k=%.1f...",
                   simulation_params$R[param_idx],
                   simulation_params$k[param_idx]))
    initial_flight_chains <- initialize_flight_chains(
      initial_chains$initial_chain[[param_idx]], 
      daily_flight_probability, 
      si_draws,
      tmax = sim_params$tmax
    )
    
    # Save the flight chains
    qs::qsave(initial_flight_chains, flight_chains_cache)
    message("Saved flight chains to cache.")
  }
  
  # Store the flight chains
  all_flight_chains[[param_idx]] <- initial_flight_chains
  
  # Calculate number of unique chains that would fly for this parameter set
  n_flying_chains[param_idx] <- initial_flight_chains %>%
    filter(will_fly == TRUE) %>%
    pull(chain) %>%
    n_distinct()
  
  message(sprintf("Found %d unique flying chains for R=%.1f, k=%.1f", 
                 n_flying_chains[param_idx],
                 simulation_params$R[param_idx],
                 simulation_params$k[param_idx]))
}

toc()

#############################
# Define Intervention Scenarios
#############################

# Common parameters across scenarios
base_scenario <- list(
  pre = 1,
  post = 0,
  pre_delay = 1,
  post_delay = NA,
  flight_duration = 0.2,
  quarantine_duration = 14
)

# Define scenarios with different intervention start times
intervention_days <- c(0, 25, 50, 75, 100)
scenarios <- lapply(seq_along(intervention_days), function(i) {
  c(
    base_scenario,
    list(
      name = sprintf("%s. Pre-flight Day %d", LETTERS[i+1], intervention_days[i]),
      interventions_enacted = intervention_days[i]
    )
  )
})

#############################
# Run Scenarios
#############################

# Process each scenario
tic()
all_results <- list()

for (scenario_idx in seq_along(scenarios)) {
  # Get current scenario
  current_scenario <- scenarios[[scenario_idx]]
  
  scenario_results <- list()
  
  # Process each parameter combination
  for(param_idx in 1:nrow(simulation_params)) {
    # Process scenario for this parameter set
    param_results <- process_scenario(
      all_flight_chains[[param_idx]],
      current_scenario,
      quarantine_duration = current_scenario$quarantine_duration,
      tmax = sim_params$tmax
    ) %>%
      # Add transmission parameters and flying chains count
      mutate(
        R = simulation_params$R[param_idx],
        k = simulation_params$k[param_idx],
        total_flying_chains = n_flying_chains[param_idx]
      )
    
    scenario_results[[param_idx]] <- param_results
  }
  
  # Combine results for all parameter sets
  all_results[[scenario_idx]] <- bind_rows(scenario_results)
  
  # Print progress
  print(sprintf("Completed scenario %s - %s", 
                LETTERS[scenario_idx+1], 
                current_scenario$name))
}

# Combine results from all scenarios
all_results <- bind_rows(all_results)
toc()
