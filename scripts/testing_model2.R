
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

daily_flight_probability*100
#############################
# Setup Model Parameters
#############################

# Define transmission parameters
simulation_params <- expand_grid(
  R = c(2),
  k = c(1000)
) %>% 
  mutate(R_k_id = row_number())

# Generate initial branching process chains
initial_chains_cache <- "processing/initial_chains.qs"

if (file.exists(initial_chains_cache) && !overwrite_initial_chains) {
  message("Loading cached initial chains...")
  initial_chains <- qs::qread(initial_chains_cache)
} else {
  message("Generating new initial chains...")
  initial_chains <- generate_initial_chains(
    simulation_params, 
    n_chains = sim_params$n_chains, 
    stat_threshold = sim_params$stat_threshold,
    tmax = sim_params$tmax
  )
  
  # Save the initial chains
  qs::qsave(initial_chains, initial_chains_cache)
  message("Saved initial chains to cache.")
}

# Initialize flight chains
tic()
# Create processing directory if it doesn't exist
if (!dir.exists("processing")) {
  dir.create("processing")
}

# Process each R and k combination
all_flight_chains <- list()
n_flying_chains <- numeric(nrow(simulation_params))

# Define base flight parameters before initialization
base_flight_params <- list(flight_duration = 0.2)

for(param_idx in 1:nrow(simulation_params)) {
  # Define the cache file path for this parameter combination
  flight_chains_cache <- sprintf("processing/initial_flight_chains_R%.1f_k%.1f.qs", 
                               simulation_params$R[param_idx],
                               simulation_params$k[param_idx])
  
  # Load or generate flight chains for this parameter set
  if (file.exists(flight_chains_cache) && !overwrite_flight_chains) {
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
      tmax = sim_params$tmax,
      scenario = base_flight_params
    )
    
    # Save the flight chains
    qs::qsave(initial_flight_chains, flight_chains_cache)
    message("Saved flight chains to cache.")
  }
  
  # Store the flight chains
  all_flight_chains[[param_idx]] <- initial_flight_chains
  
  # Calculate number of unique chains that would fly for this parameter set
  n_flying_chains[param_idx] <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
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

# Create no testing baseline scenario
no_testing <- list(
  pre = 0,
  post = 0,
  pre_delay = NA,
  post_delay = NA,
  flight_duration = 0.2,
  quarantine_duration = 14,
  name = "A. No Testing",
  interventions_enacted = 0
)

# Create other scenarios
scenarios <- c(
  list(no_testing),  # Add baseline scenario first
  lapply(seq_along(intervention_days), function(i) {
    c(
      base_scenario,
      list(
        name = sprintf("%s. Pre-flight Day %d", LETTERS[i+1], intervention_days[i]),
        interventions_enacted = intervention_days[i]
      )
    )
  })
)

# Add debug messages for testing scenarios
message("\nTesting Scenarios:")
message("----------------")
for (scenario in scenarios) {
  message(sprintf("%s:", scenario$name))
  if (scenario$pre == 1) {
    message(sprintf("  - Pre-flight testing starts on day %d", scenario$interventions_enacted))
    message(sprintf("  - Test occurs %d day(s) before flight", scenario$pre_delay))
  }
  if (scenario$post == 1) {
    message(sprintf("  - Post-flight testing starts on day %d", scenario$interventions_enacted))
    message(sprintf("  - Test occurs %d day(s) after flight", scenario$post_delay))
  }
  if (scenario$pre == 0 && scenario$post == 0) {
    message("  - No testing implemented")
  }
  message("")
}

#############################
# Run Scenarios
#############################

tic()

# Create directory for scenario results if it doesn't exist
scenario_cache_dir <- "processing/scenario_results"
if (!dir.exists(scenario_cache_dir)) {
  dir.create(scenario_cache_dir, recursive = TRUE)
}

# Function to generate cache filename for a specific scenario and parameters
get_scenario_cache_file <- function(scenario_name, R, k) {
  clean_name <- gsub("[^[:alnum:]]", "_", scenario_name)
  file.path(scenario_cache_dir, 
            sprintf("scenario_%s_R%.1f_k%.1f.qs", clean_name, R, k))
}

# Initialize empty list for all results
all_results <- list()

for (scenario_idx in seq_along(scenarios)) {
  # Get current scenario
  current_scenario <- scenarios[[scenario_idx]]
  
  # Process each parameter combination
  for(param_idx in 1:nrow(simulation_params)) {
    # Generate cache filename for this scenario-parameter combination
    cache_file <- get_scenario_cache_file(
      current_scenario$name,
      simulation_params$R[param_idx],
      simulation_params$k[param_idx]
    )
    
    # Check if we should use cached results
    if (file.exists(cache_file) && !overwrite_scenarios) {
      message(sprintf("Loading cached results for scenario '%s' (R=%.1f, k=%.1f)...",
                     current_scenario$name,
                     simulation_params$R[param_idx],
                     simulation_params$k[param_idx]))
      param_results <- qs::qread(cache_file)
    } else {
      message(sprintf("Processing scenario '%s' (R=%.1f, k=%.1f)...",
                     current_scenario$name,
                     simulation_params$R[param_idx],
                     simulation_params$k[param_idx]))
      
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
          total_flying_chains = n_flying_chains[param_idx],
          scenario_name = current_scenario$name
        )
      
      # Cache the results for this scenario-parameter combination
      qs::qsave(param_results, cache_file)
    }
    
    # Append to results list
    all_results[[paste(current_scenario$name, param_idx)]] <- param_results
    
    # Clean up to free memory
    rm(param_results)
    gc()
  }
  
  # Print progress
  message(sprintf("Completed scenario %s - %s", 
                 LETTERS[scenario_idx], 
                 current_scenario$name))
}

# Combine all results at the end
message("Combining all results...")
all_results <- bind_rows(all_results)

# Save the combined results (optional, since we have individual caches)
final_results_cache <- "processing/all_scenario_results.qs"
message("Saving combined results to cache...")
qs::qsave(all_results, final_results_cache)

toc()
