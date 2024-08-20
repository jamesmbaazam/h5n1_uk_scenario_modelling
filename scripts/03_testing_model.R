# Set seed for reproducibility
set.seed(123)

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
  generation_time = function(n) rgamma(n, shape = si_shape, scale = si_scale),
  mu = 2,
  size = 1
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
