library(tidyverse)
library(purrr)
library(gridExtra)
library(fitdistrplus)

# Set seed for reproducibility
set.seed(123)

#' Calculate proportion of runs that have controlled outbreak
#'
#' @author Joel Hellewell
#' @export
#' @inheritParams detect_extinct
extinct_prob <- function(outbreak_df_week = NULL, cap_cases  = NULL, week_range = 12:16) {
  
  n_sim <- max(outbreak_df_week$sim)
  
  extinct_runs <- detect_extinct(outbreak_df_week, cap_cases, week_range)
  out <-  sum(extinct_runs$extinct) / n_sim
  
  return(out)
}

#' Calculate whether outbreaks went extinct or not
#' @author Joel Hellewell
#' @param outbreak_df_week data.table  weekly cases produced by the outbreak model
#' @param cap_cases integer number of cumulative cases at which the branching process was terminated
#' @param week_range integer vector giving the (zero indexed) week range to test for whether an extinction occurred.
#' @importFrom data.table as.data.table fifelse
#'
#' @export
#'
detect_extinct <- function(outbreak_df_week  = NULL, cap_cases  = NULL, week_range = 12:16) {
  
  outbreak_df_week <- as.data.table(outbreak_df_week)
  outbreak_df_week <- outbreak_df_week[week %in% week_range]
  outbreak_df_week[, list(
    extinct = fifelse(all(weekly_cases == 0 & cumulative < cap_cases), 1, 0)
  ), by = sim]

}

# Define the flight_test_fun function
flight_test_fun <- function(chain_data, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay, flight_duration) {
  
  results <- chain_data %>%
    mutate(
      flight_time = sample(x = si_draws$y_rep, size = n(), replace = TRUE),
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
        pre_flight_test_sensitivity = map_dbl(pre_flight_test_t, sensitivity_function),
        pre_flight_test_res = rbinom(n(), 1, pre_flight_test_sensitivity),
        isolated_pre = if_else(pre_flight_test_res == 1, TRUE, isolated_pre)
      )
  }
  
  # Test after flight
  if (test_after_flight == 1) {
    results <- results %>%
      mutate(
        post_flight_test_t = flight_end + post_flight_test_delay,
        post_flight_test_sensitivity = map_dbl(post_flight_test_t, sensitivity_function),
        post_flight_test_res = rbinom(n(), 1, post_flight_test_sensitivity),
        isolated_post = if_else(post_flight_test_res == 1, TRUE, isolated_post)
      )
  }
  
  # Earliest isolation time
  results <- results %>%
    mutate(
      isolation_time = pmin(if_else(isolated_pre, pre_flight_test_t, Inf), if_else(isolated_post, post_flight_test_t, Inf))
    )
  
  # Indicate whether 2nd+ gen infections were averted if the infector was isolated
  # and if infection were to occur post-flight

  gen_2_averted <- results %>% 
    filter(generation  == 2) %>% 
    select(chain, infector, infectee, generation, flight_time, time) %>%
    right_join(results %>% filter(generation == 1) %>% 
                 select(chain, isolation_time), 
               by = "chain") %>%
    filter(time > isolation_time | time < flight_time)
  
  #if in the chain averted in gen 2, then remove those chains from gen 2 and above
  subsequent_chains <- results %>% 
    filter(generation >= 2) %>% 
    anti_join(gen_2_averted, by = "chain") %>% 
    select(chain, infector, infectee, generation, time, flight_time)
  
  # Combine 1st gen, 2nd gen, and gen 2+ chains
  filtered_chains <- results %>% 
    filter(generation == 1) %>% 
    select(chain, infector, infectee, generation, time, flight_time) %>% 
    bind_rows(subsequent_chains)
  
  return(filtered_chains)
}

# Function to calculate R and k using negative binomial fitting
calculate_r_and_k <- function(chain_data) {
  offspring_counts <- chain_data %>%
    mutate(across(c(chain, generation, infector), as.factor)) %>%
    group_by(chain, .drop = FALSE) %>%
    summarise(N = n(), .groups = "drop") %>%
    pull(N)
  
  if (all(offspring_counts == 0)) {
    warning("All offspring counts are zero. Unable to estimate R and k.")
    return(list(R = NA, k = NA))
  }
  
  tryCatch({
    nb_fit <- fitdistrplus::fitdist(offspring_counts, "nbinom")
    r <- nb_fit$estimate["mu"]
    k <- nb_fit$estimate["size"]
    return(list(R = r, k = k))
  }, error = function(e) {
    warning("Error in fitting negative binomial distribution: ", e$message)
    return(list(R = NA, k = NA))
  })
}

# Function to run a single scenario
run_scenario <- function(sim_chains_df, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay,  flight_duration, scenario_name) {
  filtered_results <- flight_test_fun(
    chain_data = sim_chains_df,
    test_before_flight = test_before_flight,
    test_after_flight = test_after_flight,
    pre_flight_test_delay = pre_flight_test_delay,
    post_flight_test_delay = post_flight_test_delay,
    flight_duration = flight_duration
  )
  
  return(filtered_results)
}

# Define scenarios
scenarios <- list(
  list(name = "No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_duration = 0.2),
  list(name = "Pre-flight (1 day before)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2),
  list(name = "Post-flight (1 day after)", pre = 0, post = 1, pre_delay = NA, post_delay = 1, flight_duration = 0.2),
  list(name = "Both (1 day before and after)", pre = 1, post = 1, pre_delay = 1, post_delay = 1, flight_duration = 0.2)
)

# Create a data frame with all combinations of R values and size parameters
simulation_params <- expand_grid(
  R = 1,
  size = 0.1
)

run_simulation <- function(R, size) {
  results <- map(1:100, function(i) {
    sim_chains <- simulate_chains(
      n_chains = 100,
      statistic = "length",
      offspring_dist = function(n, mu) rnbinom(n, mu = mu, size = size),
      stat_threshold = 15,
      generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE),
      mu = R
    )
    
    sim_chains_df <- as.data.frame(sim_chains)
    
    map(scenarios, ~list(
      scenario = .$name,
      chains = run_scenario(sim_chains_df, .$pre, .$post, .$pre_delay, .$post_delay, .$flight_duration, .$name)
    ))
  })
}

# Run simulations for all combinations of R and size
all_results <- simulation_params %>%
  mutate(results = map2(R, size, run_simulation)) %>%
  unnest(results)

# Function to process chains and calculate daily cases
process_chains <- function(chains, scenario) {
  chains %>%
    mutate(time = floor(time)) %>%
    filter(time > flight_time) %>%
    group_by(time) %>% 
    count() %>%
    mutate(Scenario = scenario)
}

# Process all chains
processed_chains <- all_results %>% 
  unnest(results) %>% 
  mutate(sim = row_number()) %>% 
  unnest_wider(results) %>%
  mutate(processed = map2(chains, scenario, process_chains)) %>%
  unnest(processed)

# Plot all simulations for each scenario
plot_all_simulations <- function(data) {
  ggplot(data, aes(x = time, y = n, group = interaction(Scenario, sim))) +
    geom_line(alpha = 0.1) +
    facet_wrap(~Scenario, scales = "free_y") +
    labs(x = "Time", y = "Daily Cases", title = "Daily Cases for 100 Simulations per Scenario") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Generate and display the plot
p <- plot_all_simulations(processed_chains)
print(p)
