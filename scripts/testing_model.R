library(tidyverse)
library(purrr)
library(gridExtra)
library(fitdistrplus)

# Set seed for reproducibility
set.seed(123)

tmax = 1000

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
  list(name = "A. No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_duration = 0.2),
  list(name = "B. Pre-flight (1 day before)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2)
  #list(name = "C. Post-flight (1 day after)", pre = 0, post = 1, pre_delay = NA, post_delay = 1, flight_duration = 0.2),
  #list(name = "D. Both (1 day before and after)", pre = 1, post = 1, pre_delay = 1, post_delay = 1, flight_duration = 0.2)
)

# Create a data frame with all combinations of R values and size parameters
simulation_params <- expand_grid(
  R = c(1,2),
  size = c(0.1,1)
) %>% 
  mutate(R_k_id = row_number())

run_simulation <- function(R, size) {
    sim_chains <- simulate_chains(
      n_chains = 100,
      statistic = "length",
      offspring_dist = function(n, mu) rnbinom(n, mu = R, size = size),
      stat_threshold = 20,
      generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE)
    )
    
    
    sim_chains_df <- as.data.frame(sim_chains)
    
    res <- map(scenarios, ~list(
      scenario = .$name,
      chains = run_scenario(sim_chains_df, .$pre, .$post, .$pre_delay, .$post_delay, .$flight_duration, .$name)
    ))
    
    # Convert to data.frame
    result <- res %>%
      map_df(~ .x$chains %>% 
               mutate(scenario = .x$scenario) %>%
               select(scenario, everything()))
}

# Run simulations for all combinations of R and size
tic()
all_results <- simulation_params %>%
  mutate(results = map2(R, size, run_simulation))  %>%
  unnest(results) 
toc()


# Function to aggregate chains and calculate daily cases
aggregate_chains <- function(data) {
browser()
  data %>%
    mutate(time = floor(time)) %>%
    filter(time > flight_time) %>%
    group_by(time, .drop = FALSE) %>% 
    count() %>% 
    ungroup() %>% 
    complete(time = 1:tmax, fill = list(n = 0)) %>% 
    mutate(Scenario = scenario) 
}

# Plot aggregated daily cases
all_results %>% 
  filter(time > flight_time) %>%
  mutate(day = floor(time)) %>%
  group_by(scenario, R, size, rep, day) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = day, y = n, color = scenario)) +
  geom_line(alpha = 0.2) +
  facet_grid(R~size) +
  labs(x = "Day", y = "Daily Cases", title = "Daily Cases by Scenario and R_k_id") 

#Plot cumulative cases
all_results %>% 
  filter(time > flight_time) %>%
  mutate(day = floor(time)) %>%
  group_by(scenario, R, size, rep, chain, day) %>%
  summarise(n = n()) %>%
  arrange(day) %>%
  ggplot(aes(x = day, y = cumsum(n), group = interaction(rep,chain))) +
  geom_line(alpha = 0.2) +
  facet_grid(R~size) +
  labs(x = "Day", y = "Cumulative Cases", title = "Cumulative Cases by Scenario and R_k_id")

#Plot individual chains
all_results %>% 
  filter(time > flight_time) %>%
  mutate(day = floor(time)) %>%
  group_by(scenario, R, size, rep, day, chain) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = day, y = n, group = interaction(rep,chain))) +
  geom_line(alpha = 0.2) +
  facet_grid(R~size) +
  labs(x = "Day", y = "Daily Cases", title = "Daily Cases by Scenario and R_k_id")

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
plot_all_simulations(processed_chains)

# Function to calculate the proportion of extinct simulations per day
calculate_extinct_proportion <- function(data) {
  
  data %>%
    group_by(Scenario, time) %>%
    summarise(
      total_sims = 100,
      extinct_sims = sum(n == 0),
      proportion_extinct = extinct_sims / total_sims,
      .groups = "drop"
    )
}

# Calculate the proportion of extinct simulations
extinct_proportions <- processed_chains %>%
  calculate_extinct_proportion()

# Plot the proportion of extinct simulations
plot_extinct_proportions <- function(data) {
  ggplot(data, aes(x = time, y = proportion_extinct, color = Scenario)) +
    geom_line() +
    facet_wrap(~Scenario) +
    labs(x = "Time", y = "Proportion of Extinct Simulations", 
         title = "Proportion of Extinct Simulations Over Time") +
    theme_minimal() +
    theme(legend.position = "none") +
    ylim(0, 1)
}

# Generate and display the extinction plot
p_extinct <- plot_extinct_proportions(extinct_proportions)
print(p_extinct)


