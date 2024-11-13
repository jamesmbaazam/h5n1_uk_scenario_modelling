# Updated function to initialize flight variables and filter chains
initialize_flight_chains <- function(chain_data, daily_flight_probability, si_draws) {
  results <- chain_data %>%
    mutate(
      flight_time = time + sample(x = si_draws$y_rep, size = n(), replace = TRUE),
      will_fly = rbinom(n(), 1, daily_flight_probability) == 1,
      flight_time = if_else(will_fly, flight_time, Inf)
    ) %>% 
    arrange(chain, time)
  
  # Add ancestry column
  results <- add_ancestry_column(results)
  
  # Filter out non-flying individuals and their subsequent chains
  flying_chains <- results %>%
    group_by(chain) %>%
    filter(any(will_fly)) %>%  # Keep chains where at least one person flies
    ungroup()
  
  # Add indicator for infections in destination country
  flying_chains <- flying_chains %>%
    group_by(chain) %>%
    mutate(
      flown_ancestors = list(infectee[will_fly]),
      destination_infection = map_lgl(ancestry, ~any(.x %in% flown_ancestors[[1]]))
    ) %>%
    ungroup() %>%
    select(-flown_ancestors)  # Remove the temporary column
  
  return(flying_chains)
}

# Updated flight_test_fun to include flight_duration
flight_test_fun <- function(flying_chains, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay, flight_duration, quarantine_start_day) {
  
  # Add flight_end based on scenario-specific flight_duration
  flying_chains <- flying_chains %>%
    mutate(
      flight_end = flight_time + flight_duration,
      apply_quarantine = will_fly & (flight_time >= quarantine_start_day)
    )
  
  # Implement pre-flight testing (if applicable)
  if (test_before_flight == 1) {
    flying_chains <- flying_chains %>%
      mutate(
        pre_flight_test_time = if_else(apply_quarantine, flight_time - pre_flight_test_delay, Inf),
        time_since_infection_at_pre_test = pmax(0, pre_flight_test_time - time),
        pre_flight_test_sensitivity = map_dbl(time_since_infection_at_pre_test, sensitivity_function),
        pre_flight_test_result = rbinom(n(), 1, pre_flight_test_sensitivity),
        isolated_pre_flight = pre_flight_test_result == 1
      )
  } else {
    flying_chains <- flying_chains %>%
      mutate(
        pre_flight_test_time = Inf,
        pre_flight_test_result = 0,
        isolated_pre_flight = FALSE
      )
  }
  
  # Implement post-flight testing (if applicable)
  if (test_after_flight == 1) {
    flying_chains <- flying_chains %>%
      mutate(
        post_flight_test_time = if_else(apply_quarantine, flight_end + post_flight_test_delay, Inf),
        time_since_infection_at_post_test = pmax(0, post_flight_test_time - time),
        post_flight_test_sensitivity = map_dbl(time_since_infection_at_post_test, sensitivity_function),
        post_flight_test_result = rbinom(n(), 1, post_flight_test_sensitivity),
        isolated_post_flight = post_flight_test_result == 1
      )
  } else {
    flying_chains <- flying_chains %>%
      mutate(
        post_flight_test_time = Inf,
        post_flight_test_result = 0,
        isolated_post_flight = FALSE
      )
  }
  
  # Determine the earliest isolation time for each individual
  flying_chains <- flying_chains %>%
    mutate(
      isolation_time = pmin(
        if_else(isolated_pre_flight, pre_flight_test_time, Inf),
        if_else(isolated_post_flight, post_flight_test_time, Inf)
      )
    )
  
  # Identify isolated individuals
  isolated_individuals <- flying_chains %>%
    filter(isolation_time < Inf) %>%
    select(chain, infectee, isolation_time)
  
  # Prune chains based on isolated ancestors
  pruned_chains <- flying_chains %>%
    mutate(
      should_remove = map_lgl(seq_len(n()), function(i) {
        any(isolated_individuals$infectee %in% ancestry[[i]] &
              isolated_individuals$isolation_time < time[i] &
              isolated_individuals$chain == chain[i])
      })
    ) %>%
    filter(!should_remove) %>%
    select(-should_remove)
  
  # Apply the effects of isolation
  final_chains <- pruned_chains %>%
    mutate(
      prevented_flight = isolated_pre_flight & will_fly,
      will_fly = will_fly & !prevented_flight,
      flight_time = if_else(prevented_flight, Inf, flight_time),
      flight_end = if_else(prevented_flight, Inf, flight_end)
    )
  
  return(final_chains)
}

# Updated run_scenarios_on_chains function
run_scenarios_on_chains <- function(initial_chains_df, scenarios, daily_flight_probability, si_draws) {
  # Initialize flight chains once
  flying_chains <- initialize_flight_chains(initial_chains_df, daily_flight_probability, si_draws)
  
  # Run scenarios
  res <- map(scenarios, ~list(
    scenario = .$name,
    chains = flight_test_fun(flying_chains, .$pre, .$post, .$pre_delay, .$post_delay, .$flight_duration, .$quarantine_start_day)
  ))
  
  # Convert to data.frame
  map_df(res, ~ .x$chains %>% 
           mutate(scenario = .x$scenario) %>%
           select(scenario, everything()))
}



# Modified function to generate ancestry vector for an individual, considering chain
get_ancestry <- function(individual, chain, chain_data) {
  ancestry <- c()
  current <- individual
  current_chain <- chain
  while (TRUE) {
    ancestor_row <- chain_data[chain_data$infectee == current & chain_data$chain == current_chain, ]
    if (nrow(ancestor_row) == 0 || is.na(ancestor_row$infector)) break
    infector <- ancestor_row$infector
    ancestry <- c(ancestry, infector)
    current <- infector
  }
  return(rev(ancestry))  # Reverse to get oldest ancestor first
}

# Add ancestry column to the data frame
add_ancestry_column <- function(chain_data) {
  chain_data %>%
    rowwise() %>%
    mutate(ancestry = list(get_ancestry(infectee, chain, chain_data))) %>%
    ungroup()
}

# Function to prune the chain based on an isolated ancestor
prune_chain <- function(chain_data, isolated_ancestor) {
  chain_data %>%
    filter(!map_lgl(ancestry, ~isolated_ancestor %in% .))
}