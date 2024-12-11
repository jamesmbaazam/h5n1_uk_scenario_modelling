# Modified initialize_flight_chains function
initialize_flight_chains <- function(chain_data, daily_flight_probability, si_draws) {
  
  # Get unique chains
  unique_chains <- unique(chain_data$chain)
  n_unique_chains <- length(unique_chains)
  
  # First, randomly assign flight status to all chains
  chain_flight_status <- tibble(
    chain = unique_chains,
    will_fly_chain = rbinom(n_unique_chains, 1, daily_flight_probability) == 1 
  )
  
  # Filter chains first based on flight status
  filtered_chains <- chain_data %>%
    left_join(chain_flight_status, by = "chain") %>%
    filter(will_fly_chain) %>%
    select(-will_fly_chain)
  
  # Then process only the filtered chains
  results <- filtered_chains %>%
    mutate(
      flight_time = time + sample(x = si_draws$y_rep, size = n(), replace = TRUE),
      will_fly = TRUE,  # All remaining chains will fly
      flight_time = flight_time  # No need for if_else anymore
    ) %>% 
    arrange(chain, time)
  
  # Add ancestry column
  results <- add_ancestry_column(results)
  
  # Add indicator for infections in destination country
  results <- results %>%
    group_by(chain) %>%
    mutate(
      flown_ancestors = list(infectee),  # All infectees in remaining chains will fly
      destination_infection = map_lgl(ancestry, ~TRUE)  # All infections in remaining chains are destination infections
    ) %>%
    ungroup() %>%
    select(-flown_ancestors)
  
  return(results)
}

# Modified generate_initial_chains function
generate_initial_chains <- function(simulation_params, n_chains, stat_threshold) {
  
  simulation_params %>%
    mutate(
      initial_chain = map(R_k_id, ~ simulate_chains(
        n_chains = n_chains,
        statistic = "length",
        offspring_dist = function(n, mu) rnbinom(n, mu = 1.5, size = 0.1),
        stat_threshold = stat_threshold,
        generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE),
        tf = tmax
      ))
    )
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
run_scenarios_on_chains <- function(flying_chains, scenarios, daily_flight_probability, si_draws) {
  # Run scenarios directly on the pre-initialized flight chains
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

##### Plotting functions #####
# Function to create daily cases plot
plot_daily_cases <- function(data, color_palette) {
  data %>% 
    filter(destination_infection) %>% 
    mutate(day = floor(time)) %>%
    group_by(scenario, R, size, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    ggplot(aes(x = day, y = n, color = scenario, group = interaction(chain, scenario))) +
    geom_line(alpha = 0.2) +
    scale_colour_brewer(palette = color_palette) +
    facet_grid(R ~ size) +
    labs(x = "Day", y = "Daily Cases", title = "Daily Cases by Scenario") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to create cumulative cases plot
plot_cumulative_cases <- function(data, color_palette) {
  data %>% 
    mutate(day = floor(time)) %>%
    group_by(scenario, R, size, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(day = 0:1000, scenario, R, size, chain, fill = list(n = 0)) %>%
    group_by(scenario, R, size, chain) %>%
    mutate(cumsum_n = cumsum(n)) %>% 
    ggplot(aes(x = day, y = cumsum_n, group = interaction(chain, scenario), colour = scenario)) +
    geom_line(alpha = 0.2) +
    scale_colour_brewer(palette = color_palette) +
    facet_grid(R ~ size) +
    labs(x = "Day", y = "Cumulative Cases", title = "Cumulative Cases by Scenario") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to create average cumulative cases plot
plot_avg_cumulative_cases <- function(data, color_palette) {
  avg_data <- data %>%
    mutate(day = floor(time)) %>%
    group_by(scenario, R, size, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(day = 0:1000, scenario, R, size, chain, fill = list(n = 0)) %>%
    group_by(scenario, R, size, chain) %>%
    mutate(cumsum_n = cumsum(n)) %>%
    group_by(scenario, R, size, day) %>%
    summarise(avg_cumsum = mean(cumsum_n),
              low_ci = quantile(cumsum_n, 0.025),
              high_ci = quantile(cumsum_n, 0.975),
              .groups = "drop")
  
  ggplot(avg_data, aes(x = day, y = avg_cumsum, color = scenario)) +
    geom_line() +
    #geom_ribbon(aes(ymin = low_ci, ymax = high_ci, fill = scenario), alpha = 0.2, colour = NA) +
    scale_colour_brewer(palette = color_palette) +
    scale_fill_brewer(palette = color_palette) +
    facet_grid(R ~ size) +
    labs(x = "Day", y = "Average Cumulative Cases", title = "Average Cumulative Cases by Scenario") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to calculate time to reach 100 cases with uncertainty
calculate_time_to_100_cases <- function(data) {
  # Function to calculate day of 100 cases for a chain
  calc_day_100 <- function(chain_data) {
    result <- chain_data %>%
      filter(destination_infection) %>%
      mutate(day = floor(time)) %>%
      group_by(day) %>%
      summarise(daily_cases = n(), .groups = "drop") %>%
      arrange(day) %>%
      mutate(cumulative_cases = cumsum(daily_cases)) %>%
      filter(cumulative_cases >= 100) %>%
      slice_min(day, n = 1) %>%
      pull(day)
    
    # Return NA if no day reaches 100 cases
    if (length(result) == 0) return(NA)
    return(result)
  }
  
  # Calculate days to 100 cases for each chain
  results <- data %>%
    group_by(scenario, R, size, chain) %>%
    nest() %>%
    mutate(day_100 = map_dbl(data, ~calc_day_100(.))) %>%
    group_by(scenario, R, size) %>%
    summarise(
      actual_results = median(day_100, na.rm = TRUE),
      lower_ci = quantile(day_100, probs = 0.025, na.rm = TRUE),
      upper_ci = quantile(day_100, probs = 0.975, na.rm = TRUE),
      n_chains_reached_100 = sum(!is.na(day_100)),
      total_chains = n(),
      .groups = "drop"
    )
  
  # Calculate actual cumulative cases on the median day
  results <- results %>%
    rowwise() %>%
    mutate(
      cumulative_cases = data %>%
        filter(destination_infection) %>%
        mutate(day = floor(time)) %>%
        filter(day <= actual_results) %>%
        nrow()
    )
  
  return(results)
}

# Update the plotting function to show proportion of chains reaching 100 cases
plot_time_to_100_cases <- function(data, color_palette) {
  ggplot(data, aes(x = scenario, y = actual_results, colour = scenario)) +
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), size = 1) +
    geom_text(aes(label = sprintf("Day %d\n(%d cases)\n[CI: %d-%d]\n%d/%d chains", 
                                 round(actual_results), 
                                 cumulative_cases,
                                 round(lower_ci),
                                 round(upper_ci),
                                 n_chains_reached_100,
                                 total_chains)), 
              vjust = -0.5, size = 3) +
    scale_colour_brewer(palette = color_palette) +
    facet_grid(R ~ size) +
    labs(x = "Scenario", 
         y = "Time (days)", 
         title = "Time to Reach 100 Total Cumulative Cases",
         subtitle = "Points show median day, error bars show 95% quantiles across chains") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
}

# Function to calculate extinction probability
calculate_extinction_probability <- function(data) {
  
  data %>%
    mutate(day  = floor(time)) %>%
    group_by(scenario, R, size, chain, day) %>%
    summarise(n = n(),
              extinct = n == 0,
              .groups = "drop") %>%
    group_by(scenario, R, size, day) %>%
    summarise(
      extinction_prob = mean(extinct),
      lower_ci = binom.test(sum(extinct), n())$conf.int[1],
      upper_ci = binom.test(sum(extinct), n())$conf.int[2],
      .groups = "drop"
    ) %>%
    mutate(
      outbreak_prob = 1 - extinction_prob,
      outbreak_lower_ci = 1 - upper_ci,
      outbreak_upper_ci = 1 - lower_ci
    )
}

# Function to plot outbreak probability
plot_outbreak_probability <- function(data, color_palette) {
  ggplot(data, aes(x = day, y = outbreak_prob, colour = scenario)) +
    geom_step()+
    geom_ribbon(aes(ymin = outbreak_lower_ci, ymax = outbreak_upper_ci, fill = scenario), alpha = 0.2, colour = NA) +
    scale_fill_brewer(palette = color_palette) +
    facet_grid(R ~ size+scenario) +
    labs(x = "Scenario", y = "Probability", title = "Probability of Outbreak (1 - Extinction Probability)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
}
