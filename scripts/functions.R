# Modified generate_initial_chains function
generate_initial_chains <- function(simulation_params, n_chains, stat_threshold, tmax) {
  # Convert to data.table and ensure copy
  dt <- data.table::as.data.table(simulation_params)
  data.table::setDT(dt)  # Ensure it's a data.table
  
  # Create list to store chains
  chain_list <- vector("list", nrow(dt))
  
  # Generate chains for each parameter combination
  for(i in seq_len(nrow(dt))) {
    current_R <- dt$R[i]
    current_k <- dt$k[i]
    
    chain_list[[i]] <- simulate_chains(
      n_chains = n_chains,
      statistic = "length",
      offspring_dist = function(n, mu) {
        rnbinom(n, mu = current_R, size = current_k)
      },
      stat_threshold = stat_threshold,
      generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE),
      tf = tmax
    )
  }
  
  # Add chains to data.table
  dt[, initial_chain := chain_list]
  
  # Convert back to data.frame if needed
  return(as.data.frame(dt))
}

initialize_flight_chains <- function(chain_data, daily_flight_probability, si_draws, tmax) {
  message("Initializing flight chains...")
  dt <- data.table::as.data.table(chain_data)
  n_cases <- nrow(dt)
  message(sprintf("Processing %d cases", n_cases))
  
  # Initial flight assignments
  dt[, will_fly := rbinom(n_cases, 1, daily_flight_probability) == 1]
  n_flyers <- sum(dt$will_fly)
  message(sprintf("Identified %d potential flyers (%.1f%%)", n_flyers, 100*n_flyers/n_cases))
  
  # Exit if no flights
  if(!any(dt$will_fly)) {
    message("No flights found, returning empty dataset")
    return(list(
      data = data.table(
        chain = integer(0),
        infectee = integer(0),
        will_fly = logical(0),
        time = numeric(0),
        flight_time = numeric(0),
        flight_end = numeric(0)
      ),
      n_flying_chains = 0
    ))
  }
  
  # Extract potential flyers
  flyers <- dt[will_fly == TRUE]
  
  # For each flyer, trace back to find if they have flying ancestors
  setkey(dt, chain, infectee)
  flyers[, is_first_flight := TRUE]
  
  for(i in 1:nrow(flyers)) {
    if(i %% 10 == 0) message(sprintf("Processing flyer %d of %d", i, nrow(flyers)))
    
    current_chain <- flyers[i, chain]
    current_case <- flyers[i, infectee]
    
    # Trace back until we find a flying ancestor or hit root
    while(!is.na(dt[.(current_chain, current_case), infector])) {
      current_case <- dt[.(current_chain, current_case), infector]
      if(dt[.(current_chain, current_case), will_fly]) {
        flyers[i, is_first_flight := FALSE]
        break
      }
    }
  }
  
  # Keep only first flyers and set their flight times
  first_flyers <- flyers[is_first_flight == TRUE]
  first_flyers[, ':='(
    flight_time = time + sample(si_draws$y_rep, .N, replace = TRUE),
    flight_end = time + sample(si_draws$y_rep, .N, replace = TRUE) + 0.2
  )]
  
  message(sprintf("Found %d first flyers", nrow(first_flyers)))
  
  # Create index of cases to keep
  keep_idx <- vector("logical", nrow(dt))
  
  # Process each first flyer separately
  for(i in 1:nrow(first_flyers)) {
    current_chain <- first_flyers[i, chain]
    flyer_id <- first_flyers[i, infectee]
    
    # Mark flyer
    keep_idx[dt[.(current_chain, flyer_id), which = TRUE]] <- TRUE
    
    # Mark all descendants
    current_ids <- flyer_id
    while(length(current_ids) > 0) {
      # Find next generation
      descendant_rows <- dt[chain == current_chain & infector %in% current_ids, which = TRUE]
      if(length(descendant_rows) == 0) break
      
      # Mark these descendants
      keep_idx[descendant_rows] <- TRUE
      
      # Move to next generation
      current_ids <- dt[descendant_rows, infectee]
    }
  }
  
  # Keep marked cases
  result <- dt[keep_idx]
  
  # Add flight times from first flyers
  result[, flight_time := NA_real_]
  result[, flight_end := NA_real_]
  
  # Update flight times for each chain
  for(i in 1:nrow(first_flyers)) {
    current_chain <- first_flyers[i, chain]
    result[chain == current_chain, ':='(
      flight_time = first_flyers[i, flight_time],
      flight_end = first_flyers[i, flight_end]
    )]
  }
  
  # Calculate destination infections
  result[, destination_infection := time > flight_end]
  
  # Sort
  setorder(result, chain, generation, time)
  
  return(result)
}

process_scenario <- function(flying_chains, scenario, quarantine_duration = 14, tmax) {
  message("Processing scenario: ", scenario$name)
  
  # Convert to data.table if not already
  dt <- if(data.table::is.data.table(flying_chains)) {
    data.table::copy(flying_chains)
  } else {
    data.table::as.data.table(flying_chains)
  }
  
  # Print column names for debugging
  message("Available columns: ", paste(names(dt), collapse = ", "))
  
  # Check if required columns exist
  required_cols <- c("chain", "infectee", "will_fly", "time", "flight_time", "flight_end")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Set key after ensuring columns exist
  data.table::setkey(dt, chain, infectee)
  
  # Extract flyers only
  flyers <- dt[will_fly == TRUE]
  message(sprintf("Processing %d flyers", nrow(flyers)))
  
  # Process flyers first
  flyers[, ':='(
    original_flight_time = flight_time,
    original_flight_end = flight_end,
    interventions_active = time >= scenario$interventions_enacted,
    pre_flight_test_time = Inf,
    post_flight_test_time = Inf,
    isolation_time = Inf,
    quarantine_end = Inf,
    prevented_flight = FALSE
  )]
  
  # Pre-flight testing (only for flyers)
  if(scenario$pre == 1) {
    message("Applying pre-flight testing")
    flyers[interventions_active == TRUE, ':='(
      pre_flight_test_time = pmax(0, original_flight_time - scenario$pre_delay),
      time_since_infection_pre = pmax(0, pre_flight_test_time - time)
    )]
    
    # Calculate test sensitivity
    flyers[interventions_active == TRUE, 
           pre_flight_test_sensitivity := sensitivity_function(time_since_infection_pre)]
    
    # Generate test results
    flyers[interventions_active == TRUE, ':='(
      pre_flight_test_result = rbinom(.N, 1, 
                                    pmin(1, pmax(0, pre_flight_test_sensitivity)))
    )]
    
    # Mark prevented flights
    flyers[interventions_active == TRUE & pre_flight_test_result == 1, 
           prevented_flight := TRUE]
  }
  
  # Post-flight testing (only for non-prevented flights)
  if(scenario$post == 1) {
    message("Applying post-flight testing")
    flyers[interventions_active == TRUE & prevented_flight == FALSE, ':='(
      post_flight_test_time = pmax(0, original_flight_end + scenario$post_delay),
      time_since_infection_post = pmax(0, post_flight_test_time - time)
    )]
    
    # Calculate test sensitivity
    flyers[interventions_active == TRUE & prevented_flight == FALSE,
           post_flight_test_sensitivity := sensitivity_function(time_since_infection_post)]
    
    # Generate test results
    flyers[interventions_active == TRUE & prevented_flight == FALSE, ':='(
      post_flight_test_result = rbinom(.N, 1, 
                                     pmin(1, pmax(0, post_flight_test_sensitivity)))
    )]
    
    # Set quarantine periods for positive tests
    flyers[interventions_active == TRUE & 
           prevented_flight == FALSE & 
           post_flight_test_result == 1,
           quarantine_end := post_flight_test_time + quarantine_duration]
  }
  
  # Update flight status for flyers
  flyers[prevented_flight == TRUE, ':='(
    will_fly = FALSE,
    flight_time = Inf,
    flight_end = Inf
  )]
  
  # Create index of cases to keep
  keep_idx <- vector("logical", nrow(dt))
  
  # Process each non-prevented flyer
  message("Processing descendants")
  active_flyers <- flyers[prevented_flight == FALSE]
  
  # Process descendants
  for(i in 1:nrow(active_flyers)) {
    if(i %% 100 == 0) message(sprintf("Processing flyer %d of %d", i, nrow(active_flyers)))
    
    current_chain <- active_flyers[i, chain]
    flyer_id <- active_flyers[i, infectee]
    
    # Mark flyer using direct index
    matching_rows <- dt[.(current_chain, flyer_id), which = TRUE]
    keep_idx[matching_rows] <- TRUE
    
    # Mark all descendants
    current_ids <- flyer_id
    while(length(current_ids) > 0) {
      # Find next generation
      descendant_idx <- dt[chain == current_chain & infector %in% current_ids, which = TRUE]
      if(length(descendant_idx) == 0) break
      
      # Mark these descendants
      keep_idx[descendant_idx] <- TRUE
      
      # Move to next generation
      current_ids <- dt[descendant_idx, infectee]
    }
  }
  
  # Keep marked cases
  result <- dt[keep_idx]
  
  # Join flight and quarantine information from flyers
  result[, ':='(
    flight_time = NA_real_,
    flight_end = NA_real_,
    quarantine_end = Inf
  )]
  
  # Update information for each chain
  setkey(active_flyers, chain)  # Set key for active_flyers
  for(i in 1:nrow(active_flyers)) {
    current_chain <- active_flyers[i, chain]
    result[chain == current_chain, ':='(
      flight_time = active_flyers[i, flight_time],
      flight_end = active_flyers[i, flight_end],
      quarantine_end = active_flyers[i, quarantine_end]
    )]
  }
  
 # Calculate destination infections
  result[, ':='(
    destination_infection = 
      !is.na(flight_end) & 
      time > flight_end & 
      (time < quarantine_end | quarantine_end == Inf),
    scenario = scenario$name  # Add scenario name here
  )]
  
   # Sort
  setorder(result, chain, generation, time)
  
  message(sprintf("Final dataset has %d cases", nrow(result)))
  return(as.data.frame(result))
}

##### Plotting functions #####
# Function to create daily cases plot
plot_daily_cases <- function(data, initial_flight_chains, color_palette) {
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(will_fly == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  data %>% 
    filter(destination_infection) %>% 
    mutate(day = floor(time)) %>%
    group_by(scenario, R, k, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario, R, k),
      chain = flying_chains,
      fill = list(n = 0)
    ) %>%
    ggplot(aes(x = day, y = n, color = scenario, group = interaction(chain, scenario))) +
    geom_line(alpha = 0.5) +
    scale_colour_brewer(palette = color_palette) +
    facet_grid(R ~ k) +
    labs(x = "Day", 
         y = "Daily Cases", 
         title = "Daily Cases by Scenario",
         subtitle = sprintf("Based on %d possible flying chains", total_flying_chains)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to create cumulative cases plot
plot_cumulative_cases <- function(data, initial_flight_chains, color_palette) {
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(will_fly == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  data %>% 
    filter(destination_infection) %>% 
    mutate(day = floor(time)) %>%
    group_by(scenario, R, k, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario, R, k),
      chain = flying_chains,
      fill = list(n = 0)
    ) %>%
    group_by(scenario, R, k, chain) %>%
    mutate(cumsum_n = cumsum(n)) %>%
    ggplot(aes(x = day, y = cumsum_n, group = interaction(chain, scenario), colour = scenario)) +
    geom_line(alpha = 0.2) +
    scale_colour_brewer(palette = color_palette) +
    facet_grid(R ~ k) +
    labs(x = "Day", 
         y = "Cumulative Cases", 
         title = "Cumulative Cases by Scenario",
         subtitle = sprintf("Based on %d possible flying chains", total_flying_chains)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to create average cumulative cases plot
plot_avg_cumulative_cases <- function(data, initial_flight_chains, color_palette) {
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(will_fly == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  avg_data <- data %>%
    filter(destination_infection) %>% 
    mutate(day = floor(time)) %>%
    group_by(scenario, R, k, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario, R, k),
      chain = flying_chains,
      fill = list(n = 0)  # Chains prevented from flying will have zero cases
    ) %>%
    group_by(scenario, R, k, chain) %>%
    mutate(cumsum_n = cumsum(n)) %>%
    group_by(scenario, R, k, day) %>%
    summarise(
      avg_cumsum = mean(cumsum_n),
      low_ci = quantile(cumsum_n, 0.025),
      high_ci = quantile(cumsum_n, 0.975),
      .groups = "drop"
    )
  
  ggplot(avg_data, aes(x = day, y = avg_cumsum, color = scenario)) +
    geom_line() +
    #geom_ribbon(aes(ymin = low_ci, ymax = high_ci, fill = scenario), alpha = 0.2, colour = NA) +
    scale_colour_brewer(palette = color_palette) +
    scale_fill_brewer(palette = color_palette) +
    facet_grid(R ~ k) +
    labs(x = "Day", 
         y = "Average Cumulative Cases", 
         title = "Average Cumulative Cases by Scenario",
         subtitle = sprintf("Based on %d possible flying chains", total_flying_chains)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to calculate time to reach N cumulative cases using survival analysis
calculate_time_to_n_cases <- function(data, initial_flight_chains, n_cases = 10) {
  require(survival)
  
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(will_fly == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  # Calculate time to threshold for each chain
  time_to_threshold <- data %>%
    filter(destination_infection) %>%
    mutate(day = floor(time)) %>%
    group_by(scenario, R, k, chain, day) %>%
    summarise(daily_cases = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    tidyr::complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      tidyr::nesting(scenario, R, k),
      chain = flying_chains,
      fill = list(daily_cases = 0)
    ) %>%
    group_by(scenario, R, k, chain) %>%
    arrange(day) %>%
    mutate(cum_cases = cumsum(daily_cases)) %>%
    summarise(
      time = if(max(cum_cases) >= n_cases) {
        as.numeric(min(day[cum_cases >= n_cases]))
      } else {
        as.numeric(max(day))
      },
      status = max(cum_cases) >= n_cases,
      .groups = "drop"
    )
  
  # Calculate survival statistics for each scenario/R/k combination
  results <- time_to_threshold %>%
    group_by(scenario, R, k) %>%
    summarise(
      n_chains_reached_threshold = sum(status),
      total_flying_chains = n(),
      proportion_reached = mean(status),
      cumulative_cases_threshold = n_cases,
      mean_time = mean(time[status], na.rm = TRUE),
      median_time = median(time[status], na.rm = TRUE),
      sd_time = sd(time[status], na.rm = TRUE),
      surv_fit = list(survfit(Surv(time, status) ~ 1)),
      .groups = "drop"
    ) %>%
    mutate(
      surv_summary = map(surv_fit, summary),
      actual_results = map_dbl(surv_summary, ~.x$table["median"]),
      lower_ci = map_dbl(surv_summary, ~.x$table["0.95LCL"]),
      upper_ci = map_dbl(surv_summary, ~.x$table["0.95UCL"]),
      rmst = map_dbl(surv_summary, ~.x$table["rmean"]),
      rmst_se = map_dbl(surv_summary, ~.x$table["se(rmean)"])
    ) %>%
    mutate(
      # Add time percentiles
      time_to_25pct = map_dbl(surv_fit, ~quantile(.x, probs = 0.25)$quantile),
      time_to_50pct = map_dbl(surv_fit, ~quantile(.x, probs = 0.50)$quantile),
      time_to_75pct = map_dbl(surv_fit, ~quantile(.x, probs = 0.75)$quantile),
      
      # Add comparative statistics
      delay_vs_baseline = if_else(scenario == "B. Pre-flight Day 0", 
                                  0, 
                                  actual_results - first(actual_results)),
      
      risk_reduction = if_else(scenario == "B. Pre-flight Day 0",
                               0,
                               (first(proportion_reached) - proportion_reached) / 
                                 first(proportion_reached))
    )
  
  return(list(
    results = results,
    time_to_threshold = time_to_threshold
  ))
}

# Update the plotting function to show survival curves - remove debug prints
plot_time_to_n_cases <- function(data, color_palette, time_to_threshold) {
  require(survival)
  
  # Get unique R and k combinations
  r_k_combos <- data %>%
    select(R, k) %>%
    unique()
  
  # Create plots for each R/k combination
  plots <- lapply(1:nrow(r_k_combos), function(i) {
    current_R <- r_k_combos$R[i]
    current_k <- r_k_combos$k[i]
    
    # Get data for this R/k combination
    group_data <- data %>%
      filter(R == current_R, k == current_k)
    
    # Create survival data from pre-computed fits
    surv_data <- NULL
    
    for(j in 1:nrow(group_data)) {
      # Extract pre-computed survival fit
      fit <- group_data$surv_fit[[j]]
      scenario <- group_data$scenario[j]
      
      # Create data frame for this scenario and flip the probabilities
      new_data <- data.frame(
        time = fit$time,
        surv = 1 - fit$surv,  # Flip to show probability of reaching threshold
        lower = 1 - fit$upper,  # Note: upper/lower bounds are swapped when flipping
        upper = 1 - fit$lower,
        scenario = scenario,
        R = current_R,
        k = current_k,
        cumulative_cases_threshold = unique(group_data$cumulative_cases_threshold)
      )
      
      # Add to surv_data
      surv_data <- if(is.null(surv_data)) new_data else rbind(surv_data, new_data)
    }
    
    # Create basic survival plot
    p <- ggplot(surv_data, aes(x = time, y = surv, color = scenario)) +
      geom_step() +
      geom_step(aes(y = lower), linetype = "dashed", alpha = 0.5) +
      geom_step(aes(y = upper), linetype = "dashed", alpha = 0.5) +
      scale_color_brewer(palette = color_palette) +
      labs(
        x = "Time (days)",
        y = "Probability of Reaching Threshold",
        title = sprintf(
          "Time to %d cases (R=%.1f, k=%.1f)",
          unique(surv_data$cumulative_cases_threshold),
          unique(surv_data$R),
          unique(surv_data$k)
        )
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    return(list(plot = p))
  })
  
  return(plots)
}

# Function to calculate extinction probability
calculate_extinction_probability <- function(data, initial_flight_chains) {
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(will_fly == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  # Calculate daily cases with complete time series
  data %>%
    filter(destination_infection) %>%
    mutate(day = floor(time)) %>%
    group_by(scenario, R, k, chain, day) %>%
    summarise(daily_cases = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario, R, k),
      chain = flying_chains,
      fill = list(daily_cases = 0)
    ) %>%
    # Calculate cumulative cases for each chain
    group_by(scenario, R, k, chain) %>%
    arrange(day) %>%
    mutate(cum_cases = cumsum(daily_cases)) %>%
    # Calculate proportion of chains with zero cumulative cases each day
    group_by(scenario, R, k, day) %>%
    summarise(
      extinct = sum(cum_cases == 0),
      n = n(),
      .groups = "drop")
}

# Function to plot outbreak probability
plot_outbreak_probability <- function(data, color_palette) {
  ggplot(data, aes(x = day, y = 1 - (extinct / n), colour = scenario)) +
    geom_line() +
    facet_grid(R ~ k) +
    labs(x = "Day", 
         y = "Probability", 
         title = "Cumulative Probability of Outbreak (1 - Extinction Probability)") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to calculate time to reach outbreak probability threshold
calculate_time_to_outbreak_prob <- function(data, prob_threshold = 0.8) {
  # Find first day reaching threshold for each scenario
  data %>%
    mutate(outbreak_prob = 1 - (extinct / n)) %>%
    filter(outbreak_prob >= prob_threshold) %>%
    group_by(scenario, R, k) %>%
    arrange(day) %>%
    slice(1) %>%
    summarise(
      days_to_threshold = first(day),
      outbreak_prob = first(outbreak_prob),
      n_chains = first(n),
      n_outbreaks = first(n - extinct),
      .groups = "drop"
    ) %>%
    mutate(
      prob_formatted = sprintf("%.1f%%", outbreak_prob * 100),
      proportion_outbreaks = n_outbreaks / n_chains
    )
}
