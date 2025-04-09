
##### Plotting functions #####
# Function to create daily cases plot
plot_daily_cases <- function(data, initial_flight_chains, color_palette) {

  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  data %>% 
    filter(prevented_flight == FALSE) %>%
    #filter(destination_infection) %>% 
    mutate(day = floor(time)) %>%
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario_name, R, k),
      chain = flying_chains,
      fill = list(n = 0)
    ) %>%
    ggplot(aes(x = day, y = n, color = scenario_name, group = interaction(chain, scenario_name))) +
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
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  data %>% 
    filter(prevented_flight == FALSE) %>%
    mutate(day = floor(time)) %>%
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario_name, R, k),
      chain = flying_chains,
      fill = list(n = 0)
    ) %>%
    group_by(scenario_name, R, k, chain) %>%
    mutate(cumsum_n = cumsum(n)) %>%
    ggplot(aes(x = day, y = cumsum_n, group = interaction(chain, scenario_name), colour = scenario_name)) +
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
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  avg_data <- data %>%
    filter(prevented_flight == FALSE) %>%
    mutate(day = floor(time)) %>%
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(n = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario_name, R, k),
      chain = flying_chains,
      fill = list(n = 0)  # Chains prevented from flying will have zero cases
    ) %>%
    group_by(scenario_name, R, k, chain) %>%
    mutate(cumsum_n = cumsum(n)) %>%
    group_by(scenario_name, R, k, day) %>%
    summarise(
      avg_cumsum = mean(cumsum_n),
      low_ci = quantile(cumsum_n, 0.025),
      high_ci = quantile(cumsum_n, 0.975),
      .groups = "drop"
    )
  
  ggplot(avg_data, aes(x = day, y = avg_cumsum, color = scenario_name)) +
    geom_line() +
    #geom_ribbon(aes(ymin = low_ci, ymax = high_ci, fill = scenario_name), alpha = 0.2, colour = NA) +
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
calculate_time_to_n_cases <- function(data, initial_flight_chains, n_cases) {
  # Convert inputs to data.table for efficiency
  dt <- data.table::as.data.table(data)
  
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  # Calculate cumulative cases by scenario and chain
  cases_by_chain <- dt[
    chain %in% flying_chains,  # Include all potential flying chains
    .(
      cumulative_cases = cumsum(destination_infection),
      time = time,
      prevented = any(prevented_flight)  # Track if flight was prevented
    ), 
    by = .(scenario_name, chain)
  ]
  
  # For prevented flights, set time_to_threshold to Inf (never reaches threshold)
  time_to_threshold <- cases_by_chain[
    prevented == FALSE & cumulative_cases >= n_cases,  # Only look at threshold for non-prevented
    .(time_to_threshold = min(time)), 
    by = .(scenario_name, chain)
  ]
  
  # Add prevented flights with Inf time
  prevented_flights <- unique(cases_by_chain[prevented == TRUE, .(scenario_name, chain)])
  if (nrow(prevented_flights) > 0) {
    prevented_flights[, time_to_threshold := Inf]
    time_to_threshold <- rbind(time_to_threshold, prevented_flights)
  }
  
  # Calculate summary statistics by scenario
  results <- time_to_threshold[, .(
    mean_time = mean(time_to_threshold[is.finite(time_to_threshold)]),  # Only finite times
    median_time = median(time_to_threshold[is.finite(time_to_threshold)]),
    sd_time = sd(time_to_threshold[is.finite(time_to_threshold)]),
    time_to_25pct = quantile(time_to_threshold[is.finite(time_to_threshold)], 0.25),
    time_to_50pct = quantile(time_to_threshold[is.finite(time_to_threshold)], 0.50),
    time_to_75pct = quantile(time_to_threshold[is.finite(time_to_threshold)], 0.75),
    n_chains_reached_threshold = sum(is.finite(time_to_threshold)),
    n_prevented = sum(is.infinite(time_to_threshold))
  ), by = scenario_name]
  
  # Count total chains per scenario
  total_chains <- dt[
    chain %in% flying_chains, 
    .(total_flying_chains = uniqueN(chain)), 
    by = scenario_name
  ]
  
  # Merge results
  results <- merge(results, total_chains, by = "scenario_name")
  
  # Calculate proportions and risk reduction
  results[, ':='(
    proportion_reached = n_chains_reached_threshold / total_flying_chains,
    proportion_prevented = n_prevented / total_flying_chains
  )]
  
  # Calculate delay vs baseline (using first scenario as baseline)
  baseline_time <- results[1, median_time]
  results[, delay_vs_baseline := median_time - baseline_time]
  
  # Calculate risk reduction including prevented flights
  baseline_prop <- results[1, proportion_reached]
  results[, risk_reduction := ifelse(scenario_name == results[1, scenario_name],
                                   0,  # baseline scenario gets 0
                                   100 * (1 - proportion_reached/baseline_prop))]
  
  return(list(
    results = results,
    time_to_threshold = time_to_threshold
  ))
}

# Update the plotting function to use faceting
plot_time_to_n_cases <- function(data, color_palette, time_to_threshold) {
  require(survival)
  require(ggplot2)
  
  # Convert to data.table if not already
  dt <- data.table::as.data.table(time_to_threshold)
  
  # Create survival data for each scenario independently
  surv_data <- purrr::map_dfr(unique(dt$scenario_name), function(scenario) {
    # Get data for this scenario
    scenario_data <- dt[scenario_name == scenario]
    
    # Create survival object for this scenario
    surv_obj <- survival::Surv(
      time = scenario_data$time_to_threshold,
      event = rep(1, nrow(scenario_data))  # All events are observed
    )
    
    # Fit survival curve
    surv_fit <- survival::survfit(surv_obj ~ 1)
    
    # Convert to data frame
    data.frame(
      time = surv_fit$time,
      surv = surv_fit$surv,
      scenario_name = scenario  # Use full scenario name
    )
  }, .id = NULL)  # Don't add an ID column
  
  # Ensure scenario_name is a factor with correct order
  surv_data$scenario_name <- factor(
    surv_data$scenario_name,
    levels = unique(dt$scenario_name)  # Preserve original order
  )
  
  # Create plot
  p <- ggplot(surv_data, aes(x = time, y = surv, color = scenario_name)) +
    geom_step() +
    scale_color_brewer(palette = color_palette) +
    labs(
      x = "Time (days)",
      y = "Proportion not reaching threshold",
      title = "Time to reach case threshold",
      color = "Scenario"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}

# Function to calculate extinction probability
calculate_extinction_probability <- function(data, initial_flight_chains) {
  #browser()
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  # Calculate daily cases with complete time series
  data %>%
    filter(prevented_flight == FALSE) %>%
    mutate(day = floor(time)) %>%
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(daily_cases = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario_name, R, k),
      chain = flying_chains,
      fill = list(daily_cases = 0)
    ) %>%
    # Calculate cumulative cases for each chain
    group_by(scenario_name, R, k, chain) %>%
    arrange(day) %>%
    mutate(cum_cases = cumsum(daily_cases)) %>%
    # Calculate proportion of chains with zero cumulative cases each day
    group_by(scenario_name, R, k, day) %>%
    summarise(
      extinct = sum(cum_cases == 0),
      n = n(),
      .groups = "drop"
    )
}

# Function to plot outbreak probability
plot_outbreak_probability <- function(data, color_palette) {
  ggplot(data, aes(x = day, y = 1 - (extinct / n), colour = scenario_name)) +
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
    group_by(scenario_name, R, k) %>%
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

# Function to calculate average time to reach N cases (naive approach)
calculate_naive_time_to_cases <- function(data, initial_flight_chains, case_thresholds = c(10, 50)) {
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  
  # Calculate time to threshold for each chain
  data %>%
    filter(prevented_flight == FALSE) %>%
    mutate(day = floor(time)) %>%
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(daily_cases = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario_name, R, k),
      chain = flying_chains,
      fill = list(daily_cases = 0)
    ) %>%
    group_by(scenario_name, R, k, chain) %>%
    arrange(day) %>%
    mutate(cum_cases = cumsum(daily_cases)) %>%
    group_by(scenario_name, R, k, chain) %>%
    summarise(
      time_to_cases = map_dbl(case_thresholds, ~{
        if(max(cum_cases) >= .x) {
          min(day[cum_cases >= .x])
        } else {
          NA_real_
        }
      }),
      threshold = case_thresholds,
      .groups = "drop"
    ) %>%
    group_by(scenario_name, R, k, threshold) %>%
    summarise(
      n_reached = sum(!is.na(time_to_cases)),
      mean_time = mean(time_to_cases, na.rm = TRUE),
      median_time = median(time_to_cases, na.rm = TRUE),
      lower_ci = quantile(time_to_cases, 0.025, na.rm = TRUE),
      upper_ci = quantile(time_to_cases, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Add proportion formatting
    mutate(
      total_chains = n(),
      prop_reached = n_reached / total_chains,
      prop_reached_pct = sprintf("%.1f%%", prop_reached * 100)
    )
}

# Function to plot naive time to cases
plot_naive_time_to_cases <- function(data, color_palette) {
  ggplot(data, aes(x = scenario_name, y = mean_time, color = scenario_name)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    facet_grid(threshold ~ paste("R =", R, ", k =", k)) +
    scale_color_brewer(palette = color_palette) +
    labs(
      x = "Scenario",
      y = "Days to Reach Threshold",
      title = "Time to Reach Case Thresholds",
      subtitle = "Mean with 95% CI"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Function to calculate extinction probability over time
calculate_extinction_over_time <- function(data, initial_flight_chains) {
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  # Get prevented flight information and extract intervention day from scenario name
  prevented_flights <- data %>%
    select(scenario_name, chain, prevented_flight) %>%
    distinct() %>%
    mutate(
      # Extract intervention day from scenario name
      # For "No Testing" scenario, set to Inf
      # For other scenarios, extract the number after "Day "
      interventions_enacted = case_when(
        scenario_name == "A. No Testing" ~ Inf,
        TRUE ~ as.numeric(str_extract(scenario_name, "(?<=Day )[0-9]+"))
      )
    )
  
  # Calculate daily cases with complete time series
  data %>%
    mutate(day = floor(time)) %>%
    # First get daily cases per chain
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(daily_cases = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario_name, R, k),
      chain = flying_chains,
      fill = list(daily_cases = 0)
    ) %>%
    # Join with prevented flight information
    left_join(prevented_flights, by = c("scenario_name", "chain")) %>%
    # Calculate extinction metrics per day
    group_by(scenario_name, R, k, chain, day) %>%
    summarise(
      # A chain is extinct if:
      # - it has no cases that day AND
      # - either it's not prevented OR the intervention hasn't started yet
      is_extinct = daily_cases == 0 & 
                  (!prevented_flight | day < interventions_enacted),
      interventions_enacted = first(interventions_enacted),
      .groups = "drop"
    ) %>%
    # Calculate proportion extinct per scenario/R/k/day
    group_by(scenario_name, R, k, day) %>%
    summarise(
      n_extinct = sum(is_extinct, na.rm = TRUE),
      n_total = n(),
      interventions_enacted = first(interventions_enacted),
      .groups = "drop"
    ) %>%
    mutate(
      extinction_prob = n_extinct / n_total
    )
}

# Function to plot extinction probability
plot_extinction_probability <- function(data, color_palette) {
  ggplot(data, aes(x = day, y = extinction_prob, colour = scenario_name)) +
    geom_line() +
    # Add vertical lines for intervention start times
    #geom_vline(aes(xintercept = interventions_enacted, color = scenario_name),
    #           linetype = "dashed", alpha = 0.5) +
    facet_grid(R ~ k) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_brewer(palette = color_palette) +
    labs(x = "Day", 
         y = "Probability", 
         title = "Daily Extinction Probability",
         subtitle = "Proportion of chains with zero cases on each day\nDashed lines show intervention start times") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Function to calculate detection rates over time
calculate_detection_rates <- function(data, initial_flight_chains) {
  # Get only chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain)
  
  # Calculate daily prevention counts more efficiently
  data %>%
    # First filter to only include flying chains
    filter(chain %in% flying_chains) %>%
    # Get just the columns we need
    select(scenario_name, R, k, chain, prevented_flight, time) %>%
    # Calculate day
    mutate(day = floor(time)) %>%
    # Summarize by day
    group_by(scenario_name, R, k, day) %>%
    summarise(
      n_prevented = sum(prevented_flight, na.rm = TRUE),
      n_chains = n_distinct(chain),
      .groups = "drop"
    ) %>%
    # Complete missing days
    complete(
      day = 0:max(day, na.rm = TRUE),
      nesting(scenario_name, R, k),
      fill = list(n_prevented = 0, n_chains = 0)
    ) %>%
    # Calculate rate
    mutate(
      prevention_rate = if_else(n_chains > 0, n_prevented / n_chains, 0)
    )
}

# Function to plot detection rates
plot_detection_rates <- function(data, color_palette) {
  ggplot(data, aes(x = day, y = prevention_rate, color = scenario_name)) +
    geom_line() +
    facet_grid(R ~ k) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_brewer(palette = color_palette) +
    labs(x = "Day", 
         y = "Prevention Rate", 
         title = "Daily Prevention Rate Over Time",
         subtitle = "Proportion of chains prevented from flying each day") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# Function to calculate overall prevention statistics
calculate_detection_summary <- function(data) {
  data %>%
    group_by(scenario_name, R, k) %>%
    summarise(
      total_chains = sum(n_chains),
      total_prevented = sum(n_prevented),
      mean_daily_rate = mean(prevention_rate),
      .groups = "drop"
    ) %>%
    mutate(
      overall_rate = total_prevented / total_chains,
      overall_rate_pct = scales::percent(overall_rate, accuracy = 0.1),
      mean_rate_pct = scales::percent(mean_daily_rate, accuracy = 0.1)
    )
}

# Function to plot split cumulative cases
plot_split_cumulative_cases <- function(data, initial_flight_chains, color_palette) {
  #browser()
  # Get all chains that would have flown in baseline
  flying_chains <- initial_flight_chains %>%
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  # Calculate daily cases split by flyer/descendant
  split_cases <- data %>%
    filter(
      chain %in% flying_chains
    ) %>%
    mutate(
      day = floor(time),
      case_type = ifelse(potential_flyer, "Imported", "Domestic transmission")
    ) %>%
    # Count cases by day and type
    group_by(scenario_name, day, case_type) %>%
    summarise(
      daily_cases = n(),
      .groups = "drop"
    ) %>%
    # Complete missing days/types with zeros
    complete(
      day = 0:max(day),
      nesting(scenario_name),
      case_type = c("Imported", "Domestic transmission"),
      fill = list(daily_cases = 0)
    ) %>%
    # Calculate cumulative cases
    group_by(scenario_name, case_type) %>%
    arrange(day) %>%
    mutate(cumulative_cases = cumsum(daily_cases)) %>%
    ungroup()
  
  # Create stacked area plot
  ggplot(split_cases, aes(x = day, y = cumulative_cases, fill = case_type)) +
    geom_area(alpha = 0.8, position = "stack") +
    facet_wrap(~scenario_name) +
    scale_fill_brewer(palette = color_palette) +
    labs(
      x = "Day",
      y = "Cumulative Cases",
      title = "Cumulative Cases by Source",
      subtitle = "Split by imported cases and domestic transmission",
      fill = "Source"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white"),
      panel.spacing = unit(1, "lines")
    )
}
