
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
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  data %>% 
    filter(prevented_flight == FALSE) %>%
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
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  avg_data <- data %>%
    filter(prevented_flight == FALSE) %>%
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
    filter(potential_flyer == TRUE) %>%
    pull(chain) %>%
    unique()
  
  total_flying_chains <- length(flying_chains)
  
  # Calculate time to threshold for each chain
  time_to_threshold <- data %>%
    filter(destination_infection) %>%
    filter(prevented_flight == FALSE) %>%  # Keep only cases where flight wasn't prevented
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
    # For prevented flights, set status to FALSE
    left_join(
      data %>% 
        select(scenario, chain, prevented_flight) %>% 
        distinct(),
      by = c("scenario", "chain")
    ) %>%
    summarise(
      time = if(max(cum_cases) >= n_cases && !any(prevented_flight, na.rm = TRUE)) {
        as.numeric(min(day[cum_cases >= n_cases]))
      } else {
        as.numeric(max(day))
      },
      status = max(cum_cases) >= n_cases && !any(prevented_flight, na.rm = TRUE),
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
      surv_summary = map(surv_fit, ~summary(.x)),
      # Handle cases where median is NA (no events)
      actual_results = map_dbl(surv_fit, ~{
        med <- quantile(.x, probs = 0.5)$quantile
        if(is.na(med)) max(time_to_threshold$time) else med
      }),
      # Handle confidence intervals using quantile method
      lower_ci = map_dbl(surv_fit, ~{
        q <- quantile(.x, conf.int = TRUE)
        if(length(q$lower) == 0) NA else q$lower[1]
      }),
      upper_ci = map_dbl(surv_fit, ~{
        q <- quantile(.x, conf.int = TRUE)
        if(length(q$upper) == 0) NA else q$upper[1]
      }),
      rmst = map_dbl(surv_summary, ~.x$table["rmean"]),
      rmst_se = map_dbl(surv_summary, ~.x$table["se(rmean)"])
    ) %>%
    mutate(
      # Add time percentiles with NA handling
      time_to_25pct = map_dbl(surv_fit, ~{
        q <- quantile(.x, probs = 0.25)$quantile
        if(length(q) == 0) NA else q
      }),
      time_to_50pct = map_dbl(surv_fit, ~{
        q <- quantile(.x, probs = 0.50)$quantile
        if(length(q) == 0) NA else q
      }),
      time_to_75pct = map_dbl(surv_fit, ~{
        q <- quantile(.x, probs = 0.75)$quantile
        if(length(q) == 0) NA else q
      })
    ) %>%
    group_by(R, k) %>%
    mutate(
      # Calculate comparative statistics within each R/k group
      baseline_result = actual_results[scenario == "B. Pre-flight Day 0"],
      baseline_proportion = proportion_reached[scenario == "B. Pre-flight Day 0"],
      
      delay_vs_baseline = if_else(
        scenario == "B. Pre-flight Day 0",
        0,
        actual_results - baseline_result
      ),
      
      risk_reduction = if_else(
        scenario == "B. Pre-flight Day 0",
        0,
        (baseline_proportion - proportion_reached) / baseline_proportion
      )
    ) %>%
    ungroup() %>%
    select(-baseline_result, -baseline_proportion)
  
  return(list(
    results = results,
    time_to_threshold = time_to_threshold
  ))
}

# Update the plotting function to use faceting
plot_time_to_n_cases <- function(data, color_palette, time_to_threshold) {
  require(survival)
  
  # Create survival data from pre-computed fits for all R/k combinations
  surv_data <- data %>%
    rowwise() %>%
    mutate(
      # Extract survival data for each fit
      surv_points = list(data.frame(
        time = surv_fit[[1]]$time,
        surv = 1 - surv_fit[[1]]$surv,  # Flip to show probability of reaching threshold
        lower = 1 - surv_fit[[1]]$upper, # Note: upper/lower bounds are swapped when flipping
        upper = 1 - surv_fit[[1]]$lower
      ))
    ) %>%
    unnest(surv_points)
  
  # Create faceted survival plot
  ggplot(surv_data, aes(x = time, y = surv, color = scenario)) +
    geom_step() +
   #geom_step(aes(y = lower), linetype = "dashed", alpha = 0.5) +
    #geom_step(aes(y = upper), linetype = "dashed", alpha = 0.5) +
    facet_grid(R ~ k, labeller = label_both) +
    scale_color_brewer(palette = color_palette) +
    labs(
      x = "Time (days)",
      y = "Probability of Reaching Threshold",
      title = sprintf("Time to %d cases", unique(data$cumulative_cases_threshold))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
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
    group_by(scenario, R, k, chain, day) %>%
    summarise(daily_cases = n(), .groups = "drop") %>%
    # Complete with ALL chains that would have flown in baseline
    complete(
      day = 0:max(floor(data$time), na.rm = TRUE),
      nesting(scenario, R, k),
      chain = flying_chains,
      fill = list(daily_cases = 0)
    ) %>%
    group_by(scenario, R, k, chain) %>%
    arrange(day) %>%
    mutate(cum_cases = cumsum(daily_cases)) %>%
    group_by(scenario, R, k, chain) %>%
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
    group_by(scenario, R, k, threshold) %>%
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
  ggplot(data, aes(x = scenario, y = mean_time, color = scenario)) +
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
