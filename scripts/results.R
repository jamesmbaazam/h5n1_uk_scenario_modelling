#############################
# Generate Primary Plots
#############################

# Set plot parameters based on scenario set
if(SCENARIO_SET == 1) {
  plot_title_suffix <- paste(tools::toTitleCase(FLU_TYPE), "Flu - Testing by Time")
  color_scale <- "Set1"
  intervention_lines <- data.frame(
    scenario = c("B. Pre-flight Day 0", "C. Pre-flight Day 25", 
                "D. Pre-flight Day 50", "E. Pre-flight Day 75", 
                "F. Pre-flight Day 100"),
    start_day = c(0, 25, 50, 75, 100)
  )
} else {
  plot_title_suffix <- paste(tools::toTitleCase(FLU_TYPE), "Flu - Testing Regimes")
  color_scale <- "Dark2"
  intervention_lines <- data.frame(
    scenario = unique(all_results$scenario_name),
    start_day = rep(50, 6)  # All interventions start at day 50
  )
}

# Create base plots for different metrics
plots <- list(
  daily = plot_daily_cases(all_results, initial_flight_chains, color_scale),
  cumulative = plot_cumulative_cases(all_results, initial_flight_chains, color_scale),
  avg_cumulative = plot_avg_cumulative_cases(all_results, initial_flight_chains, color_scale)
)

# Update plot titles to include analysis type
plots <- lapply(plots, function(p) {
  p + labs(title = paste(tools::toTitleCase(FLU_TYPE), "Flu -", p$labels$title),
           subtitle = paste("Results from", if(SCENARIO_SET == 1) "Testing by Time" else "Testing Regimes"))
})

# Save primary plots in both PNG and PDF formats
plot_specs <- list(
  daily = list(name = "daily_cases", plot = plots$daily),
  cumulative = list(name = "cumulative_cases", plot = plots$cumulative),
  avg_cumulative = list(name = "avg_cumulative_cases", plot = plots$avg_cumulative)
)

# Save each plot with consistent parameters
for (spec in plot_specs) {
  for (ext in c("png", "pdf")) {
    ggsave(
      filename = file.path(ANALYSIS_RESULTS_DIR, paste0(spec$name, ".", ext)),
      plot = spec$plot,
      width = 12, height = 8,
      bg = "white", dpi = 600
    )
  }
}

#############################
# Calculate Additional Metrics
#############################

# Calculate outbreak probability
outbreak_prob <- calculate_extinction_probability(all_results, initial_flight_chains)
p_outbreak_prob <- plot_outbreak_probability(outbreak_prob, color_scale)

# Save outbreak probability plots
ggsave(file.path(ANALYSIS_RESULTS_DIR, "outbreak_probability.png"), p_outbreak_prob, 
       width = 12, height = 8, bg = "white", dpi = 600)
ggsave(file.path(ANALYSIS_RESULTS_DIR, "outbreak_probability.pdf"), p_outbreak_prob, 
       width = 12, height = 8, bg = "white", dpi = 600)

# Calculate time to reach case thresholds
case_thresholds <- list(
  list(n = 10, name = "time_to_10_cases"),
  list(n = 100, name = "time_to_100_cases")
)

# Add debugging before processing thresholds
print("Starting threshold calculations...")
print(paste("Number of scenarios:", length(unique(all_results$scenario_name))))

# Calculate and plot for each threshold using purrr::map
threshold_results <- purrr::map(case_thresholds, function(threshold) {
  print(paste("Processing threshold:", threshold$n, "cases"))
  
  # Calculate results
  results_data <- calculate_time_to_n_cases(
    all_results,
    initial_flight_chains = initial_flight_chains,
    n_cases = threshold$n
  )
  
  # Create plot
  plot_result <- plot_time_to_n_cases(
    data = results_data$results,
    color_palette = color_scale,
    time_to_threshold = results_data$time_to_threshold
  )
  
  # Save survival plot
  plot_name <- file.path(ANALYSIS_RESULTS_DIR, sprintf("%s_survival", threshold$name))
  
  # Save plots
  ggsave(paste0(plot_name, ".png"), plot_result,
         width = 12, height = 10, bg = "white", dpi = 600)
  ggsave(paste0(plot_name, ".pdf"), plot_result,
         width = 12, height = 10, bg = "white")
  
  return(list(
    results = results_data$results,
    plot = plot_result
  ))
})

# Check for failed calculations
failed_thresholds <- purrr::map_lgl(threshold_results, is.null)
if (any(failed_thresholds)) {
  warning("Some threshold calculations failed - check debug output above")
  print(paste("Failed thresholds:", 
              paste(case_thresholds[failed_thresholds], collapse = ", ")))
}

# Calculate and plot extinction probability over time
extinction_prob <- calculate_extinction_over_time(all_results, initial_flight_chains)
p_extinction_prob <- plot_extinction_probability(extinction_prob, color_scale)

# Save extinction probability plots
ggsave(file.path(ANALYSIS_RESULTS_DIR, "extinction_probability.png"), p_extinction_prob, 
       width = 12, height = 8, bg = "white", dpi = 600)
ggsave(file.path(ANALYSIS_RESULTS_DIR, "extinction_probability.pdf"), p_extinction_prob, 
       width = 12, height = 8, bg = "white", dpi = 600)

#############################
# Generate Summary Tables
#############################

# Combine results from different case thresholds
time_to_cases_table <- bind_rows(
  threshold_results[[1]]$results %>% mutate(threshold = "10 cases"),
  threshold_results[[2]]$results %>% mutate(threshold = "100 cases")
) %>%
  select(
    scenario_name, 
    threshold,
    mean_time,
    median_time,
    sd_time,
    time_to_25pct,
    time_to_50pct,
    time_to_75pct,
    delay_vs_baseline,
    risk_reduction,
    n_chains_reached_threshold,
    total_flying_chains,
    proportion_reached
  ) %>%
  rename(
    scenario = scenario_name,
    mean_days = mean_time,
    median_days = median_time,
    sd_days = sd_time,
    days_to_25pct = time_to_25pct,
    days_to_50pct = time_to_50pct,
    days_to_75pct = time_to_75pct,
    delay = delay_vs_baseline,
    risk_reduction_pct = risk_reduction,
    chains_reached = n_chains_reached_threshold,
    total_simulations = total_flying_chains,
    proportion = proportion_reached
  ) %>%
  arrange(threshold, scenario)

# Calculate outbreak threshold table
outbreak_threshold_table <- calculate_time_to_outbreak_prob(outbreak_prob, prob_threshold = 0.8)

# Calculate and plot detection rates
detection_rates <- calculate_detection_rates(all_results, initial_flight_chains)
p_detection_rates <- plot_detection_rates(detection_rates, color_scale) +
  # Add vertical lines for intervention start times
  geom_vline(data = intervention_lines, 
             aes(xintercept = start_day, color = scenario),
             linetype = "dashed", alpha = 0.5)

# Calculate summary statistics
detection_summary <- calculate_detection_summary(detection_rates)

#############################
# Save Results
#############################

# Save CSV files
write.csv(time_to_cases_table, 
          file.path(ANALYSIS_RESULTS_DIR, "time_to_cases_summary.csv"), 
          row.names = FALSE)
write.csv(outbreak_threshold_table, 
          file.path(ANALYSIS_RESULTS_DIR, "outbreak_threshold_summary.csv"), 
          row.names = FALSE)
write.csv(detection_summary, 
          file.path(ANALYSIS_RESULTS_DIR, "detection_summary.csv"), 
          row.names = FALSE)

# Create and save formatted tables with flextable
if (require(flextable)) {
  # Function to create and save flextable
  save_flex_table <- function(data, caption, filename) {
    ft <- flextable(data) %>%
      theme_vanilla() %>%
      set_caption(caption) %>%
      autofit()
    
    # Try to save as HTML and docx only if pandoc is available
    if (rmarkdown::pandoc_available()) {
      tryCatch({
        # Save as HTML
        save_as_html(ft, path = file.path(ANALYSIS_RESULTS_DIR, paste0(filename, ".html")))
        
        # Save as docx
        save_as_docx(ft, path = file.path(ANALYSIS_RESULTS_DIR, paste0(filename, ".docx")))
      }, error = function(e) {
        warning("Could not save formatted tables: ", e$message)
      })
    } else {
      warning("Pandoc not available - skipping HTML and Word output. Tables saved as CSV only.")
    }
    
    return(ft)
  }
  
  # Save tables
  time_to_cases_ft <- save_flex_table(
    time_to_cases_table,
    "Time to Reach Case Thresholds by Scenario",
    "time_to_cases_summary"
  )
  
  outbreak_threshold_ft <- save_flex_table(
    outbreak_threshold_table,
    "Days Until 80% Outbreak Probability by Scenario",
    "outbreak_threshold_summary"
  )
  
  detection_summary_ft <- save_flex_table(
    detection_summary,
    "Overall Case Detection Rates by Scenario",
    "detection_summary"
  )
} else {
  warning("flextable package not available - tables saved as CSV only")
}

# Save detection rate plots
ggsave(file.path(ANALYSIS_RESULTS_DIR, "detection_rates.png"), p_detection_rates, 
       width = 12, height = 8, bg = "white", dpi = 600)
ggsave(file.path(ANALYSIS_RESULTS_DIR, "detection_rates.pdf"), p_detection_rates, 
       width = 12, height = 8, bg = "white", dpi = 600)

# Create split cumulative cases plot
p_split_cases <- plot_split_cumulative_cases(all_results, initial_flight_chains, color_scale)

# Save the plot
ggsave(file.path(ANALYSIS_RESULTS_DIR, "split_cumulative_cases.png"), p_split_cases,
       width = 12, height = 8, bg = "white", dpi = 600)
ggsave(file.path(ANALYSIS_RESULTS_DIR, "split_cumulative_cases.pdf"), p_split_cases,
       width = 12, height = 8, bg = "white")
