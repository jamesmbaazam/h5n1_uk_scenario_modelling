#############################
# Generate Primary Plots
#############################

# Create base plots for different metrics
plots <- list(
  daily = plot_daily_cases(all_results, initial_flight_chains, "Set1"),
  cumulative = plot_cumulative_cases(all_results, initial_flight_chains, "Set1"),
  avg_cumulative = plot_avg_cumulative_cases(all_results, initial_flight_chains, "Set1")
)

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
      filename = file.path("results", paste0(spec$name, ".", ext)),
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
p_outbreak_prob <- plot_outbreak_probability(outbreak_prob, "Set1")

# Calculate time to reach case thresholds
case_thresholds <- list(
  list(n = 10, name = "time_to_10_cases"),
  list(n = 100, name = "time_to_100_cases")
)

# Calculate and plot for each threshold
threshold_results <- lapply(case_thresholds, function(threshold) {
  results_data <- calculate_time_to_n_cases(all_results, initial_flight_chains, n_cases = threshold$n)
  plots <- plot_time_to_n_cases(
    data = results_data$results, 
    color_palette = "Set1", 
    time_to_threshold = results_data$time_to_threshold
  )
  
  # Save survival plots quietly
  if(length(plots) > 0) {
    for(i in seq_along(plots)) {
      if(!is.null(plots[[i]]$plot)) {
        plot_name <- sprintf("results/%s_survival_R%.1f_k%.1f", 
                           threshold$name,
                           results_data$results$R[i],
                           results_data$results$k[i])
        
        # Save plots without messages
        suppressMessages({
          ggsave(paste0(plot_name, ".png"), plots[[i]]$plot,
                 width = 12, height = 10, bg = "white", dpi = 600)
          ggsave(paste0(plot_name, ".pdf"), plots[[i]]$plot,
                 width = 12, height = 10, bg = "white")
        })
      }
    }
  }
  
  return(list(
    results = results_data$results,
    plots = plots
  ))
})

# Save outbreak probability plots
ggsave("results/outbreak_probability.png", p_outbreak_prob, 
       width = 12, height = 8, bg = "white", dpi = 600)
ggsave("results/outbreak_probability.pdf", p_outbreak_prob, 
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
    scenario, 
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

#############################
# Save Results
#############################

# Save CSV files
write.csv(time_to_cases_table, 
          "results/time_to_cases_summary.csv", 
          row.names = FALSE)
write.csv(outbreak_threshold_table, 
          "results/outbreak_threshold_summary.csv", 
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
        save_as_html(ft, path = file.path("results", paste0(filename, ".html")))
        
        # Save as docx
        save_as_docx(ft, path = file.path("results", paste0(filename, ".docx")))
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
} else {
  warning("flextable package not available - tables saved as CSV only")
}
