# Load required libraries if not already loaded
library(ggplot2)
library(gridExtra)
library(dplyr)
# Plot daily cases
plot_daily_cases <- function(data, color_palette) {
  ggplot(data, aes(x = time, fill = scenario)) +
    geom_histogram(binwidth = 1, position = "dodge") +
    scale_fill_brewer(palette = color_palette) +
    labs(title = "Daily Cases by Scenario",
         x = "Time (days)",
         y = "Number of Cases") +
    theme_minimal()
}

# Plot cumulative cases
plot_cumulative_cases <- function(data, color_palette) {
  data %>%
    group_by(scenario, time) %>%
    summarise(cases = n(), .groups = "drop") %>%
    arrange(scenario, time) %>%
    group_by(scenario) %>%
    mutate(cumulative_cases = cumsum(cases)) %>%
    ggplot(aes(x = time, y = cumulative_cases, color = scenario)) +
    geom_line() +
    scale_color_brewer(palette = color_palette) +
    labs(title = "Cumulative Cases by Scenario",
         x = "Time (days)",
         y = "Cumulative Cases") +
    theme_minimal()
}

# Plot average cumulative cases
plot_avg_cumulative_cases <- function(data, color_palette) {
  data %>%
    group_by(scenario, time) %>%
    summarise(cases = n(), .groups = "drop") %>%
    arrange(scenario, time) %>%
    group_by(scenario) %>%
    mutate(cumulative_cases = cumsum(cases)) %>%
    group_by(scenario) %>%
    summarise(avg_cases = mean(cumulative_cases)) %>%
    ggplot(aes(x = scenario, y = avg_cases, fill = scenario)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = color_palette) +
    labs(title = "Average Cumulative Cases by Scenario",
         x = "Scenario",
         y = "Average Cumulative Cases") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot time to 100 cases
plot_time_to_100_cases <- function(data, color_palette) {
  ggplot(data, aes(x = scenario, y = time_to_100, fill = scenario)) +
    geom_boxplot() +
    scale_fill_brewer(palette = color_palette) +
    labs(title = "Time to 100 Cases by Scenario",
         x = "Scenario",
         y = "Time (days)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot outbreak probability
plot_outbreak_probability <- function(data, color_palette) {
  ggplot(data, aes(x = scenario, y = probability, fill = scenario)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = color_palette) +
    labs(title = "Outbreak Probability by Scenario",
         x = "Scenario",
         y = "Probability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Calculate time to 100 cases
calculate_time_to_100_cases <- function(data) {
  data %>%
    group_by(scenario, chain) %>%
    arrange(time) %>%
    mutate(cumulative = row_number()) %>%
    filter(cumulative == 100) %>%
    select(scenario, chain, time) %>%
    rename(time_to_100 = time)
}

# Calculate extinction probability
calculate_extinction_probability <- function(data) {
  data %>%
    group_by(scenario) %>%
    summarise(
      total_chains = n_distinct(chain),
      extinct_chains = sum(max(time) < Inf),
      probability = extinct_chains / total_chains
    )
}

# Calculate summary statistics
calculate_summary_statistics <- function(data) {
  data %>%
    group_by(scenario) %>%
    summarise(
      R = mean(will_fly),
      k = var(will_fly) / mean(will_fly)^2
    )
}
# Function to create all plots and save them
create_and_save_plots <- function(all_results) {
  # Create plots using combined results
  p_daily <- plot_daily_cases(all_results, "Set1")
  p_cumulative <- plot_cumulative_cases(all_results, "Set1")
  p_avg_cumulative <- plot_avg_cumulative_cases(all_results, "Set1")
  
  # Calculate additional metrics
  time_to_100_cases <- calculate_time_to_100_cases(all_results)
  outbreak_prob <- calculate_extinction_probability(all_results)
  
  p_time_to_100 <- plot_time_to_100_cases(time_to_100_cases, "Set1")
  p_outbreak_prob <- plot_outbreak_probability(outbreak_prob, "Set1")
  
  # Create results directory if it doesn't exist
  dir.create("results", showWarnings = FALSE)
  
  # Save all plots
  ggsave("results/daily_cases.png", p_daily, width = 12, height = 8, bg = "white", dpi = 600)
  ggsave("results/daily_cases.pdf", p_daily, width = 12, height = 8, bg = "white", dpi = 600)
  
  ggsave("results/cumulative_cases.png", p_cumulative, width = 12, height = 8, bg = "white", dpi = 600)
  ggsave("results/cumulative_cases.pdf", p_cumulative, width = 12, height = 8, bg = "white", dpi = 600)
  
  ggsave("results/avg_cumulative_cases.png", p_avg_cumulative, width = 12, height = 8, bg = "white", dpi = 600)
  ggsave("results/avg_cumulative_cases.pdf", p_avg_cumulative, width = 12, height = 8, bg = "white", dpi = 600)
  
  ggsave("results/time_to_100_cases.png", p_time_to_100, width = 12, height = 8, bg = "white", dpi = 600)
  ggsave("results/time_to_100_cases.pdf", p_time_to_100, width = 12, height = 8, bg = "white", dpi = 600)
  
  ggsave("results/outbreak_probability.png", p_outbreak_prob, width = 12, height = 8, bg = "white", dpi = 600)
  ggsave("results/outbreak_probability.pdf", p_outbreak_prob, width = 12, height = 8, bg = "white", dpi = 600)
  
  # Calculate summary statistics
  results_summary <- calculate_summary_statistics(all_results)
  
  # Create summary plots
  p1 <- ggplot(results_summary, aes(x = Scenario, y = R)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    geom_text(aes(label = sprintf("%.2f", R)), vjust = -0.5) +
    ylim(0, max(results_summary$R, na.rm = TRUE) * 1.1) +
    ggtitle("Estimated R for Exported Infections by Testing Scenario") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(results_summary, aes(x = Scenario, y = k)) +
    geom_bar(stat = "identity", fill = "darkred") +
    geom_text(aes(label = sprintf("%.2f", k)), vjust = -0.5) +
    ylim(0, max(results_summary$k, na.rm = TRUE) * 1.1) +
    ggtitle("Estimated k for Exported Infections by Testing Scenario") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add sensitivity function plot
  t_values <- seq(0, 7, by = 0.1)
  sensitivity_values <- sapply(t_values, sensitivity_function)
  
  p3 <- ggplot(data.frame(t = t_values, sensitivity = sensitivity_values), aes(x = t, y = sensitivity)) +
    geom_line() +
    labs(title = "Time-varying Test Sensitivity", x = "Days since infection", y = "Test Sensitivity") +
    theme_minimal()
  
  # Calculate and plot relative risk
  no_intervention_R <- results_summary$R[results_summary$Scenario == "No testing"]
  results_summary$RelativeRisk <- results_summary$R / no_intervention_R
  
  p4 <- ggplot(results_summary, aes(x = Scenario, y = RelativeRisk)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.2f", RelativeRisk)), vjust = -0.5) +
    ylim(0, max(results_summary$RelativeRisk, na.rm = TRUE) * 1.1) +
    ggtitle("Relative Risk of Exported Infections by Testing Scenario") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red")
  
  # Save summary plot
  summary_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave("results/summary_plots.png", summary_plot, width = 15, height = 12, bg = "white", dpi = 600)
  ggsave("results/summary_plots.pdf", summary_plot, width = 15, height = 12, bg = "white", dpi = 600)
  
  # Return summary statistics
  return(results_summary)
}
