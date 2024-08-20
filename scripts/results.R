# Calculate relative risk
no_intervention_R <- results_summary$R[results_summary$Scenario == "No testing"]
results_summary$RelativeRisk <- results_summary$R / no_intervention_R

cat("Summary of results:\n")
print(results_summary)

# Visualise the results
# Plot R values
p1 <- ggplot(results_summary, aes(x = Scenario, y = R)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = sprintf("%.2f", R)), vjust = -0.5) +
  ylim(0, max(results_summary$R, na.rm = TRUE) * 1.1) +
  ggtitle("Estimated R for Exported Infections by Testing Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot k values
p2 <- ggplot(results_summary, aes(x = Scenario, y = k)) +
  geom_bar(stat = "identity", fill = "darkred") +
  geom_text(aes(label = sprintf("%.2f", k)), vjust = -0.5) +
  ylim(0, max(results_summary$k, na.rm = TRUE) * 1.1) +
  ggtitle("Estimated k for Exported Infections by Testing Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add a plot to visualise the sensitivity function
t_values <- seq(0, 7, by = 0.1)
sensitivity_values <- sapply(t_values, sensitivity_function)

p3 <- ggplot(data.frame(t = t_values, sensitivity = sensitivity_values), aes(x = t, y = sensitivity)) +
  geom_line() +
  labs(title = "Time-varying Test Sensitivity", x = "Days since infection", y = "Test Sensitivity") +
  theme_minimal()

# Plot Relative Risk
p4 <- ggplot(results_summary, aes(x = Scenario, y = RelativeRisk)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f", RelativeRisk)), vjust = -0.5) +
  ylim(0, max(results_summary$RelativeRisk, na.rm = TRUE) * 1.1) +
  ggtitle("Relative Risk of Exported Infections by Testing Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")

# Arrange plots
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
