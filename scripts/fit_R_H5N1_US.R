# Code adapted from https://github.com/epiverse-trace/howto/blob/h5_example/analyses/quantify_transmission/reproduction_number_cluster_size.qmd

if(!require("pak")) install.packages("pak")
pak::pak("epiverse-trace/epiparameter")

# Load required packages
library(epichains)
library(MCMCpack)
library(epiparameter)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)

# Get prior for R based on Aditama et al, PLOS ONE, 2012
get_prior <- extract_param(
  type = "percentiles",
  values = c(0.009, 0.315),
  distribution = "gamma",
  percentiles = c(0.025,0.975)
)


h5_prior_r <- function(x){dgamma(x,shape = get_prior[["shape"]], scale = get_prior[["scale"]])}

# Define scenarios
scenario_1 <- c(rep(1, 67), c(2,3)) #67 single spillover, 3x cluster of 2 (cali x2, Miss and sources)
scenario_2 <- c(rep(1, 67), c(2,2), 3) #67 single spillover , 1x cluster of 3 (Miss, source and household contact), 2x cluster of 2 (Cali cases and source)

# Function to perform MCMC and summarize results
run_mcmc <- function(h5_clusters, n_iter = 1e4, n_burn = 1e3) {
  # Likelihood function (as before)
  lik_function <- function(param) {
    if (any(param <= 0)) return(-Inf)
    r_val <- as.numeric(param[1])
    k_val <- as.numeric(param[2])
    log_likelihood <- likelihood(
      chains = h5_clusters,
      statistic = "size",
      offspring_dist = rnbinom,
      size = k_val,
      mu = r_val
    )
    log_prior <- h5_prior_r(r_val)
    return(log_likelihood + log_prior)
  }
  
  # Run MCMC
  result <- MCMCmetrop1R(
    lik_function,
    theta.init = c(R = 0.1, k = 1.0),
    burnin = n_burn,
    mcmc = n_iter,
    thin = 1
  )
  
  # Summarise posterior
  posterior_samples_R <- result[, 1]
  summary_stats <- list(
    mean = mean(posterior_samples_R),
    median = median(posterior_samples_R),
    ci95 = quantile(posterior_samples_R, c(0.025, 0.975))
  )
  
  return(list(result = result, posterior_samples_R = posterior_samples_R, summary_stats = summary_stats))
}

# Run for both scenarios
mcmc_scenario_1 <- run_mcmc(scenario_1)
mcmc_scenario_2 <- run_mcmc(scenario_2)


# Create prior density data
prior_density <- data.frame(
  x = seq(0, 5, length.out = 1000),  # Adjust range as needed
  y = h5_prior_r(seq(0, 5, length.out = 1000)),
  Scenario = "Prior"
)

# Combine posterior densities with prior density
posterior_density_1 <- density(mcmc_scenario_1$posterior_samples_R)
posterior_density_2 <- density(mcmc_scenario_2$posterior_samples_R)

density_df <- rbind(
  data.frame(x = posterior_density_1$x, y = posterior_density_1$y, Scenario = "Scenario 1"),
  data.frame(x = posterior_density_2$x, y = posterior_density_2$y, Scenario = "Scenario 2"),
  prior_density
)

# Plot prior and posterior distributions
R0_US <- ggplot(density_df, aes(x = x, y = y, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1) +
  labs(
    title = "",
    x = "R",
    y = "Density",
    color = "Distribution",
    linetype = "Distribution"
  ) +
  scale_color_manual(values = c("Scenario 1" = "blue", "Scenario 2" = "darkgreen", "Prior" = "red")) +
  scale_linetype_manual(values = c("Scenario 1" = "solid", "Scenario 2" = "solid", "Prior" = "dashed")) +
  scale_x_continuous(
    limits = c(0, 0.15),
    expand = c(0, 0)
  )+
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"), 
    axis.text.x = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "top"
  )

# Save plot 
ggsave(
  filename = file.path("plots", "R0_US.png"), 
  plot = R0_US, 
  width = 8, 
  height = 6, 
  dpi = 300
)

# Summarise results for both scenarios
summary_df <- data.frame(
  Scenario = c("Scenario 1", "Scenario 2"),
  Mean_R = c(mcmc_scenario_1$summary_stats$mean, mcmc_scenario_2$summary_stats$mean),
  Median_R = c(mcmc_scenario_1$summary_stats$median, mcmc_scenario_2$summary_stats$median),
  CI95 = c(
    paste0("(", signif(mcmc_scenario_1$summary_stats$ci95[1], 3), ", ",
           signif(mcmc_scenario_1$summary_stats$ci95[2], 3), ")"),
    paste0("(", signif(mcmc_scenario_2$summary_stats$ci95[1], 3), ", ",
           signif(mcmc_scenario_2$summary_stats$ci95[2], 3), ")")
  )
)

# Display summary table
knitr::kable(summary_df, caption = "Summary of Posterior Estimates for Both Scenarios")


#Dispersion summary 
posterior_samples_k_scenario_1 <- mcmc_scenario_1$result[, 2]
posterior_samples_k_scenario_2 <- mcmc_scenario_2$result[, 2]

summary_k_scenario_1 <- list(
  mean = mean(posterior_samples_k_scenario_1),
  median = median(posterior_samples_k_scenario_1),
  ci95 = quantile(posterior_samples_k_scenario_1, c(0.025, 0.975))
)

summary_k_scenario_2 <- list(
  mean = mean(posterior_samples_k_scenario_2),
  median = median(posterior_samples_k_scenario_2),
  ci95 = quantile(posterior_samples_k_scenario_2, c(0.025, 0.975))
)

results_table <- data.frame(
  Scenario = c("Scenario 1", "Scenario 2"),
  Mean_R = c(mcmc_scenario_1$summary_stats$mean, mcmc_scenario_2$summary_stats$mean),
  Median_R = c(mcmc_scenario_1$summary_stats$median, mcmc_scenario_2$summary_stats$median),
  CI95_R = c(
    paste0("(", signif(mcmc_scenario_1$summary_stats$ci95[1], 3), ", ", 
           signif(mcmc_scenario_1$summary_stats$ci95[2], 3), ")"),
    paste0("(", signif(mcmc_scenario_2$summary_stats$ci95[1], 3), ", ", 
           signif(mcmc_scenario_2$summary_stats$ci95[2], 3), ")")
  ),
  Mean_k = c(summary_k_scenario_1$mean, summary_k_scenario_2$mean),
  Median_k = c(summary_k_scenario_1$median, summary_k_scenario_2$median),
  CI95_k = c(
    paste0("(", signif(summary_k_scenario_1$ci95[1], 3), ", ", 
           signif(summary_k_scenario_1$ci95[2], 3), ")"),
    paste0("(", signif(summary_k_scenario_2$ci95[1], 3), ", ", 
           signif(summary_k_scenario_2$ci95[2], 3), ")")
  )
)

knitr::kable(results_table, caption = "Posterior Estimates for R and k Across Scenarios")


# Calculate Effective Sample Size  for both parameters
ess_scenario_1 <- effectiveSize(mcmc_scenario_1$result)
ess_scenario_2 <- effectiveSize(mcmc_scenario_2$result)

# Create ESS summary table
ess_table <- data.frame(
  Scenario = c("Scenario 1", "Scenario 2"),
  ESS_R = c(ess_scenario_1[1], ess_scenario_2[1]),
  ESS_k = c(ess_scenario_1[2], ess_scenario_2[2])
)

knitr::kable(ess_table, caption = "Effective Sample Size (ESS) for R and k in Each Scenario")


#MCMC race plots

mcmc_obj_1 <- as.mcmc(mcmc_scenario_1$result)
mcmc_obj_2 <- as.mcmc(mcmc_scenario_2$result)


# Basic trace plots 
par(mfrow = c(2, 2))

traceplot(mcmc_obj_1[, 1], main = "Scenario 1 - R")
traceplot(mcmc_obj_1[, 2], main = "Scenario 1 - k")
traceplot(mcmc_obj_2[, 1], main = "Scenario 2 - R")
traceplot(mcmc_obj_2[, 2], main = "Scenario 2 - k")



# Combine MCMC samples into long format
make_trace_df <- function(samples, scenario_label) {
  df <- as.data.frame(samples)
  df$Iteration <- 1:nrow(df)
  df$Scenario <- scenario_label
  
  pivot_longer(
    df,
    cols = -c(Iteration, Scenario),  # exclude these from pivoting
    names_to = "Parameter",
    values_to = "Value"
  )
}

trace_df <- bind_rows(
  make_trace_df(mcmc_scenario_1$result, "Scenario 1"),
  make_trace_df(mcmc_scenario_2$result, "Scenario 2")
)

# ggplot
ggplot(trace_df, aes(x = Iteration, y = Value, color = Scenario)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Parameter + Scenario, scales = "free_y") +
  labs(title = "Trace Plots for R and k", y = "Value") +
  theme_classic() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )
