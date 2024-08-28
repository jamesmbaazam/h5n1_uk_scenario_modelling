library(tidyverse)
library(purrr)
library(gridExtra)
library(fitdistrplus)

# Set seed for reproducibility
set.seed(123)

#' Calculate proportion of runs that have controlled outbreak
#'
#' @author Joel Hellewell
#' @export
#' @inheritParams detect_extinct
extinct_prob <- function(outbreak_df_week = NULL, cap_cases  = NULL, week_range = 12:16) {
  
  n_sim <- max(outbreak_df_week$sim)
  
  extinct_runs <- detect_extinct(outbreak_df_week, cap_cases, week_range)
  out <-  sum(extinct_runs$extinct) / n_sim
  
  return(out)
}

#' Calculate whether outbreaks went extinct or not
#' @author Joel Hellewell
#' @param outbreak_df_week data.table  weekly cases produced by the outbreak model
#' @param cap_cases integer number of cumulative cases at which the branching process was terminated
#' @param week_range integer vector giving the (zero indexed) week range to test for whether an extinction occurred.
#' @importFrom data.table as.data.table fifelse
#'
#' @export
#'
detect_extinct <- function(outbreak_df_week  = NULL, cap_cases  = NULL, week_range = 12:16) {
  
  outbreak_df_week <- as.data.table(outbreak_df_week)
  outbreak_df_week <- outbreak_df_week[week %in% week_range]
  outbreak_df_week[, list(
    extinct = fifelse(all(weekly_cases == 0 & cumulative < cap_cases), 1, 0)
  ), by = sim]

}

# Define the flight_test_fun function
flight_test_fun <- function(chain_data, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay, flight_duration) {
  
  results <- chain_data %>%
    mutate(
      flight_time = sample(x = si_draws$y_rep, size = n(), replace = TRUE),
      flight_end = flight_time + flight_duration,
      pre_flight_test_t = Inf,
      pre_flight_test_res = 0,
      post_flight_test_t = Inf,
      post_flight_test_res = 0,
      isolated_pre = FALSE,
      isolated_post = FALSE,
      isolation_time = Inf
    )
  
  # Test before flight
  if (test_before_flight == 1) {
    results <- results %>%
      mutate(
        pre_flight_test_t = flight_time - pre_flight_test_delay,
        pre_flight_test_sensitivity = map_dbl(pre_flight_test_t, sensitivity_function),
        pre_flight_test_res = rbinom(n(), 1, pre_flight_test_sensitivity),
        isolated_pre = if_else(pre_flight_test_res == 1, TRUE, isolated_pre)
      )
  }
  
  # Test after flight
  if (test_after_flight == 1) {
    results <- results %>%
      mutate(
        post_flight_test_t = flight_end + post_flight_test_delay,
        post_flight_test_sensitivity = map_dbl(post_flight_test_t, sensitivity_function),
        post_flight_test_res = rbinom(n(), 1, post_flight_test_sensitivity),
        isolated_post = if_else(post_flight_test_res == 1, TRUE, isolated_post)
      )
  }
  
  # Earliest isolation time
  results <- results %>%
    mutate(
      isolation_time = pmin(if_else(isolated_pre, pre_flight_test_t, Inf), if_else(isolated_post, post_flight_test_t, Inf))
    )
  
  # Indicate whether 2nd+ gen infections were averted if the infector was isolated
  # and if infection were to occur post-flight

  gen_2_averted <- results %>% 
    filter(generation  == 2) %>% 
    select(chain, infector, infectee, generation, flight_time, time) %>%
    right_join(results %>% filter(generation == 1) %>% 
                 select(chain, isolation_time), 
               by = "chain") %>%
    filter(time > isolation_time | time < flight_time)
  
  #if in the chain averted in gen 2, then remove those chains from gen 2 and above
  subsequent_chains <- results %>% 
    filter(generation >= 2) %>% 
    anti_join(gen_2_averted, by = "chain") %>% 
    select(chain, infector, infectee, generation, time, flight_time)
  
  # Combine 1st gen, 2nd gen, and gen 2+ chains
  filtered_chains <- results %>% 
    filter(generation == 1) %>% 
    select(chain, infector, infectee, generation, time, flight_time) %>% 
    bind_rows(subsequent_chains)
  
  return(filtered_chains)
}

# Function to calculate R and k using negative binomial fitting
calculate_r_and_k <- function(chain_data) {
  offspring_counts <- chain_data %>%
    mutate(across(c(chain, generation, infector), as.factor)) %>%
    group_by(chain, .drop = FALSE) %>%
    summarise(N = n(), .groups = "drop") %>%
    pull(N)
  
  if (all(offspring_counts == 0)) {
    warning("All offspring counts are zero. Unable to estimate R and k.")
    return(list(R = NA, k = NA))
  }
  
  tryCatch({
    nb_fit <- fitdistrplus::fitdist(offspring_counts, "nbinom")
    r <- nb_fit$estimate["mu"]
    k <- nb_fit$estimate["size"]
    return(list(R = r, k = k))
  }, error = function(e) {
    warning("Error in fitting negative binomial distribution: ", e$message)
    return(list(R = NA, k = NA))
  })
}

# Function to run a single scenario
run_scenario <- function(sim_chains_df, test_before_flight, test_after_flight, pre_flight_test_delay, post_flight_test_delay,  flight_duration, scenario_name) {
  filtered_results <- flight_test_fun(
    chain_data = sim_chains_df,
    test_before_flight = test_before_flight,
    test_after_flight = test_after_flight,
    pre_flight_test_delay = pre_flight_test_delay,
    post_flight_test_delay = post_flight_test_delay,
    flight_duration = flight_duration
  )
  
  r_k <- calculate_r_and_k(filtered_results)
  
  return(list(filtered_results = filtered_results, 
              R = r_k$R,
              k = r_k$k))
}

# Define scenarios
scenarios <- list(
  list(name = "No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_duration = 0.2),
  list(name = "Pre-flight (1 day before)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2),
  list(name = "Post-flight (1 day after)", pre = 0, post = 1, pre_delay = NA, post_delay = 1, flight_duration = 0.2),
  list(name = "Both (1 day before and after)", pre = 1, post = 1, pre_delay = 1, post_delay = 1, flight_duration = 0.2)
)

# Create a data frame with all combinations of R values and size parameters
simulation_params <- expand_grid(
  R = 2,#seq(0.5, 2, by = 0.1),
  size = 100#c(0.1, 1, 10, 100)
)

run_simulation <- function(R, size) {
  sim_chains <- simulate_chains(
    n_chains = 300,
    statistic = "size",
    offspring_dist = function(n, mu) rnbinom(n, mu = mu, size = size),
    stat_threshold = 10000,
    generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE),
    mu = R
  )
  
  sim_chains_df <- as.data.frame(sim_chains)
  
  results <- map(scenarios, ~run_scenario(sim_chains_df, .$pre, .$post, .$pre_delay, .$post_delay, .$flight_duration, .$name))
  
  results_summary <- tibble(
    Scenario = map_chr(scenarios, "name"),
    results_R = map_dbl(results, "R"),
    results_k = map_dbl(results, "k"),
    Input_R = R,
    Input_Size = size
  )
  
  chains <- tibble(chains = map(results, "filtered_results"),
                   Scenario = map_chr(scenarios, "name"))

  return(list("results_summary" = results_summary,
              "chains" = chains))
}

# Run simulations for all combinations of R and size
all_results <- simulation_params %>%
  mutate(results = map2(R, size, run_simulation)) %>%
  unnest_wider(results)

# plot chains 
plot_chains <- function(chains, Scenario) {
#browser()
  chains %>%
    #aggregate by day
    mutate(time = floor(time)) %>%
    filter(time > flight_time) %>%
    group_by(time) %>% 
    count() %>%
    ggplot(aes(x = time, y = n)) +
          geom_point() +
          geom_line() +
          xlab("Time") +
          ylab("Cases") +
    #name chains by scenario
          ggtitle(Scenario) +
          theme_minimal() +
          theme(legend.position = "bottom")
}

# plot them all together
all_results %>% 
  unnest(chains) %>% 
  mutate(plot = map2(chains, Scenario,  plot_chains)) %>% 
  pull(plot) %>% 
  wrap_plots()&
  coord_cartesian(xlim = c(0, 100))&
  scale_y_continuous(limits = c(0, 10000))
  
# 
# # Create plots
# create_plot_R <- function(data) {
#   ggplot(data, aes(x = Input_R, y = results_R, color = Scenario, group = Scenario)) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~ Input_Size, scales = "free_y", 
#                labeller = labeller(Input_Size = function(x) paste("Size =", x))) +
#     xlab("Input R") +
#     ylab("Estimated R") +
#     ggtitle("Estimated R for Different Testing Scenarios") +
#     theme_minimal() +
#     theme(legend.position = "bottom") +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
# }
# 
# create_plot_k <- function(data) {
#   ggplot(data, aes(x = Input_R, y = results_k, color = Scenario, group = Scenario)) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~ Input_Size, scales = "free_y", 
#                labeller = labeller(Input_Size = function(x) paste("Size =", x))) +
#     xlab("Input R") +
#     ylab("Estimated k") +
#     ggtitle("Estimated k for Different Testing Scenarios") +
#     theme_minimal() +
#     theme(legend.position = "bottom") +
#     geom_hline(aes(yintercept = Input_Size), linetype = "dashed", color = "black")
# }
# 
# # Generate plots
# plot_R <- create_plot_R(all_results)
# plot_k <- create_plot_k(all_results)
# 
# # Display plots
# print(plot_R)
# print(plot_k)
# 
# # Add a plot to visualise the sensitivity function
# t_values <- seq(0, 7, by = 0.1)
# sensitivity_values <- map_dbl(t_values, sensitivity_function)
# 
# p_sensitivity <- tibble(t = t_values, sensitivity = sensitivity_values) %>%
#   ggplot(aes(x = t, y = sensitivity)) +
#   geom_line() +
#   labs(title = "Time-varying Test Sensitivity", x = "Days since infection", y = "Test Sensitivity") +
#   theme_minimal()
# 
# # Display the sensitivity function plot
# print(p_sensitivity)
