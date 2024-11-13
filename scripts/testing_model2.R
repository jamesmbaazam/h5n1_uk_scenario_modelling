library(tidyverse)
library(purrr)
library(gridExtra)
library(fitdistrplus)

# Set seed for reproducibility
set.seed(123)

tmax = 1000
daily_flight_probability = 0.5 # 0.1% chance of flying each day

# Define scenarios
scenarios <- list(
  list(name = "A. No testing", pre = 0, post = 0, pre_delay = NA, post_delay = NA, flight_duration = 0.2, quarantine_start_day = Inf),
  list(name = "B. Pre-flight (1 day before)", pre = 1, post = 0, pre_delay = 1, post_delay = NA, flight_duration = 0.2, quarantine_start_day = 10),
  list(name = "C. Post-flight (1 day after)", pre = 0, post = 1, pre_delay = NA, post_delay = 1, flight_duration = 0.2, quarantine_start_day = 10),
  list(name = "D. Both (1 day before and after)", pre = 1, post = 1, pre_delay = 1, post_delay = 1, flight_duration = 0.2, quarantine_start_day = 10)
)

# Create a data frame with all combinations of R values and size parameters
simulation_params <- expand_grid(
  R = c(2),
  size = c(0.1)
) %>% 
  mutate(R_k_id = row_number())

# Generate initial chains
initial_chains <- generate_initial_chains(simulation_params, n_chains = 50, stat_threshold = 10)

# Run scenarios on the initial chains
tic()
all_results <- initial_chains %>%
  mutate(results = map(initial_chains, ~ run_scenarios_on_chains(., scenarios, daily_flight_probability, si_draws))) %>%
  unnest(results)
toc()



# Plot aggregated daily cases
all_results %>% 
  filter(destination_infection) %>% 
  #filter(time > flight_time) %>%
  mutate(day = floor(time)) %>%
  group_by(scenario, R, size, chain, day) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = day, y = n, color = scenario, group = chain)) +
  geom_line() +
  facet_grid(R~size+scenario) +
  labs(x = "Day", y = "Daily Cases", title = "Daily Cases by Scenario and R_k_id") 

#Plot cumulative cases
all_results %>% 
  #filter(time > flight_time) %>%
  mutate(day = floor(time)) %>%
  group_by(scenario, R, size, chain, day) %>%
  summarise(n = n()) %>%
  group_by(scenario, R, size, chain) %>%
  mutate(cumsum_n = cumsum(n)) %>% 
  ggplot(aes(x = day, y = cumsum_n, group = chain, colour = scenario)) +
  geom_line() +
  facet_grid(R~size+scenario) +
  labs(x = "Day", y = "Cumulative Cases", title = "Cumulative Cases by Scenario and R_k_id")

#Plot individual chains
all_results %>% 
  mutate(day = floor(time)) %>%
  group_by(scenario, R, size, day, chain) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = day, y = n, group = interaction(chain,scenario),colour = scenario)) +
  geom_line(alpha = 0.2) +
  facet_grid(R~size) +
  labs(x = "Day", y = "Daily Cases", title = "Daily Cases by Scenario and R_k_id")

