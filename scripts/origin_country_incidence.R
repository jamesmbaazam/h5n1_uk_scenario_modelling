sim_chains <- simulate_chains(
  n_chains = 100,
  statistic = "length",
  offspring_dist = function(n, mu) rnbinom(n, mu = 1.5, size = 0.1),
  stat_threshold = 20,
  generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE)
)

sim_chains_agg <- aggregate(sim_chains)

sim_chains_df %>% 
  mutate(day = floor(