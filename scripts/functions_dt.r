library(data.table)

# Generate initial chains function
generate_initial_chains <- function(simulation_params, n_chains, stat_threshold) {
  setDT(simulation_params)
  
  simulation_params[, .(
    initial_chain = list(simulate_chains(
      n_chains = n_chains,
      statistic = "length",
      offspring_dist = function(n, mu) rnbinom(n, mu = 1.5, size = 0.1),
      stat_threshold = stat_threshold,
      generation_time = function(n) sample(x = si_draws$y_rep, size = n, replace = TRUE),
      tf = tmax
    ))
  ), by = R_k_id]
}

# Add ancestry column function
add_ancestry_column <- function(dt) {
  # Convert to data.table if not already
  setDT(dt)
  
  # Create ancestry list for each individual
  dt[, ancestry := {
    # Initialize empty list for ancestry
    anc_list <- vector("list", .N)
    
    # For each individual
    for(i in 1:.N) {
      current_anc <- integer(0)
      current_id <- infector[i]
      
      # Trace back through ancestors
      while(!is.na(current_id)) {
        current_anc <- c(current_anc, current_id)
        current_id <- infector[match(current_id, infectee)]
      }
      
      anc_list[[i]] <- current_anc
    }
    anc_list
  }, by = chain]
  
  return(dt)
}

# Determine which individuals fly and prune chains
prune_flight_chains <- function(chain_data, daily_flight_probability, si_draws) {
  # Convert to data.table if not already
  setDT(chain_data)
  
  # For each individual, determine potential flyers
  chain_data[, `:=`(
    potential_flyer = rbinom(.N, 1, daily_flight_probability) == 1,
    will_fly = FALSE
  )]
  
  # Only the first potential flyer in each chain actually flies
  chain_data[, will_fly := potential_flyer & (cumsum(potential_flyer) == 1), by = chain]
  chain_data[will_fly == TRUE, flight_time := time + sample(si_draws$y_rep, .N, replace = TRUE)]
  chain_data[will_fly == FALSE, flight_time := Inf]
  
  # Sort by chain and time
  setorder(chain_data, chain, time)
  
  # Print flight summary
  flight_summary <- chain_data[, .(
    total_flyers = sum(will_fly),
    total_chains = uniqueN(chain),
    total_individuals = .N
  )]
  
  print(paste0("Chain ", chain_data$chain[1], " flight summary:"))
  print(paste0("- Number of people who flew: ", flight_summary$total_flyers))
  print(paste0("- Total individuals: ", flight_summary$total_individuals))
  
  # Add ancestry column
  chain_data <- add_ancestry_column(chain_data)
  
  # Identify chains with a flyer
  flying_chains <- chain_data[, .SD[any(will_fly)], by = chain]
  
  if(nrow(flying_chains) == 0) return(data.table())
  
  # For each chain with a flyer, identify relevant individuals
  flying_chains[, flyer := infectee[will_fly][1], by = chain]
  
  flying_chains[, `:=`(
    is_flyer = infectee == flyer,
    is_ancestor = sapply(ancestry, function(anc) any(flyer %in% anc, na.rm = TRUE)),
    is_descendant = sapply(infectee, function(x) {
      any(sapply(flying_chains[will_fly == TRUE, ancestry], 
                function(anc) any(x %in% anc, na.rm = TRUE)))
    })
  ), by = chain]
  
  # Keep only relevant individuals
  pruned_results <- flying_chains[is_flyer == TRUE | is_ancestor == TRUE | is_descendant == TRUE]
  
  # Clean up temporary columns
  pruned_results[, c("flyer", "potential_flyer", "is_flyer", "is_ancestor", "is_descendant") := NULL]
  
  # Add destination infection indicator
  pruned_results[, destination_infection := sapply(ancestry, function(anc) {
    any(will_fly[match(anc, infectee)], na.rm = TRUE)
  }), by = chain]
  
  return(pruned_results)
}

# Apply testing scenario
apply_testing_scenario <- function(pruned_chains, test_before_flight, test_after_flight, 
                                 pre_flight_test_delay, post_flight_test_delay, 
                                 flight_duration, quarantine_start_day) {
  # Convert to data.table if not already
  setDT(pruned_chains)
  
  # Initialize columns that will be used later
  pruned_chains[, `:=`(
    flight_end = flight_time + flight_duration,
    apply_quarantine = will_fly & (flight_time >= quarantine_start_day),
    prevented_flight = FALSE,
    pre_flight_test_time = Inf,
    pre_flight_test_result = 0,
    isolated_pre_flight = FALSE,
    post_flight_test_time = Inf,
    post_flight_test_result = 0,
    isolated_post_flight = FALSE
  )]
  
  # Pre-flight testing
  if(test_before_flight == 1) {
    pruned_chains[, pre_flight_test_time := fifelse(apply_quarantine, flight_time - pre_flight_test_delay, Inf)]
    pruned_chains[, time_since_infection_at_pre_test := pmax(0, pre_flight_test_time - time)]
    pruned_chains[, pre_flight_test_sensitivity := sapply(time_since_infection_at_pre_test, sensitivity_function)]
    pruned_chains[, pre_flight_test_result := rbinom(.N, 1, pre_flight_test_sensitivity)]
    pruned_chains[, isolated_pre_flight := pre_flight_test_result == 1]
  }
  
  # Post-flight testing
  if(test_after_flight == 1) {
    pruned_chains[, post_flight_test_time := fifelse(apply_quarantine, flight_end + post_flight_test_delay, Inf)]
    pruned_chains[, time_since_infection_at_post_test := pmax(0, post_flight_test_time - time)]
    pruned_chains[, post_flight_test_sensitivity := sapply(time_since_infection_at_post_test, sensitivity_function)]
    pruned_chains[, post_flight_test_result := rbinom(.N, 1, post_flight_test_sensitivity)]
    pruned_chains[, isolated_post_flight := post_flight_test_result == 1]
  }

  # Apply isolation effects after all testing is complete
  pruned_chains[, prevented_flight := isolated_pre_flight & will_fly]
  pruned_chains[, `:=`(
    will_fly = will_fly & !prevented_flight,
    flight_time = fifelse(prevented_flight, Inf, flight_time),
    flight_end = fifelse(prevented_flight, Inf, flight_end)
  )]
  
  # Determine isolation time
  pruned_chains[, isolation_time := pmin(
    fifelse(isolated_pre_flight, pre_flight_test_time, Inf),
    fifelse(isolated_post_flight, post_flight_test_time, Inf)
  )]
  
  # Identify isolated individuals
  isolated_individuals <- pruned_chains[isolation_time < Inf, .(chain, infectee, isolation_time)]
  
  # Prune chains based on isolated ancestors
  pruned_chains[, should_remove := sapply(seq_len(.N), function(i) {
    any(isolated_individuals$infectee %in% ancestry[[i]] &
        isolated_individuals$isolation_time < time[i] &
        isolated_individuals$chain == chain[i])
  })]
  
  # Filter using data.table syntax
  final_chains <- pruned_chains[should_remove == FALSE]
  
  return(final_chains)
}