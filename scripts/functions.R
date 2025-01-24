# Modified generate_initial_chains function
generate_initial_chains <- function(simulation_params, n_chains, stat_threshold, tmax) {
  # Convert to data.table and ensure copy
  dt <- data.table::as.data.table(simulation_params)
  data.table::setDT(dt)  # Ensure it's a data.table
  
  # Create list to store chains
  chain_list <- vector("list", nrow(dt))
  
  # Generate chains for each parameter combination
  for(i in seq_len(nrow(dt))) {
    current_R <- dt$R[i]
    current_k <- dt$k[i]
    
    chain_list[[i]] <- simulate_chains(
      n_chains = n_chains,
      statistic = "length",
      offspring_dist = function(n, mu) {
        rnbinom(n, mu = current_R, size = current_k)
      },
      stat_threshold = stat_threshold,
      generation_time = function(n) sample(x = si_draws, size = n, replace = TRUE),
      tf = tmax
    )
  }
  
  # Add chains to data.table
  dt[, initial_chain := chain_list]
  
  # Convert back to data.frame if needed
  return(as.data.frame(dt))
}

initialize_flight_chains <- function(chain_data, daily_flight_probability, si_draws, tmax, 
                                   scenario = list(flight_duration = 1)) {
  message("Initializing flight chains...")
  dt <- data.table::as.data.table(chain_data)
  n_cases <- nrow(dt)
  message(sprintf("Processing %d cases", n_cases))
  
  # Initial flight assignments
  dt[, potential_flyer := rbinom(n_cases, 1, daily_flight_probability) == 1]
  n_flyers <- sum(dt$potential_flyer)
  message(sprintf("Identified %d potential flyers (%.1f%%)", n_flyers, 100*n_flyers/n_cases))
  
  # Exit if no flights
  if(!any(dt$potential_flyer)) {
    message("No flights found, returning empty dataset")
    return(list(
      data = data.table(
        chain = integer(0),
        infectee = integer(0),
        potential_flyer = logical(0),
        time = numeric(0),
        flight_time = numeric(0),
        flight_end = numeric(0)
      ),
      n_flying_chains = 0
    ))
  }
  
  # Extract potential flyers
  flyers <- dt[potential_flyer == TRUE]
  
  # For each flyer, trace back to find if they have flying ancestors
  setkey(dt, chain, infectee)
  flyers[, is_first_flight := TRUE]
  
  for(i in 1:nrow(flyers)) {
    if(i %% 10 == 0) message(sprintf("Processing flyer %d of %d", i, nrow(flyers)))
    
    current_chain <- flyers[i, chain]
    current_case <- flyers[i, infectee]
    
    # Trace back until we find a flying ancestor or hit root
    while(!is.na(dt[.(current_chain, current_case), infector])) {
      current_case <- dt[.(current_chain, current_case), infector]
      if(dt[.(current_chain, current_case), potential_flyer]) {
        flyers[i, is_first_flight := FALSE]
        break
      }
    }
  }
  
  # Keep only first flyers and set their flight times
  first_flyers <- flyers[is_first_flight == TRUE]
  first_flyers[, flight_delay := sample(si_draws, .N, replace = TRUE)]
  first_flyers[, ':='(
    flight_time = time + flight_delay,
    flight_end = time + flight_delay + scenario$flight_duration
  )]
  
  message(sprintf("Found %d first flyers", nrow(first_flyers)))
  
  # Create index of cases to keep
  keep_idx <- vector("logical", nrow(dt))
  
  # Process each first flyer separately
  for(i in 1:nrow(first_flyers)) {
    current_chain <- first_flyers[i, chain]
    flyer_id <- first_flyers[i, infectee]
    
    # Mark flyer
    keep_idx[dt[.(current_chain, flyer_id), which = TRUE]] <- TRUE
    
    # Mark all descendants
    current_ids <- flyer_id
    while(length(current_ids) > 0) {
      # Find next generation
      descendant_rows <- dt[chain == current_chain & infector %in% current_ids, which = TRUE]
      if(length(descendant_rows) == 0) break
      
      # Mark these descendants
      keep_idx[descendant_rows] <- TRUE
      
      # Move to next generation
      current_ids <- dt[descendant_rows, infectee]
    }
  }
  
  # Keep marked cases
  result <- dt[keep_idx]
  
  # Add flight times from first flyers
  result[, ':='(
    flight_time = NA_real_,
    flight_end = NA_real_,
    potential_flyer = FALSE
  )]
  
  # Update flight times and potential_flyer status for each chain
  for(i in 1:nrow(first_flyers)) {
    current_chain <- first_flyers[i, chain]
    result[chain == current_chain, ':='(
      flight_time = first_flyers[i, flight_time],
      flight_end = first_flyers[i, flight_end],
      potential_flyer = first_flyers[i, potential_flyer]
    )]
  }
  
  # Calculate destination infections
  result[, destination_infection := time > flight_end]
  
  # Sort
  setorder(result, chain, generation, time)
  
  return(result)
}

process_scenario <- function(flying_chains, scenario, quarantine_duration = 14, tmax) {
  message("Processing scenario: ", scenario$name)
  
  # Convert to data.table if not already
  dt <- if(data.table::is.data.table(flying_chains)) {
    data.table::copy(flying_chains)
  } else {
    data.table::as.data.table(flying_chains)
  }
  
  # Print column names for debugging
  message("Available columns: ", paste(names(dt), collapse = ", "))
  
  # Check if required columns exist
  required_cols <- c("chain", "infectee", "potential_flyer", "time", "flight_time", "flight_end")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Set key after ensuring columns exist
  data.table::setkey(dt, chain, infectee)
  
  # Extract flyers and their original flight times
  flyer_info <- dt[potential_flyer == TRUE, .(
    chain,
    infectee,
    flight_time,
    flight_end
  )]
  
  # Set interventions_active and base fields for ALL cases
  dt[, ':='(
    interventions_active = time >= scenario$interventions_enacted,
    original_flight_time = flight_time,
    original_flight_end = flight_end,
    prevented_flight = FALSE  # Initialize for all cases
  )]
  
  # Process flyers for testing
  flyers <- dt[potential_flyer == TRUE]
  
  # Initialize testing columns for flyers
  flyers[, ':='(
    pre_flight_test_time = Inf,
    post_flight_test_time = Inf,
    quarantine_end = Inf,
    pre_flight_test_sensitivity = 0,
    post_flight_test_sensitivity = 0,
    pre_flight_test_result = 0,
    post_flight_test_result = 0,
    time_since_infection_pre = NA_real_,
    time_since_infection_post = NA_real_
  )]
  
  # Pre-flight testing (only for flyers)
  if(scenario$pre == 1) {
    message("Applying pre-flight testing")
    
    # First set test times
    flyers[interventions_active == TRUE, 
           pre_flight_test_time := pmax(0, original_flight_time - scenario$pre_delay)]
    
    # Then calculate time since infection
    flyers[interventions_active == TRUE, 
           time_since_infection_pre := pmax(0, pre_flight_test_time - time)]
    
    # Debug values
    message("\nSample of test times and infections:")
    print(flyers[interventions_active == TRUE, 
                 .(time, pre_flight_test_time, time_since_infection_pre)] %>% 
          head(10))
    
    # Calculate test sensitivity using new function with uncertainty
    flyers[interventions_active == TRUE, 
           pre_flight_test_sensitivity := sensitivity_function(time_since_infection_pre, sample_uncertainty = TRUE)]
    
    # Debug message
    message("Test sensitivity summary:")
    print(summary(flyers[interventions_active == TRUE, pre_flight_test_sensitivity]))
    
    # Generate test results
    flyers[interventions_active == TRUE, ':='(
      pre_flight_test_result = rbinom(.N, 1, 
                                    pmin(1, pmax(0, pre_flight_test_sensitivity)))
    )]
    
    # Debug message
    message("Number of positive tests: ", 
           sum(flyers[interventions_active == TRUE, pre_flight_test_result]))
    
    # Mark prevented flights
    flyers[interventions_active == TRUE & pre_flight_test_result == 1, 
           prevented_flight := TRUE]
    
    # Debug final counts
    message("Number of prevented flights: ",
           sum(flyers[interventions_active == TRUE, prevented_flight]))
  }
  
  # Update main dataset with prevention status
  dt[flyers, on = .(chain, infectee), ':='(
    prevented_flight = i.prevented_flight,
    potential_flyer = i.potential_flyer,
    quarantine_end = i.quarantine_end
  )]
  
  # Keep only flyers (including prevented ones) and their descendants
  keep_idx <- vector("logical", nrow(dt))
  
  # Process each flyer (both prevented and non-prevented)
  message("Processing descendants")
  flyers_to_process <- dt[potential_flyer == TRUE]
  
  for(i in 1:nrow(flyers_to_process)) {
    if(i %% 100 == 0) message(sprintf("Processing flyer descendants (flyer %d of %d)", i, nrow(flyers_to_process)))
    
    current_chain <- flyers_to_process[i, chain]
    flyer_id <- flyers_to_process[i, infectee]
    
    # Mark flyer and set their flight times
    matching_rows <- dt[.(current_chain, flyer_id), which = TRUE]
    keep_idx[matching_rows] <- TRUE
    
    # Only set flight times if not prevented
    if (!flyers_to_process[i, prevented_flight]) {
      dt[matching_rows, ':='(
        flight_time = flyers_to_process[i, flight_time],
        flight_end = flyers_to_process[i, flight_end]
      )]
      
      # Mark all descendants and set their flight times
      current_ids <- flyer_id
      while(length(current_ids) > 0) {
        # Find next generation
        descendant_idx <- dt[chain == current_chain & infector %in% current_ids, which = TRUE]
        if(length(descendant_idx) == 0) break
        
        # Mark these descendants and set their flight times
        keep_idx[descendant_idx] <- TRUE
        dt[descendant_idx, ':='(
          flight_time = flyers_to_process[i, flight_time],
          flight_end = flyers_to_process[i, flight_end]
        )]
        
        # Move to next generation
        current_ids <- dt[descendant_idx, infectee]
      }
    }
  }
  
  # Keep marked cases
  result <- dt[keep_idx]
  
  # Join with flyers to get all testing information
  result[flyers, on = .(chain, infectee), ':='(
    pre_flight_test_time = i.pre_flight_test_time,
    post_flight_test_time = i.post_flight_test_time,
    pre_flight_test_sensitivity = i.pre_flight_test_sensitivity,
    post_flight_test_sensitivity = i.post_flight_test_sensitivity,
    pre_flight_test_result = i.pre_flight_test_result,
    post_flight_test_result = i.post_flight_test_result,
    time_since_infection_pre = i.time_since_infection_pre,
    time_since_infection_post = i.time_since_infection_post
  )]
  
  # Calculate destination infections based on flyer times
  result[, ':='(
    destination_infection = 
      !is.na(flight_end) & 
      time > flight_end & 
      (time > quarantine_end | quarantine_end == Inf),
    scenario = scenario$name
  )]
  
  setorder(result, chain, generation, time)
  return(as.data.frame(result))
}
