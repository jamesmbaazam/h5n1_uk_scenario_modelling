create_intervals <- function(dt) {
  list(
    zoonotic = rep(dt$days, dt$zoonotic),
    nonzoonotic = rep(dt$days, dt$nonzoonotic)
  )
}

create_stan_data <- function(intervals) {
  lapply(intervals, function(int) {
    list(N = length(int), N_rep = 1000, y = int)
  })
}

fit_models <- function(model, stan_data) {
  lapply(stan_data, function(data) {
    model$sample(data = data, chains = 4, iter_warmup = 1000, iter_sampling = 2000)
  })
}

save_fits <- function(fits, prefix) {
  lapply(names(fits), function(name) {
    saveRDS(fits[[name]], paste0("fits/", prefix, "_", name, ".rds"))
  })
}

extract_draws <- function(fit) {
  spread_draws(fit, y_rep[i]) |> data.table()
}

add_metadata <- function(draws, exposure, case_type, distribution) {
  lapply(draws, function(draw) {
    draw[, `Exposure type` := exposure][, `Case type` := case_type][, Distribution := distribution]
  })
}

combine_draws <- function(case_source, posterior_draws, distribution) {
  exposure_cases <- list(
    "index" = list("Non-Zoonotic", "Index"),
    "serial" = list("Non-Zoonotic", "Serial")
  )
  
  distribution_capitalized <- switch(
    distribution,
    "gamma" = "Gamma",
    "lognormal" = "Lognormal",
    distribution # default case if other distributions are added
  )
  
  result_list <- lapply(c("nonzoonotic", "zoonotic"), function(exposure_type) {
    case_details <- if (exposure_type == "nonzoonotic") {
      list("Non-Zoonotic", exposure_cases[[case_source]][2])
    } else {
      list("Zoonotic", exposure_cases[[case_source]][2])
    }
    
    draw_dt <- posterior_draws[[distribution]][[case_source]][[exposure_type]]
    
    # Initialize and assign all columns explicitly with capitalized terms
    draw_dt[, `Exposure Type` := case_details[[1]]]
    draw_dt[, `Case Source` := case_details[[2]]]
    draw_dt[, `Distribution Type` := distribution_capitalized]
    
    return(draw_dt)
  })
  
  # Combine the list of data.tables into a single data.table
  combined_dt <- rbindlist(result_list, fill = TRUE) # `fill = TRUE` to ensure all columns are present
  
  return(combined_dt)
}