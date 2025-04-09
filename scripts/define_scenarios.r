#############################
# Define Intervention Scenarios
#############################

# Common parameters across scenarios
base_scenario <- list(
  pre = 0,
  post = 0,
  pre_delay = NA,
  post_delay = NA,
  flight_duration = 0.2,
  quarantine_duration = 0,
  name = "A. No Interventions",
  interventions_enacted = 0
)

# Define intervention start times
intervention_days <- c(0, 25, 50, 75, 100)

# Create scenario helper function
create_scenario <- function(name, pre=0, post=0, quarantine=FALSE, start_day=0) {
  modifiers <- c()
  if(pre == 1) modifiers <- c(modifiers, "Pre")
  if(post == 1) modifiers <- c(modifiers, "Post")
  if(quarantine) modifiers <- c(modifiers, "Q")
  
  scenario_name <- if(length(modifiers) > 0) {
    paste0(name, ". ", paste(modifiers, collapse="+"), " Day ", start_day)
  } else {
    paste0(name, ". No Interventions")
  }
  
  return(list(
    pre = pre,
    post = post,
    pre_delay = if(pre == 1) 3 else NA,
    post_delay = if(post == 1) 1 else NA,
    flight_duration = 0.2,
    quarantine_duration = if(quarantine) 14 else 0,
    name = scenario_name,
    interventions_enacted = start_day
  ))
}

if(SCENARIO_SET == 1) {
  #############################
  # Scenario Set 1: Pre-flight Testing Only
  #############################
  
  # Create no testing baseline scenario
  no_testing <- list(
    pre = 0,
    post = 0,
    pre_delay = NA,
    post_delay = NA,
    flight_duration = 0.2,
    quarantine_duration = 0,
    name = "A. No Testing",
    interventions_enacted = 0
  )
  
  # Create scenarios with only pre-flight testing
  scenarios <- c(
    list(no_testing),  # Add baseline scenario first
    lapply(seq_along(intervention_days), function(i) {
      list(
        pre = 1,
        post = 0,
        pre_delay = 3,
        post_delay = NA,
        flight_duration = 0.2,
        quarantine_duration = 0,
        name = sprintf("%s. Pre-flight Day %d", LETTERS[i+1], intervention_days[i]),
        interventions_enacted = intervention_days[i]
      )
    })
  )
  
} else {
  #############################
  # Scenario Set 2: Comprehensive Comparison
  #############################
  
  # Define scenario combinations
  scenario_types <- list(
    list(pre=0, post=0, quarantine=FALSE, name="No Interventions"),
    list(pre=1, post=0, quarantine=FALSE, name="Pre-flight"),
    list(pre=0, post=1, quarantine=FALSE, name="Post-flight"),
    list(pre=1, post=1, quarantine=FALSE, name="Pre+Post"),
    list(pre=1, post=0, quarantine=TRUE,  name="Pre+Q"),
    list(pre=1, post=1, quarantine=TRUE,  name="Pre+Post+Q")
  )
  
  # Generate scenarios
  scenarios <- list()
  
  # Generate all combinations for day 50
  for(i in seq_along(scenario_types)) {
    type <- scenario_types[[i]]
    scenarios[[i]] <- create_scenario(
      LETTERS[i],
      pre = type$pre,
      post = type$post,
      quarantine = type$quarantine,
      start_day = 50
    )
  }
}

# Add debug messages for testing scenarios
message("\nTesting Scenarios:")
message("----------------")
for (scenario in scenarios) {
  message(sprintf("%s:", scenario$name))
  interventions <- c()
  if (scenario$pre == 1) {
    interventions <- c(interventions, 
                      sprintf("Pre-flight test %d day(s) before flight", scenario$pre_delay))
  }
  if (scenario$post == 1) {
    interventions <- c(interventions, 
                      sprintf("Post-flight test %d day(s) after flight", scenario$post_delay))
  }
  if (scenario$quarantine_duration > 0) {
    interventions <- c(interventions, 
                      sprintf("%d day quarantine", scenario$quarantine_duration))
  }
  if (length(interventions) == 0) {
    message("  - No interventions")
  } else {
    message(sprintf("  - Interventions start on day %d:", scenario$interventions_enacted))
    for (intervention in interventions) {
      message(sprintf("    * %s", intervention))
    }
  }
  message("")
}