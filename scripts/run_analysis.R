require("here")

# Create processing directory if it doesn't exist
dir.create(here("processing"), showWarnings = FALSE)

###################
# Key parameters  #
###################

# Choose flu type for analysis
FLU_TYPES <- "H1N1" #c("H5N1 HPAI", "H1N1")  # Run both types
SCENARIO_SETS <- c(2)  # Run both scenario sets

###################
# Cache settings  #
###################

# Set overwrite flags for different caching stages
overwrite_initial_chains <- TRUE   # Controls regeneration of initial chains
overwrite_flight_chains <- TRUE   # Controls regeneration of flight chains
overwrite_scenarios <- TRUE       # Controls regeneration of scenario results

# Loop through each combination
for(FLU_TYPE in FLU_TYPES) {
  for(SCENARIO_SET in SCENARIO_SETS) {
    message(sprintf("\nProcessing %s - Scenario Set %d (%s)", 
                   FLU_TYPE, 
                   SCENARIO_SET,
                   if(SCENARIO_SET == 1) "Testing by Time" else "Testing Regimes"))
    
    # Create scenario-specific directories
    scenario_name <- if(SCENARIO_SET == 1) "testing_by_time" else "testing_regimes"
    
    # Set up directories for this analysis
    output_dir <- file.path("output", FLU_TYPE, scenario_name)
    processing_dir <- file.path("processing", FLU_TYPE, scenario_name)
    results_dir <- file.path("results", FLU_TYPE, scenario_name)
    cache_dir <- file.path("cache", FLU_TYPE, scenario_name)
    
    # Create all directories
    for(dir in c(output_dir, processing_dir, results_dir, cache_dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Store directory paths in environment for other scripts to use
    assign("ANALYSIS_OUTPUT_DIR", output_dir, envir = .GlobalEnv)
    assign("ANALYSIS_PROCESSING_DIR", processing_dir, envir = .GlobalEnv)
    assign("ANALYSIS_RESULTS_DIR", results_dir, envir = .GlobalEnv)
    assign("ANALYSIS_CACHE_DIR", cache_dir, envir = .GlobalEnv)
    
    #Source analysis scripts
    source(here("scripts", "load_packages.R"))
    source(here("scripts", "load_params.R"))
    source(here("scripts", "functions.R"))
    source(here("scripts", "define_scenarios.R"))
    source(here("scripts", "testing_model2.R"))
    
    source(here("scripts", "plotting_functions.r"))
    source(here("scripts", "results.R"))
    
    message(sprintf("Completed %s - Scenario Set %d", FLU_TYPE, SCENARIO_SET))
    beepr::beep(sound = 1)
    
    # Clean up large objects
    rm(list = c(
      "initial_flight_chains",
      "all_results",
      "plots",
      "outbreak_prob",
      "extinction_prob",
      "detection_rates",
      "threshold_results"
    ))
    
    # Force garbage collection
    gc()
  }
}

message("\nAll analyses complete!")
beepr::beep(sound = 2)
