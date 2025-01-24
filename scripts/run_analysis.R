require("here")

# Set overwrite flags for different caching stages
overwrite_initial_chains <- TRUE   # Controls regeneration of initial chains
overwrite_flight_chains <- TRUE   # Controls regeneration of flight chains
overwrite_scenarios <- TRUE       # Controls regeneration of scenario results

#Source analysis scripts
source(here("scripts","load_packages.R"))
source(here("scripts","load_params.R"))
source(here("scripts","functions.R"))
source(here("scripts","testing_model2.R"))
source(here("scripts","plotting_functions.r"))
source(here("scripts","results.R"))
