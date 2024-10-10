# Install and load required packages
if (!require("pak")) install.packages("pak")
#pak::pak(c("epiverse-trace/epichains", "epiverse-trace/epiparameter", "tidyverse", "truncdist", "MASS", "fitdistrplus", "ggplot2", "gridExtra"))

# Load required libraries
library(tidyverse)
#library(tidytable)
library(epichains)
library(truncdist)
library(epiparameter)
library(MASS)
library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(readxl)
library(tictoc)

select <- dplyr::select
