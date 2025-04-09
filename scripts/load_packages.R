# Install and load required packages
if (!require("pak")) install.packages("pak")
#pak::pak(c("tidyverse", "tidytable", "epichains", "truncdist", "epiparameter", "MASS", "fitdistrplus", "ggplot2", "gridExtra", "readxl", "flextable", "qs", "survival"))

# Load required libraries
library(tidyverse)
library(tidytable)
library(epichains)
library(truncdist)
library(epiparameter)
library(epiparameter) 
library(MASS)
library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(readxl)
library(tictoc)
library(data.table)
library(patchwork)
library(here)
library(flextable)
library(qs)
library(survival)
library(tictoc)


select <- dplyr::select
