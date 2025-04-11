
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modelling scenarios of importation of H5N1 in the UK and the impact of interventions

<!-- badges: start -->

<!-- badges: end -->

This repository contains the code to reproduce the analyses of: [Ward et
al. (2024) Estimates of epidemiological parameters for H5N1 influenza in
humans: a rapid review.
https://doi.org/10.1101/2024.12.11.24318702.](https://doi.org/10.1101/2024.12.11.24318702).

## Package structure

The package is organised into analysis scripts in the `scripts/` folder,
the plots included in the paper are in the `plots/` folder.

``` r
library(fs)
fs::dir_tree()
#> .
#> ├── LICENSE
#> ├── R
#> │   └── utils.R
#> ├── README.Rmd
#> ├── README.md
#> ├── data
#> │   ├── H5N1pptdat.xlsx
#> │   ├── si_param_summary.csv
#> │   ├── si_posteriors.csv
#> │   └── si_raw_data
#> │       ├── index.csv
#> │       └── serial.csv
#> ├── h5n1_uk_scenario_modelling.Rproj
#> ├── plots
#> │   ├── CFR_review.png
#> │   ├── H7N7_R.png
#> │   ├── R0_US.png
#> │   ├── R0_review.png
#> │   ├── inc_review.png
#> │   ├── lat_inf_review.png
#> │   ├── outbreak_length.png
#> │   ├── outbreak_size.png
#> │   ├── serial_review.png
#> │   └── sero_review.png
#> ├── posterior_predictive
#> │   └── dt_draws.rds
#> └── scripts
#>     ├── fit_R_H5N1_US.R
#>     ├── fit_R_H7N7.R
#>     ├── fit_si_distributions.R
#>     ├── outbreak_distribution.R
#>     └── rapid_review_forest_plots.R
```

## Analyses

### H5N1 Epidemiological Parameter Rapid Review

The data for the epidemiological parameters collected from the rapid
review are stored in the `data/` folder in the `H5N1pptdat.xlsx` file.

These parameters are used to produce the forest plots that can be found
in the `plots/` folder. These plots have the name `*_review.png`,
e.g. `inc_review.png` for the incubation periods.

The R script to produce all of the forest plots is the
`rapid_review_forest_plots.R` file in the `scripts/` folder.

### Reproduction Number Estimation

The scripts to estimate the posterior distribution H5N1 data and H7N7
data are in `fit_R_H5N1_US.R` and `fit_R_H7N7.R`, respectively. Both in
the `scripts/` folder.

The posterior distribution plots for the reproduction number estimates
are both in `plots/`: `R0_US.png` and `H7N7_R.png`.

### Outbreak size distribution

The script to run the outbreak size and length distribution simulation
and generate the plots is in `scripts/outbreak_distribution.R`. The
plots are included in the `plot` folder, including the The
[`outbreak_size` plot](plots/outbreak_size.png) and [`outbreak_length`
plot](plots/outbreak_length.png).

## Other branches

This repository contains other branches which contain code that was
written in the process of developing and writing [this
manuscript](https://doi.org/10.1101/2024.12.11.24318702). Here we
briefly describe these branches:

- `bp_report`: This branch contains a report that was written to develop
  the outbreak distribution analysis and plots while the number of cases
  reported in the US was increasing during the second half of 2024. It
  was used as a basis for the outbreak distribution analysis included in
  the paper but is not itself included as part of the publication. It is
  not actively updated and the
  [paper](https://doi.org/10.1101/2024.12.11.24318702) should be
  referred to for the most up-to-date information.

- `deprecated-models`: This branch includes scripts and functions to run
  a travel testing model. This was developed while the project was being
  formulated, but was not included in the final version of the
  manuscript. The code remains on this branch, to potentially be used
  for a future analysis.
