
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
#> │   ├── Flu_params.xlsx
#> │   ├── si_param_summary.csv
#> │   ├── si_posteriors.csv
#> │   └── si_raw_data
#> │       ├── index.csv
#> │       └── serial.csv
#> ├── h5n1_uk_scenario_modelling.Rproj
#> ├── plots
#> │   ├── outbreak_length.png
#> │   └── outbreak_size.png
#> ├── posterior_predictive
#> │   └── dt_draws.rds
#> └── scripts
#>     ├── fit_si_distributions.R
#>     ├── functions.R
#>     ├── load_packages.R
#>     ├── load_params.R
#>     ├── origin_country_incidence.R
#>     ├── outbreak_distribution.R
#>     ├── results.R
#>     ├── run_analysis.R
#>     ├── testing_model.R
#>     └── testing_model2.R
```

## Analyses

### Outbreak size distribution

The script to run the outbreak size and length distribution simulation
and generate the plots is in `scripts/outbreak_distribution.R`. The
plots are included in the `plot` folder, including the The
[`outbreak_size` plot](plots/outbreak_size.png) and [`outbreak_length`
plot](plots/outbreak_length.png).
