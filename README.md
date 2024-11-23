<!-- badges: start -->
  [![R-CMD-check](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/AnestisTouloumis/geer/graph/badge.svg)](https://app.codecov.io/gh/AnestisTouloumis/geer)
  <!-- badges: end -->  
  

## Installation

You can install the development version of `geer`:

``` r
# install.packages('devtools')
devtools::install_github("AnestisTouloumis/geer")
```

The source code for the development version of `geer` is available on
github at:

- <https://github.com/AnestisTouloumis/geer>

To use `geer`, you should load the package as follows:

``` r
library("geer")
#> Loading required package: gnm
```

## Usage

This package provides a generalized estimating equations (GEE) solver
for fitting marginal regression models with or without an adjustment 
vector.

There are two core functions to fit GEE models:

- `geewa` for fitting GEE models with correlated responses. Options for
  estimation process include the ordinary GEE, bias-reduced or -corrected
  GEE and penalized GEE,
- `geewa_binary` for fitting GEE models with correlated binary responses. Options for
  estimation process include the ordinary GEE, bias-reduced or -corrected
  GEE and penalized GEE.

The main arguments in both functions are:

- an optional data frame (`data`),
- a model formula (`formula`),
- a cluster identifier variable (`id`),
- an optional vector that identifies the order of the observations
  within each cluster (`repeated`).


There are also five useful utility functions:

- `confint` for obtaining Wald–type confidence intervals for the
  regression parameters using the standard errors of the sandwich
  or of the bias-corrected or of the model–based covariance matrix.
  The default option is the sandwich covariance matrix,
- `waldts` for assessing the goodness-of-fit of two nested GEE models
  based on a Wald test statistic,
- `score_test` for assessing the goodness-of-fit of two nested GEE models
  based on a score test statistic,
- `vcov` for obtaining the sandwich, bias-corrected or model–based
  covariance matrix of the regression parameters,
- `gee_criteria` for reporting commonly used criteria to select
  variables and/or association structure for GEE models.
