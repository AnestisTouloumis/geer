<!-- badges: start -->
[![R-CMD-check](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/AnestisTouloumis/geer/graph/badge.svg)](https://app.codecov.io/gh/AnestisTouloumis/geer)
<!-- badges: end -->

## Overview

`geer` fits marginal models for independent, repeated, or clustered
responses using Generalized Estimating Equations (GEE). Supported
estimation methods include the traditional GEE, bias-reducing GEE,
bias-corrected GEE, and Jeffreys-prior penalized GEE. Continuous and
count responses are handled by `geewa`, while binary responses are
handled by `geewa_binary` through an odds-ratio parameterization.

## Installation

You can install the development version of `geer` from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("AnestisTouloumis/geer")
```

## Usage

Load the package:

``` r
library("geer")
```

### Fitting models

There are two core fitting functions:

- `geewa()` for continuous and count responses (Gaussian, Poisson,
  binomial, Gamma, inverse Gaussian, and quasi families).
- `geewa_binary()` for binary responses via a marginalized odds-ratio
  parameterization.

Both functions support the following estimation methods via the
`method` argument:

| Method | Description |
|---|---|
| `"gee"` | Traditional GEE |
| `"brgee-naive"`, `"brgee-robust"`, `"brgee-empirical"` | Bias-reducing GEE |
| `"bcgee-naive"`, `"bcgee-robust"`, `"bcgee-empirical"` | Bias-corrected GEE |
| `"pgee-jeffreys"` | Fully iterated Jeffreys-prior penalized GEE |
| `"opgee-jeffreys"` | One-step penalized GEE |
| `"hpgee-jeffreys"` | Hybrid one-step GEE |

Working correlation structures for `geewa()` are controlled by
`corstr`: `"independence"`, `"exchangeable"`, `"ar1"`,
`"m-dependent"`, `"unstructured"`, and `"fixed"`. Working odds-ratio
structures for `geewa_binary()` are controlled by `orstr`:
`"independence"`, `"exchangeable"`, `"unstructured"`, and `"fixed"`.

Convergence and fitting options are set via `geer_control()`.

### Inference

Standard S3 methods are available for fitted `geer` objects:

- `summary()`, `print()` — coefficient table and model summary.
- `coef()`, `vcov()`, `confint()` — estimates, covariance matrices,
  and confidence intervals.
- `fitted()`, `residuals()`, `predict()` — fitted values and
  predictions.
- `model.matrix()` — design matrix.
- `tidy()`, `glance()` — tidy summaries following
  [broom](https://broom.tidymodels.org/) conventions.

The `cov_type` argument controls the covariance estimator used for
inference: `"robust"` (sandwich, default), `"bias-corrected"`,
`"df-adjusted"`, or `"naive"` (model-based).

### Model building and selection

- `anova()` — sequential or multi-model hypothesis test tables.
- `add1()`, `drop1()` — single-term additions and deletions with
  hypothesis tests and CIC.
- `step_p()` — stepwise model selection by hypothesis testing.
- `geecriteria()` — QIC, CIC, RJC, QICu, GESSC, and GPC model
  selection criteria.

### emmeans support

Fitted `geer` objects are compatible with the
[emmeans](https://cran.r-project.org/package=emmeans) package for
estimated marginal means.

## Datasets

The package includes seven example datasets: `cerebrovascular`,
`cholecystectomy`, `depression`, `epilepsy`, `leprosy`, `respiratory`,
and `rinse`.

## References

Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
generalized linear models. *Biometrika*, **73**, 13--22.

Touloumis, A. (2026) Bias reduction in generalized estimating
equations. Preprint.

Touloumis, A. (2026) Jeffreys-prior penalized GEE for correlated
binary data with an odds-ratio parameterization. Preprint.
