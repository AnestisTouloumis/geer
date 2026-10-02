<!-- badges: start -->
[![R-CMD-check](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml)
[![R release](https://img.shields.io/badge/R%20release-4.6.1-276DC3?logo=R&logoColor=white)](https://www.r-project.org/)
[![Codecov test coverage](https://codecov.io/gh/AnestisTouloumis/geer/graph/badge.svg)](https://app.codecov.io/gh/AnestisTouloumis/geer)
<!-- badges: end -->

## Overview

`geer` fits marginal models for independent, repeated, or clustered
responses using Generalized Estimating Equations (GEE). Supported
estimation methods include the traditional GEE, bias-reducing GEE,
bias-corrected GEE, and Jeffreys-type penalized GEE. Continuous,
binary, and count responses are handled by `geewa`, while binary
responses can also be handled by `geewa_binary` through an odds-ratio
parameterization.

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

### Quick example

Fit a bias-reducing GEE with an exchangeable working correlation to the
epilepsy seizure count data:

``` r
data("epilepsy", package = "geer")

fit <- geewa(
  formula = seizures ~ treatment + lnbaseline + lnage,
  family = poisson(link = "log"),
  data = epilepsy,
  id = id,
  corstr = "exchangeable",
  method = "brgee-robust"
)
summary(fit, cov_type = "bias-corrected")
```

For binary responses, use `geewa_binary()` with an odds-ratio
parameterization:

``` r
data("cerebrovascular", package = "geer")

fit_bin <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  link = "logit",
  data = cerebrovascular,
  id = id,
  orstr = "exchangeable",
  method = "brgee-robust"
)
summary(fit_bin, cov_type = "bias-corrected")
```

### Fitting models

There are two core fitting functions:

- `geewa()` for continuous, binary, and count responses (Gaussian, Poisson,
  binomial, Gamma, inverse Gaussian, quasi, quasibinomial, and
  quasipoisson families).
- `geewa_binary()` for binary responses via a marginalized odds-ratio
  parameterization.

Both functions support the following estimation methods via the
`method` argument:

| Method | Description |
|---|---|
| `"gee"` | Traditional GEE |
| `"brgee-robust"`, `"brgee-naive"`, `"brgee-empirical"` | Bias-reducing GEE (differing in the bias adjustment used: robust, model-based, or empirical) |
| `"bcgee-robust"`, `"bcgee-naive"`, `"bcgee-empirical"` | Bias-corrected GEE (one-step correction; same three variants) |
| `"pgee-jeffreys"` | Fully iterated Jeffreys-type penalized GEE |
| `"opgee-jeffreys"` | One-step penalized GEE |
| `"hpgee-jeffreys"` | Hybrid one-step GEE |

The working correlation structure for `geewa()` is controlled by
`corstr`: `"independence"`, `"exchangeable"`, `"ar1"`,
`"m-dependent"`, `"unstructured"`, `"toeplitz"`, and `"fixed"`. The working
odds-ratio structure for `geewa_binary()` is controlled by `orstr`:
`"independence"`, `"exchangeable"`, `"unstructured"`, and `"fixed"`.

Convergence and fitting options are set via `geer_control()`.

### Inference

Standard S3 methods are available for fitted `geer` objects:

- `summary()`, `print()` — coefficient table and model summary.
- `coef()`, `vcov()`, `confint()` — estimates, covariance matrices,
  and confidence intervals.
- `fitted()`, `residuals()`, `predict()` — fitted values, observation-level
  residuals, cluster-level Mahalanobis residuals, and predictions.
- `runs_test()` — Wald-Wolfowitz test for non-random residual sign
  sequences, with natural, fitted-value, and covariate-based ordering.
- `mcar_little_test()` — Little's test for assessing whether repeated-response
  missingness is compatible with MCAR.
- `mcar_homoscedasticity_test()` — Jamshidian-Jalal screening diagnostic for
  MCAR based on covariance homogeneity across missingness-pattern groups, with
  modified Hawkins and nonparametric Anderson-Darling options.
- `mcar_logistic_test()` — Ridout-style longitudinal MCAR diagnostic that
  models response missingness from the previous observed response and covariates
  using `geewa_binary()`.
- `frechet_bounds_cor()` — Frechet bounds on the working correlations of a
  binomial `geewa()` fit, with a count of clusters violating them.
- `model.matrix()` — design matrix.
- `tidy()`, `glance()` — tidy summaries following
  [broom](https://broom.tidymodels.org/) conventions.

The `cov_type` argument controls the covariance estimator used for
inference: `"bias-corrected"` (default), `"robust"` (sandwich),
`"df-adjusted"`, `"jackknife"` (full leave-one-cluster jackknife), or
`"naive"` (model-based). The jackknife option computes each
leave-one-cluster estimate by refitting the model in full rather than by a
one-step approximation. The working association structure is unchanged and
the applicable association parameters are held fixed at their full-data
estimates in every deletion fit; it requires at least two clusters and
signals an error if any deletion refit fails to converge. The same
`"jackknife"` choice is accepted by every public function that exposes
`cov_type`, including model-comparison, selection-criterion, stepwise-selection,
and MCAR-regression diagnostics.

``` r
vcov(fit, cov_type = "jackknife")
summary(fit, cov_type = "jackknife")
```

The Wald-Wolfowitz residual runs test can be used as a quantitative
check for non-random residual sign patterns:

``` r
runs_test(fit)
runs_test(fit, order_by = "fitted")
```

Cluster-level Mahalanobis residuals can be used to identify clusters whose
response profiles fit the working mean/covariance model poorly:

``` r
residuals(fit, type = "mahalanobis")
```

For longitudinal data with missing responses, three diagnostics are available
for assessing whether the missingness is compatible with MCAR. They need a
fitted model, or wide-format data, that actually contains missing responses,
so none is run here:

- `mcar_little_test()` implements Little's (1988) test and accepts a fitted
  `geer` object or a numeric wide-format matrix or data frame. It uses the
  exact normal-theory F reference in the bivariate monotone case and the
  chi-squared reference otherwise; `reference = "asymptotic"` forces the
  latter.
- `mcar_homoscedasticity_test()` is the Jamshidian-Jalal (2010) screening
  diagnostic, following the MissMech implementation of Jamshidian, Jalal and
  Jansen (2014). It accepts the same inputs as `mcar_little_test()`. The
  `method` argument selects the modified Hawkins test, the nonparametric
  Anderson-Darling test, or the default `"auto"` logic; `imputation` selects
  residual-resampling or normal-theory imputation; `n_imputations` repeats the
  tests over several imputations; and `imputed_data` accepts a completed data
  set from another imputation method.
- `mcar_logistic_test()` is a Ridout-style (1991) regression diagnostic that
  models missingness at occasion `t` from the response at `t - 1` and,
  optionally, covariates, using `geewa_binary()`. Its `test` argument accepts
  the same five procedures as `anova.geer()` and `pmethod` selects the
  Rao-Scott or Satterthwaite approximation for the working tests.

None of these tests can establish MCAR; see the help pages for details.

### Model building and selection

- `anova()` — sequential or multi-model hypothesis test tables.
- `add1()`, `drop1()` — single-term additions and deletions with
  hypothesis tests and CIC.
- `step_p()` — stepwise model selection by hypothesis testing.
- `geecriteria()` — QIC, QICHH, QICC, CIC, RJC, QICu, EQIC, GESSC, GPC,
  AGPC, SGPC, GHYC, and PAC model selection criteria.
  For `geecriteria()`, `cov_type = "robust"` is the default so the classical
  covariance-based definitions are returned unless another covariance estimator
  is requested explicitly.

### Post-estimation support

Fitted `geer` objects are compatible with the
[emmeans](https://cran.r-project.org/package=emmeans) package for
estimated marginal means.

They are also compatible with the
[marginaleffects](https://marginaleffects.com/) package for predictions,
average predictions, comparisons, and slopes (marginal effects):

``` r
# install.packages("marginaleffects")
marginaleffects::avg_predictions(fit_bin)
marginaleffects::avg_comparisons(fit_bin, variables = "treatment")
marginaleffects::avg_slopes(fit, variables = "lnage")
```

By default, `marginaleffects` uses the bias-corrected covariance matrix
for `geer` models. Character values such as `vcov = "robust"`,
`vcov = "df-adjusted"`, `vcov = "jackknife"`, and `vcov = "naive"`
select the corresponding `geer` covariance estimator.

## Datasets

The package includes seven example datasets: `cerebrovascular`,
`cholecystectomy`, `depression`, `epilepsy`, `leprosy`, `respiratory`,
and `rinse`.

## References

Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
generalized linear models. *Biometrika*, **73**, 13--22.

Vanegas, L.H., Rondon, L.M. and Paula, G.A. (2023) Generalized Estimating
Equations using the new R package glmtoolbox. *The R Journal*, **15**, 105--133.

Rubin, D.B. (1976) Inference and missing data. *Biometrika*, **63**,
581--592.

Little, R.J.A. (1988) A test of missing completely at random for multivariate
data with missing values. *Journal of the American Statistical Association*,
**83**, 1198--1202.

Jamshidian, M. and Jalal, S. (2010) Tests of homoscedasticity, normality, and
missing completely at random for incomplete multivariate data. *Psychometrika*,
**75**, 649--674.

Jamshidian, M., Jalal, S. and Jansen, C. (2014) MissMech: An R package for
testing homoscedasticity, multivariate normality, and missing completely at
random (MCAR). *Journal of Statistical Software*, **56**, 1--31.

Ridout, M.S. (1991) Testing for random dropouts in repeated measurement data.
*Biometrics*, **47**, 1617--1619.

Fitzmaurice, G.M., Heath, A.F. and Clifford, P. (1996) Logistic regression
models for binary panel data with attrition. *Journal of the Royal Statistical
Society: Series A*, **159**, 249--263.

Carpenter, J.R. and Smuk, M. (2021) Missing data: A statistical framework for
practice. *Biometrical Journal*, **63**, 915--947.

Lu, P. and Shelley, M. (2023) Testing the missingness mechanism in longitudinal
surveys: a case study using the Health and Retirement Study. *International
Journal of Social Research Methodology*, **26**, 439--452.

Chang, Y.-C. (2000) Residuals analysis of the generalized linear models
for longitudinal data. *Statistics in Medicine*, **19**, 1277--1293.

Hardin, J.W. and Hilbe, J.M. (2013) *Generalized Estimating Equations*,
2nd Edition. Chapman and Hall/CRC, Boca Raton.

Touloumis, A. (2026) [Bias-reduced GEE via adjusted estimating equations, with odds-ratio extensions.](https://arxiv.org/abs/2606.16043) *Preprint*.

Touloumis, A. (2026) [Jeffreys-type penalized GEE for correlated binary data with an odds-ratio parameterization.](https://arxiv.org/abs/2606.16058) *Preprint*.
