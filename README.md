<!-- badges: start -->
[![R-CMD-check](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml)
[![R release](https://img.shields.io/badge/R%20release-4.6.1-276DC3?logo=R&logoColor=white)](https://www.r-project.org/)
[![Codecov test coverage](https://codecov.io/gh/AnestisTouloumis/geer/graph/badge.svg)](https://app.codecov.io/gh/AnestisTouloumis/geer)
<!-- badges: end -->

## Overview

`geer` fits marginal models for independent, repeated, or clustered
responses using Generalized Estimating Equations (GEE). Supported
estimation methods include the traditional GEE, bias-reducing GEE,
bias-corrected GEE, and Jeffreys-prior penalized GEE. Continuous,
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
| `"pgee-jeffreys"` | Fully iterated Jeffreys-prior penalized GEE |
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
- `little_mcar_test()` — Little's test for assessing whether repeated-response
  missingness is compatible with MCAR.
- `mcar_homoscedasticity_test()` — Jamshidian-Jalal screening diagnostic for
  MCAR based on covariance homogeneity across missingness-pattern groups, with
  modified Hawkins and nonparametric Anderson-Darling options.
- `mcar_logistic_test()` — Ridout-style longitudinal MCAR diagnostic that
  models response missingness from the previous observed response and covariates
  using `geewa_binary()`.
- `model.matrix()` — design matrix.
- `tidy()`, `glance()` — tidy summaries following
  [broom](https://broom.tidymodels.org/) conventions.

The `cov_type` argument controls the covariance estimator used for
inference: `"bias-corrected"` (default), `"robust"` (sandwich),
`"df-adjusted"`, or `"naive"` (model-based).

The Wald-Wolfowitz residual runs test can be used as a quantitative
check for non-random residual sign patterns:

``` r
runs_test(fit)
runs_test(fit, type = "deviance", order_by = "fitted")
```

Cluster-level Mahalanobis residuals can be used to identify clusters whose
response profiles fit the working mean/covariance model poorly:

``` r
residuals(fit, type = "mahalanobis")
```

For longitudinal data with missing responses, Little's MCAR test can be
applied directly to a fitted model. The test reconstructs the wide response
matrix from the original `id` and `repeated` variables:

``` r
little_mcar_test(fit)
```

The same function can also be used with a numeric wide-format matrix or data
frame whose rows are independent units and columns are repeated measurements.
The implementation uses Little's degrees-of-freedom-corrected covariance
matrix. For the special bivariate monotone pattern, the exact normal-theory
F reference from Little (1988) is used automatically; otherwise the usual
large-sample chi-squared reference is used. Use `reference = "asymptotic"` to
force the chi-squared reference. The procedure is intended primarily for
quantitative variables and warns when binary variables are detected.

A distribution-robust screening diagnostic based on Jamshidian and Jalal
(2010) is also available:

``` r
mcar_homoscedasticity_test(fit)
mcar_homoscedasticity_test(fit, method = "nonparametric")
mcar_homoscedasticity_test(fit, method = "hawkins", imputation = "normal")
```

Cases are grouped by their original missingness pattern and the incomplete
responses are imputed before covariance homogeneity is assessed. The default
`method = "auto"` reports the modified Hawkins normality/homoscedasticity test
and uses the nonparametric k-sample Anderson-Darling component to distinguish
nonnormality from covariance heterogeneity when Hawkins rejects.
Distribution-free residual-resampling imputation is the default when there are
at least 10 complete cases and at least `2 * p` complete cases; otherwise the
function warns and falls back to conditional normal imputation. Patterns with
fewer than seven cases are omitted by default. The diagnostic assumes that,
apart from missingness patterns, cases arise from a common population, so known
groups with genuinely different covariance matrices can trigger rejection even
under MCAR. The reported test uses one imputed data set; multiple imputation is
best treated as a sensitivity diagnostic for imputation variability. This test
ignores the fitted regression structure and is therefore a screening device
that complements, rather than replaces, the regression-based MCAR diagnostic
below.

A complementary Ridout-style diagnostic models longitudinal missingness
on a transition risk set. At occasion `t`, a row is included when the response
at `t - 1` is observed; the binary outcome records whether the response at `t`
is missing. The previous response is included automatically, while occasion
effects are treated as nuisance terms:

``` r
out <- mcar_logistic_test(fit)
out$tests
mcar_logistic_test(fit, formula = ~ treatment + age, test = "score")
mcar_logistic_test(fit, test = "working-score", pmethod = "satterthwaite")
```

The missingness model is fitted with `geewa_binary(..., link = "logit",
method = "gee")`. The default odds-ratio structure is independence and the
default covariance is the bias-corrected estimator. The `test` argument
accepts the same five procedures as `anova.geer()`: `"wald"`, `"score"`,
`"working-wald"`, `"working-score"`, and `"working-lrt"`. The modified
working tests support the Rao-Scott and Satterthwaite p-value approximations
through `pmethod`; the modified working LRT requires independence. The printed
`htest` result is the response-history test, and `out$tests` also reports the
covariate and overall tests using the same selected procedure. A significant response-history effect is
evidence against covariate-dependent MCAR/random dropout. Covariate dependence
without response-history dependence is compatible with covariate-dependent
MCAR, although not strict MCAR relative to those covariates. Main occasion
effects are nuisance parameters and covariate columns aliased with them are
omitted automatically. Failure to reject does not establish MCAR or rule out
MNAR.

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
for `geer` models. Use `vcov = "robust"` for the sandwich covariance
matrix. Other `geer` covariance estimators can be supplied as a function,
for example `vcov = function(x) vcov(x, cov_type = "naive")`.

## Datasets

The package includes seven example datasets: `cerebrovascular`,
`cholecystectomy`, `depression`, `epilepsy`, `leprosy`, `respiratory`,
and `rinse`.

## References

Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
generalized linear models. *Biometrika*, **73**, 13--22.

Rubin, D.B. (1976) Inference and missing data. *Biometrika*, **63**,
581–592.

Little, R.J.A. (1988) A test of missing completely at random for multivariate
data with missing values. *Journal of the American Statistical Association*,
**83**, 1198–1202.

Jamshidian, M. and Jalal, S. (2010) Tests of homoscedasticity, normality, and
missing completely at random for incomplete multivariate data. *Psychometrika*,
**75**, 649–674.

Jamshidian, M., Jalal, S. and Jansen, C. (2014) MissMech: An R package for
testing homoscedasticity, multivariate normality, and missing completely at
random (MCAR). *Journal of Statistical Software*, **56**(6), 1–31.

Ridout, M.S. (1991) Testing for random dropouts in repeated measurement data.
*Biometrics*, **47**, 1617–1619.

Fitzmaurice, G.M., Heath, A.F. and Clifford, P. (1996) Logistic regression
models for binary panel data with attrition. *Journal of the Royal Statistical
Society: Series A*, **159**, 249–263.

Carpenter, J.R. and Smuk, M. (2021) Missing data: A statistical framework for
practice. *Biometrical Journal*, **63**, 915–947.

Lu, P. and Shelley, M. (2023) Testing the missingness mechanism in longitudinal
surveys: a case study using the Health and Retirement Study. *International
Journal of Social Research Methodology*, **26**, 439–452.

Chang, Y.-C. (2000) Residuals analysis of the generalized linear models
for longitudinal data. *Statistics in Medicine*, **19**, 1277--1293.

Hardin, J.W. and Hilbe, J.M. (2013) *Generalized Estimating Equations*,
2nd ed. Chapman and Hall/CRC.

Touloumis, A. (2026) [Bias-Reduced GEE via Adjusted Estimating Equations, with Odds-Ratio Extensions.](https://arxiv.org/abs/2606.16043) *Preprint*.

Touloumis, A. (2026) [Jeffreys-Type Penalized GEE for Correlated Binary Data with an Odds-Ratio Parameterization.](https://arxiv.org/abs/2606.16058) *Preprint*.
