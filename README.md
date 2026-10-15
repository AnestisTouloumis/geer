
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# geer: Bias-Reduced and Penalized Generalized Estimating Equations

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/geer)](https://CRAN.R-project.org/package=geer)
[![R-CMD-check](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AnestisTouloumis/geer/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/AnestisTouloumis/geer/graph/badge.svg)](https://app.codecov.io/gh/AnestisTouloumis/geer)
<!-- badges: end -->

`geer` fits marginal models for clustered or repeated responses by
generalized estimating equations (GEE; Liang and Zeger, 1986). Besides
ordinary GEE, it provides

- bias-reduced and bias-corrected GEE, which reduce the finite-sample
  bias of the regression estimates when the number of clusters is small
  (Touloumis, 2026a), and
- Jeffreys-type penalized GEE, which, for binary responses fitted by
  `geewa_binary()`, also guard against infinite estimates under
  separation (Touloumis, 2026b).

`geewa()` fits continuous, binary, and count responses with a working
correlation structure, and `geewa_binary()` fits binary responses with a
working odds-ratio structure. By default, standard errors use the
bias-corrected covariance matrix of Morel, Bokossa and Neerchal (2003).

## Installation

Install the released version from CRAN:

``` r
install.packages("geer")
```

or the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("AnestisTouloumis/geer")
```

## Example

Ordinary and bias-reduced GEE for the epilepsy seizure counts:

``` r
library("geer")

data("epilepsy", package = "geer")
fit <- geewa(
  formula = seizures ~ treatment + lnbaseline + lnage,
  family = poisson(link = "log"),
  data = epilepsy,
  id = id,
  repeated = visit,
  corstr = "exchangeable",
  method = "gee"
)
fit_br <- update(fit, method = "brgee-robust")
summary(fit_br)
#> 
#> Call:
#> geewa(formula = seizures ~ treatment + lnbaseline + lnage, family = poisson(link = "log"), 
#>     data = epilepsy, id = id, repeated = visit, corstr = "exchangeable", 
#>     method = "brgee-robust")
#> 
#> Estimating Method   : brgee-robust 
#> Number of iterations: 4 
#> Algorithm converged : TRUE 
#> 
#> Number of observations: 236 
#> Number of clusters    : 59 
#> Cluster sizes         : 4 to 4 
#> 
#> Marginal Model
#> Family       : poisson 
#> Link Function: log 
#> 
#> Coefficients:
#>                     Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)        -2.370155   1.133494 -2.0910  0.03653 *  
#> treatmentprogabide -0.013753   0.195643 -0.0703  0.94396    
#> lnbaseline          1.229923   0.153520  8.0115  <0.0001 ***
#> lnage               0.577250   0.312313  1.8483  0.06456 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Std. Errors are taken from the bias-corrected covariance matrix. 
#> 
#> Dispersion Parameter: 4.7466 
#> 
#> Association Structure: exchangeable 
#> Association Parameter: 0.3903
```

Penalized GEE for a binary response with an exchangeable odds-ratio
structure:

``` r
data("cerebrovascular", package = "geer")
fit_bin <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  link = "logit",
  data = cerebrovascular,
  id = id,
  repeated = period,
  orstr = "exchangeable",
  method = "pgee-jeffreys"
)
summary(fit_bin)
#> 
#> Call:
#> geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit", 
#>     data = cerebrovascular, id = id, repeated = period, orstr = "exchangeable", 
#>     method = "pgee-jeffreys")
#> 
#> Estimating Method   : pgee-jeffreys 
#> Number of iterations: 4 
#> Algorithm converged : TRUE 
#> 
#> Number of observations: 134 
#> Number of clusters    : 67 
#> Cluster sizes         : 2 to 2 
#> 
#> Marginal Model
#> Family       : binomial 
#> Link Function: logit 
#> 
#> Coefficients:
#>                 Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)      1.18391    0.30443  3.8890 0.0001007 ***
#> treatmentactive -0.54334    0.23944 -2.2692 0.0232566 *  
#> factor(period)2 -0.28103    0.23828 -1.1794 0.2382282    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Std. Errors are taken from the bias-corrected covariance matrix. 
#> 
#> Dispersion Parameter: 1 
#> 
#> Association Structure: exchangeable 
#> Association Parameter: 16.4118
```

## Main features

**Estimation methods** (argument `method` of `geewa()` and
`geewa_binary()`):

| Method | Description |
|----|----|
| `"gee"` | Ordinary GEE |
| `"brgee-naive"`, `"brgee-robust"`, `"brgee-empirical"` | Bias-reduced GEE |
| `"bcgee-naive"`, `"bcgee-robust"`, `"bcgee-empirical"` | Bias-corrected GEE |
| `"pgee-jeffreys"` | Jeffreys-type penalized GEE |
| `"opgee-jeffreys"`, `"hpgee-jeffreys"` | One-step and hybrid one-step penalized GEE |

**Working association structures:** `corstr` in `geewa()`
(`"independence"`, `"exchangeable"`, `"ar1"`, `"m-dependent"`,
`"unstructured"`, `"toeplitz"`, or `"fixed"`) and `orstr` in
`geewa_binary()` (`"independence"`, `"exchangeable"`, `"unstructured"`,
or `"fixed"`).

**Inference:** `summary()`, `vcov()`, `confint()`, `predict()`,
`tidy()`, and `glance()`, with the covariance matrix estimator selected
by `cov_type`: `"bias-corrected"` (the default), `"robust"`,
`"df-adjusted"`, `"jackknife"`, or `"naive"`.

**Model building:** `anova()`, `add1()`, `drop1()`, and `step_p()`, with
Wald, generalized score, and modified working Wald, score, and
likelihood ratio tests (Rotnitzky and Jewell, 1990), and
`geecriteria()`, with 16 model selection criteria, including QIC (Pan,
2001) and CIC (Hin and Wang, 2009).

**Diagnostics:**

- `residuals()`, including cluster-level Mahalanobis residuals.
- `runs_test()`, the Wald-Wolfowitz runs test of residual signs (Chang,
  2000, and Hardin and Hilbe, 2013).
- `frechet_bounds_cor()`, which checks the working correlations of
  binary `geewa()` fits.
- `mcar_little_test()` (Little, 1988), `mcar_logistic_test()` (in the
  spirit of Ridout, 1991), and `mcar_homoscedasticity_test()`
  (Jamshidian and Jalal, 2010, and Jamshidian, Jalal and Jansen, 2014),
  three tests of missing completely at random (MCAR).

**Integration with other packages:** fitted models work with
[emmeans](https://CRAN.R-project.org/package=emmeans) and
[marginaleffects](https://CRAN.R-project.org/package=marginaleffects).

**Data:** `cerebrovascular`, `cholecystectomy`, `depression`,
`epilepsy`, `leprosy`, `respiratory`, and `rinse`.

## Citation

To cite `geer` and the methods it implements, use

``` r
citation("geer")
```

## Getting help

Please report bugs and request features at
<https://github.com/AnestisTouloumis/geer/issues>.

## References

Chang, Y.-C. (2000) Residuals analysis of the generalized linear models
for longitudinal data. *Statistics in Medicine*, **19**, 1277–1293.

Hardin, J.W. and Hilbe, J.M. (2013) *Generalized Estimating Equations*,
2nd Edition. Chapman and Hall/CRC, Boca Raton.

Hin, L.Y. and Wang, Y.G. (2009) Working-correlation-structure
identification in generalized estimating equations. *Statistics in
Medicine*, **28**, 642–658.

Jamshidian, M. and Jalal, S. (2010) Tests of homoscedasticity,
normality, and missing completely at random for incomplete multivariate
data. *Psychometrika*, **75**, 649–674.

Jamshidian, M., Jalal, S. and Jansen, C. (2014) MissMech: An R package
for testing homoscedasticity, multivariate normality, and missing
completely at random (MCAR). *Journal of Statistical Software*, **56**,
1–31.

Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
generalized linear models. *Biometrika*, **73**, 13–22.

Little, R.J.A. (1988) A test of missing completely at random for
multivariate data with missing values. *Journal of the American
Statistical Association*, **83**, 1198–1202.

Morel, J.G., Bokossa, M.C. and Neerchal, N.K. (2003) Small sample
correction for the variance of GEE estimators. *Biometrical Journal*,
**45**, 395–409.

Pan, W. (2001) Akaike’s information criterion in generalized estimating
equations. *Biometrics*, **57**, 120–125.

Ridout, M.S. (1991) Testing for random dropouts in repeated measurement
data. *Biometrics*, **47**, 1617–1619.

Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
parameters in semiparametric generalized linear models for cluster
correlated data. *Biometrika*, **77**, 485–497.

Touloumis, A. (2026a) Bias-reduced GEE via adjusted estimating
equations, with odds-ratio extensions. *arXiv preprint*.
[doi:10.48550/arXiv.2606.16043](https://doi.org/10.48550/arXiv.2606.16043)

Touloumis, A. (2026b) Jeffreys-type penalized GEE for correlated binary
data with an odds-ratio parameterization. *arXiv preprint*.
[doi:10.48550/arXiv.2606.16058](https://doi.org/10.48550/arXiv.2606.16058)
