# geer 0.1.1

## New features

* Added `cov_type = "jackknife"` throughout the package wherever `cov_type` is
  accepted, including `vcov()`, summaries, confidence intervals, prediction,
  tidiers, `geecriteria()`, `add1()`, `drop1()`, `anova()`, `step_p()`, and
  `mcar_logistic_test()`. The estimator refits the regression parameters after
  deleting each cluster in turn while holding the working association structure
  and association-parameter vector fixed at their full-data values; it uses the
  centered Quenouille-Tukey jackknife covariance, including the `(K - 1) / K`
  finite-sample factor, where `K` is the number of clusters, and each
  leave-one-cluster estimate is obtained by a full refit rather than by a
  one-step approximation. It is available for all estimation methods.
  Score-based procedures retain their null-model score/information calculation
  and use the larger model's full-refit jackknife covariance as the covariance
  component.

* Added optional integration with the `marginaleffects` package. `geer` model
  objects can now be used with `marginaleffects` workflows for predictions,
  comparisons, slopes, and their averaged counterparts. The integration uses
  delayed S3 registration, so `marginaleffects` remains an optional dependency,
  and includes model-data recovery through `insight`.

* Expanded `geecriteria()` with additional criteria for working association
  structure and model selection: `QICHH`, `QICC`, `EQIC`, `GHYC`, `PAC`, `AGPC`,
  and `SGPC`, which complement the existing `QIC`, `CIC`, `RJC`, `QICu`,
  `GESSC`, and `GPC` criteria, and the three eigenvalue-based criteria of Jang
  (2011): `PT` (Pillai trace type), `WR` (Wilks ratio type), and `RMR` (Roy
  maximum root type). All three are functions of the generalized eigenvalues of
  the covariance estimate with respect to the model-based covariance under
  working independence, respond to `cov_type`, and prefer smaller values. The
  new `criteria` argument selects which criteria to report. The default,
  `"all"`, returns every criterion; otherwise a character vector of criterion
  names is accepted, matched ignoring case, with the columns returned in the
  order requested and `Parameters` always included. Only the criteria named in
  `criteria` are evaluated, so restricting the selection also avoids the work
  the others would require: the working-independence refit behind `QICHH`, the
  per-cluster loop behind `GHYC` and `PAC`, the eigendecomposition behind `PT`,
  `WR`, and `RMR`, and, when none of the covariance-based criteria is requested,
  the sandwich covariance itself. Criteria that cannot be evaluated for a
  particular fit, such as `RJC`, `QICHH`, and `EQIC` for a degenerate candidate
  model, are reported as `NA` instead of aborting the whole table; errors from
  argument validation are unaffected. The `EQIC` dispersion is not floored at
  machine precision, and an unusable value yields `NA`. `GHYC` and `PAC` are
  reported only for balanced designs, that is when every cluster observes each
  repeated position exactly once, and are `NA` otherwise, because both sum
  cluster-level covariance matrices that are conformable only in that case and
  Gosho, Hamada and Yoshimura (2011) and Pardo and Alonso (2019) define them
  under a common cluster size. Supplying the same model twice no longer fails on
  duplicate row names.

* Added `runs_test()` for the Wald-Wolfowitz nonparametric runs test of GEE
  residual signs described by Chang (2000) and Hardin and Hilbe (2013). The test
  can assess the natural cluster/repeated order, fitted-value order, or
  covariate-based orderings. It has no residual-type argument, because the
  working, Pearson, and deviance residuals share a sign sequence and therefore
  give an identical statistic. The `alternative` argument selects a two-sided
  test (the default) or a one-sided test against too few or too many runs.
  Results are returned in the standard `"htest"` slots, so the number of runs,
  its null expectation, and the two sign counts are all shown by the default
  print method. A warning is issued when either sign count is 15 or fewer,
  because the normal approximation may then be unreliable. A character
  `order_by` naming a variable that is not a model-matrix column is resolved
  against the data used to fit the model, so residuals can be ordered on a
  covariate omitted from the model. The returned object also inherits from
  `"geer_runs_test"` and carries the tested sign sequence, so that `plot()`
  draws the residual runs figure of Hardin and Hilbe (2013, Section 4.2.1),
  coloring the points by run and, under the natural ordering, marking the
  cluster boundaries.

* Added `mcar_little_test()` implementing Little's (1988) test for assessing
  whether repeated-response missingness is compatible with MCAR. The function
  works directly with fitted `geer` objects or numeric wide-format data,
  estimates the common multivariate-normal mean and covariance by an internal EM
  algorithm, and applies Little's `n / (n - 1)` covariance correction in the
  test statistic. The exact normal-theory F reference is used automatically for
  the bivariate monotone case; otherwise the large-sample chi-squared reference
  is used. A warning is issued when binary variables are detected because Little
  recommends the procedure primarily for quantitative variables. Following
  Little (1988) exactly, rows with no observed values belong to no missing-data
  pattern; they are removed with a warning rather than being counted in `n`,
  where they would inflate the degrees-of-freedom correction and the EM
  denominators. Each pattern contribution is evaluated through a Cholesky
  factorization, so the statistic is nonnegative by construction and is not
  clamped at zero; a covariance submatrix that is not positive definite is
  reported instead.

* Added `mcar_homoscedasticity_test()` implementing the Jamshidian-Jalal (2010)
  MCAR screening framework as implemented in `MissMech` (Jamshidian, Jalal and
  Jansen, 2014). Cases are grouped by their original missingness pattern,
  completed with either distribution-free residual resampling or conditional
  normal imputation, and assessed using the modified Hawkins
  normality/homoscedasticity test and/or the nonparametric k-sample
  Anderson-Darling test. The default `method = "auto"` follows the published
  diagnostic logic, while small pattern groups are omitted using the same
  seven-case default threshold as `MissMech`. The function works with fitted
  `geer` objects or numeric wide-format data and is documented as a screening
  diagnostic that does not condition on the fitted GEE regression structure.
  `n_imputations` repeats the tests on several completed data sets and returns
  per-imputation statistics, p-values, pattern-specific Neyman p-values and
  Anderson-Darling contributions for the exploratory assessment described by
  Jamshidian and Jalal (2010); the primary result is based on the first
  imputation. `imputed_data` accepts a completed data set from another
  imputation method, which also allows a test of covariance homogeneity across
  known groups in complete data. The Anderson-Darling p-value uses the simulated
  reference quantiles and interpolation of the `kSamples` package rather than
  the five-point table of Scholz and Stephens (1987), so p-values can differ
  from `MissMech`, particularly below 0.01. Simulated Neyman p-values use the
  `(1 + b) / (nrep + 1)` form and the simulated null distribution is reused
  across imputations. A warning is issued when normal-theory imputation is
  requested for the nonparametric test, whose size it can inflate for nonnormal
  data.

* Added cluster-level Mahalanobis residuals through `residuals(..., type =
  "mahalanobis")`, using the fitted working covariance matrix.

## Changes

* Revised `mcar_logistic_test()` to use a Ridout-style longitudinal risk-set
  formulation. At occasion `t`, the missingness indicator is modeled only when
  the response at `t - 1` is observed; the immediately previous response is
  included automatically and occasion effects are nuisance parameters. The
  binary model is now fitted with `geewa_binary()` using a logit link and a
  working odds-ratio structure. The primary test assesses dependence on the
  previous response, with additional tests for observed-covariate effects and
  the joint effect of covariates and response history. All five hypothesis-test
  procedures implemented for `geer` models are available: Wald, generalized
  score, modified working Wald, modified working score, and modified working
  LRT. The Rao-Scott and Satterthwaite approximations are available for the
  modified working tests. The default covariance estimator for this diagnostic
  is now the bias-corrected covariance matrix, consistent with the package's
  default inferential emphasis. This distinguishes evidence against
  covariate-dependent MCAR from covariate-only departures from strict MCAR.
  Intermittent missingness is supported as a local transition diagnostic with a
  warning that the interpretation is no longer a pure dropout test. References
  now include Ridout (1991) and Fitzmaurice, Heath and Clifford (1996).

* Deviance residuals are now scaled by the fitted dispersion parameter,
  consistent with the GEE residual definition used by `glmtoolbox`.

* Changed the default covariance estimator in `geecriteria()` to `"robust"` so
  the classical forms of covariance-based GEE criteria are returned by default.
  Other covariance estimators remain available through `cov_type`.

* All complexity penalties in `geecriteria()` now count only the
  working-association parameters that were estimated from the data. `GESSC`
  previously divided its weighted error sum of squares by `N - p - q` with `q`
  the full length of the association parameter vector, while `QICC`, `AGPC`, and
  `SGPC` already used the estimated count. The two disagreed for `corstr =
  "fixed"` and `orstr = "fixed"`, where a supplied structure costs no degrees of
  freedom; `GESSC` values change for those fits only.

* P-values from the Wald, score, working Wald, working score and working LRT
  tests (used by `anova()`, `add1()`, `drop1()`, `step_p()` and
  `mcar_logistic_test()`) are now computed from the upper tail of the
  chi-square distribution. Previously they were computed as `1 - pchisq()`,
  which returns exactly 0 for test statistics larger than about 38 instead of
  the correct, very small p-value.

* Fixed how the fitting function (`geewa()` or `geewa_binary()`) and a fixed
  dispersion are recognized in score tests, `anova()`, `add1()`, `drop1()`,
  `step_p()`, `geecriteria()`, `glance()`, and Mahalanobis residuals. They
  previously parsed the stored call, so a model fitted as `geer::geewa(...)`, or
  with `phi_fixed` supplied as a variable, was handled as an odds-ratio fit or
  as having an estimated dispersion. The recorded `fit_function` and `phi_fixed`
  are now used, and an undeterminable fitting function is an error.

* `print()` of `summary()` and of the `anova`-type tables returned by `anova()`,
  `add1()`, `drop1()` and `step_p()` now shows p-values below 0.0001 as
  `<0.0001`. `stats`, `utils` and `grDevices` are declared in `Imports`.

* Internal refactor: the two-pass logic of `bcgee-*`, `hpgee-jeffreys` and
  `opgee-jeffreys` is now implemented once (`run_geer_estimation_passes()`)
  and shared by `geewa()`, `geewa_binary()` and the jackknife refits, so the
  jackknife can no longer drift from the main fit. Results are unchanged.

* Nested-model comparisons (`anova()`, score and Wald tests) now stop when the two
  fits differ in fitting function, estimation method, working association
  structure, dispersion handling, `use_p`, m-dependence order, fixed
  association parameters or contrasts, instead of comparing them silently.

* Standard errors from a negative or non-finite variance (the bias-corrected
  covariance is not guaranteed to be positive semi-definite) are now `NA` with a
  warning naming the affected coefficients in `summary()`, `tidy()` and
  `predict(se.fit = TRUE)`, instead of silent `NA`/`NaN`. `confint()` still stops.
* The leave-one-cluster jackknife covariance is cached on the fit, so repeated
  calls to `summary()`, `confint()`, `tidy()`, `predict()` and the tests no longer
  refit the model once per cluster each time. The cache is ignored when the
  coefficients change.
* `frechet_bounds_cor()` uses the same fit-function detection as the rest of the
  package; "Fréchet" is now typeset with `\enc{}` in the help page.
* Internal: input preparation shared by `geewa()` and `geewa_binary()` now lives
  in `prepare_geer_inputs()`.

* `print()` of `runs_test()`, `mcar_little_test()`, `mcar_logistic_test()` and
  `mcar_homoscedasticity_test()` results now shows p-values below 0.0001 as
  `p-value < 0.0001`. These objects gain the class `"geer_htest"` ahead of
  `"htest"`; stored p-values are unchanged.

* `step_p(direction = "both")` now requires `p_enter <= p_remove` and stops when a
  move would return to a model already visited, instead of cycling until `steps`
  is exhausted.

* Fixed the sign of the inverse-Gaussian quasi-log-likelihood. It entered `QIC`,
  `QICu`, `QICC`, `QICHH` and the modified working LRT with the wrong sign for
  `inverse.gaussian` fits; the other families were correct. The quasi-
  log-likelihood of every family is now tested against minus half the unit
  deviance.
* `geecriteria()` now also warns when the supplied models differ in distribution
  family or in response values, not only in the number of observations.
