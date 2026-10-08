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

* The dispersion estimate is no longer rejected when it is merely small
  (previously below `.Machine$double.eps`); only zero or non-finite estimates
  stop the fit, so responses measured in very small units can be fitted.
* Unstructured and fixed working correlation structures no longer fail with an
  obscure internal error when every cluster has a single observation.
* The solvers treat a non-finite Newton step as a numerical failure (reverting
  to the last accepted iterate with a warning) instead of accepting it.
* `get_vcov()` (marginaleffects) now calls a function supplied as `vcov` on the
  fitted model, as documented, instead of ignoring it; the result is validated
  like a matrix.
* The `Step` column of the `step_p()` table is now an integer (`NA` for the
  initial model) instead of character, so `print()` no longer shows factor codes
  or sorts steps lexicographically.
* `emmeans` support: the prior weights passed to `recover_data()` are now put
  back in the original row order of the data (fits store their rows sorted by
  cluster and occasion), which matters for unsorted data with weights or
  binomial trials. Fits gain an optional `row_order` component for this.

* `glance()` computes only QIC, QICu and CIC, so that a failure in an unrelated
  criterion no longer turns them into `NA`.

* `step_p(direction = "forward")` without `scope` now warns that there are no
  candidate terms to add.

* The refits used by `add1()`, `drop1()`, `anova()` and `step_p()` no longer
  reuse a `beta_start` of the original fit, whose length does not match a
  model with added or dropped terms.

* `residuals.geer()` gains `type = "response"`, an alias of `"working"` (the raw
  residuals), and the documentation of the score tests now states that
  `alpha` and `phi` are those of the larger model, not re-estimated under the
  null model.

* Removed the build-time fields `Author`, `Maintainer` and `Packaged` from
  `DESCRIPTION`, removed `skip_if_not_installed("brglm2")` calls that could
  never skip (`brglm2` is in Imports), and added the minimum version to the
  `emmeans` skip in the tests.

* `geewa()` now rejects a factor, character or multi-column matrix response
  for families that are not binomial-type, with a clear message. Before, a
  factor was silently fitted as its level codes and a matrix produced a
  misleading "response variable and 'id' are not of same length" error.

* For `method = "opgee-jeffreys"` and `"hpgee-jeffreys"` in `geewa()`, the
  working correlation parameters and the dispersion are now estimated once, at
  the converged independence penalized estimate, and held fixed in the one-step
  update, as for the bias-corrected methods. Before, the one-step regression
  estimate used these values but the returned `alpha`, `phi` and covariance
  matrices were re-estimated at the one-step estimate. The numerical results
  of these two methods therefore change slightly. The same applies to the
  leave-one-cluster refits of the jackknife covariance. The documentation of
  `converged` now states that it is `FALSE` after a numerical failure in the
  correction step, and the documentation describes how `alpha` and `phi` are
  obtained for each one-step method.

* A scalar `offset` argument in `geewa()` and `geewa_binary()`, which was
  documented but failed with a `model.frame()` length error, is now recycled to
  the number of rows of the data.

* `geewa()` and `geewa_binary()` now reject a `beta_start` containing `NA` or
  non-finite values, name the linearly dependent columns in the
  rank-deficiency error, and give a clear error for a model matrix without
  columns (for example `y ~ 0`).

* `geewa()` and `geewa_binary()` with an `"unstructured"` or `"fixed"`
  association structure no longer fail after fitting when every cluster has a
  single observation, and a missing convergence criterion is treated as
  non-convergence rather than producing an obscure `if` error.

* `mcar_logistic_test()` no longer fails on data with two occasions: when the
  null model has no terms, an intercept-only formula is used instead of
  `stats::reformulate(character(0))`.

* The refits used by `anova()`, `add1()`, `drop1()` and `step_p()` now stop,
  naming the refit formula, when the solver did not converge. Previously only
  `geewa()` warned and the test table was built from the last accepted
  iterate.

* The Wald, score and working tests used by `anova()`, `add1()`, `drop1()` and
  `step_p()` no longer abort when the test statistic is negative, which can
  happen with covariance estimates that are not positive semi-definite (the
  default bias-corrected one in small samples). A warning is issued and the
  statistic and p-value of that row are `NA`; `step_p()` already skips `NA`
  p-values. `mcar_logistic_test()` still stops in this case, because it needs a
  statistic.

* `geewa()` and `geewa_binary()` gain the `subset` and `na.action` arguments,
  which were previously listed internally but never reachable: with the default
  `control` they produced an "unused argument" error, and with explicit
  `control` and `control_glm` they were silently ignored. Arguments left in
  `...` are now rejected when both `control` and `control_glm` are supplied.

* The modified working likelihood-ratio test (`test = "working-lrt"`) no longer
  requires a unit dispersion for Poisson and binomial models. Each model's
  working log-likelihood is scaled by that model's own dispersion, whether
  estimated from the data or fixed, for every family and link.

* Fixed the marginalized odds-ratio solver's `dV/dmu` matrix for clusters of
  size one: the `-mu` terms and the unit diagonal were skipped for singleton
  clusters, giving a wrong derivative of the variance function.

* `geewa()` with `corstr = "fixed"` and a one-step estimator (`opgee-*`,
  `hpgee-*`) now uses the supplied association parameters in the second pass
  instead of resetting them to zero.

* `anova.geer()` now keeps formula `offset()` terms when refitting the null
  model, which affected score tests (and Wald-free comparisons) for models
  with an offset in the formula.

* The design matrix of a one-column model keeps its `assign` attribute and
  column name; the former one-column reshaping was removed.

* The bias-corrected covariance estimator (Morel et al., 2003) is now set to
  `NA` when the number of clusters does not exceed the number of regression
  parameters, and `geewa()` and `geewa_binary()` warn once. Previously, with
  fewer clusters than parameters the shrinkage term was negative and the
  "variance" could be negative, and with exactly as many clusters as parameters
  the shrinkage was silently capped at 0.5. The robust, naive, and
  `"df-adjusted"` covariances are unaffected (`"df-adjusted"` already required
  more clusters than parameters), so use `cov_type = "robust"` or
  `"jackknife"` in such cases.

* Numerical failures inside the C++ estimation routines now report where they
  occurred. An error raised while processing a cluster names the routine and the
  cluster (its position and id value), and failures of the final information
  matrix solves name the estimator. Checks of the linear predictor and fitted
  values no longer copy the vectors at every trial step.

* The `step_multiplier` control of `geer_control()` is now a strictly positive
  real number instead of a positive integer. Values below 1 damp the proposed
  step, which can help when the algorithm diverges from poor starting values
  (it plays the role of `slowit` in `brglm2::brglm_control()`); values above 1
  enlarge it, as before. The default is unchanged.

* A numerical failure while the iterative algorithm evaluates a trial point
  (a singular matrix, an association parameter outside its admissible range,
  or step-halving attempts that all leave the valid region) no longer aborts
  the fit with an error. The algorithm stops, reverts to the last accepted
  iterate, and `geewa()` and `geewa_binary()` warn with the reason and report
  `converged = FALSE`; for the one-step estimators (`"bcgee-*"`,
  `"hpgee-jeffreys"`, `"opgee-jeffreys"`) the returned estimates are then the
  preliminary first-stage estimates. Failures at the starting values still
  signal an error, and jackknife refits and the first stage of the one-step
  estimators still stop with an error that now includes the reason.
  Deviation from `brglm2::brglmFit()`, which reverts to the previous inner
  trial: the last accepted iterate is used, because rejected trials can be
  worse than the accepted point.

* The C++ odds-ratio solver now rejects a response that is not finite or lies
  outside `[0, 1]` with an informative error instead of fitting it. The R
  interface already enforced this domain, so `geewa_binary()` is unaffected;
  the check guards direct calls to the internal routine.

* The C++ solvers are faster and use less memory. The Newton step computed at
  the accepted candidate is reused as the starting step of the next iteration
  instead of being recomputed, the Kronecker-product contractions in the
  bias-reduced and penalized updates no longer form `m^2 x p^2` matrices, each
  working covariance matrix is factorized once per cluster, and the iteration
  history is stored in a growable container rather than a
  `p x (maxiter + 1)` matrix. Estimates, covariance matrices, and the iteration
  path are unchanged up to floating-point rounding.

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
