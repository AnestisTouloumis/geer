# geer 0.1.1

## New features

* Added `cov_type = "jackknife"`, the leave-one-cluster-out jackknife covariance
  (each estimate a full refit with the association parameters held at their
  full-data values), available for all estimation methods and wherever
  `cov_type` is accepted. It is computed once and cached on the fit.

* Added nine small-sample sandwich estimators to `cov_type`, those of
  Kauermann and Carroll (2001), Mancl and DeRouen (2001), Ford and Westgate
  (2017), Fay and Graubard (2001), Fan, Zhang and Zhang (2013), Pan (2001),
  Gosho, Sato and Takeuchi (2014), Wang and Long (2011) and Westgate and
  Burchett (2016), for `geewa()` and `geewa_binary()` fits with any family and
  link. They are available wherever `cov_type` is accepted, including the score
  tests, where they are evaluated at the fit of the smaller model.

* Added `step_criterion()`, stepwise model selection that moves to the
  addition or deletion with the smallest QICu, QIC, QICC, AGPC or SGPC. The
  candidate table of every step is returned, and printed with `trace = TRUE`.

* Added `hatvalues()`, `dfbeta()`, `dfbetas()` and `cooks.distance()` methods,
  for clusters or single observations, with the one-step deletion estimates of
  Preisser and Qaqish (1996) or, for clusters, exact refits without each
  cluster. The one-step estimates for an observation condition on the rest of
  its cluster. Deletions that leave coefficients that cannot be estimated give
  `NA` with a warning.

* Added `local_influence()`, the total local influence and the direction of
  maximum curvature (Cook, 1986) for case-weight perturbations of the clusters
  or the observations and for additive perturbations of the responses.

* Added goodness-of-fit tests for binary responses. `gof_score_test()` gives
  the deciles-of-risk test of Horton et al. (1999) and the covariate-region
  test of Barnhart and Williamson (1998) as score tests at the fit.
  `gof_residual_test()` gives the Pearson and unweighted sum of squares tests
  of Pan (2002), and `gof_tests()` collects them in one table. The robust
  covariate-region statistic is the robust score statistic, which differs from
  the robust statistic printed by Barnhart and Williamson (1998). The tests
  were developed for `method = "gee"`, and fits with another method give a
  warning.

* Fitted objects store the cluster identifiers and the values of `repeated`
  supplied, as `cluster_labels` and `repeated_labels`. Cluster-level
  diagnostics are named by these identifiers, and observation-level
  diagnostics by `"id(repeated)"`.

* Added optional support for the `marginaleffects` package (predictions,
  comparisons and slopes).

* Expanded `geecriteria()` with `QICHH`, `QICC`, `EQIC`, `GHYC`, `PAC`, `AGPC`,
  `SGPC` and the eigenvalue-based `PT`, `WR` and `RMR` of Jang (2011), and a
  `criteria` argument to select which are computed.

* Added `runs_test()`, the Wald-Wolfowitz runs test of residual signs (Chang,
  2000, and Hardin and Hilbe, 2013), with a `plot()` method.

* Added `mcar_little_test()`, Little's (1988) MCAR test.

* Added `mcar_logistic_test()`, a Ridout (1991)-type logistic regression test of
  whether dropout depends on the previous response or on covariates, fitted with
  `geewa_binary()` and available with all five hypothesis tests.

* Added `mcar_homoscedasticity_test()`, the MCAR test of Jamshidian and Jalal
  (2010) as implemented in `MissMech`.

* Added cluster-level Mahalanobis residuals, `residuals(type = "mahalanobis")`,
  named by the cluster identifiers supplied.

## Bug fixes

* The heading of the `step_p()` table split long formulas over several pieces
  joined by runs of spaces. Formulas are now printed on one line.

* `step_p()` with a scope that does not nest the initial model failed with the
  `stats::factor.scope()` message "upper scope has term ... not included in
  model". `step_p()` and `step_criterion()` now name the scope and the missing
  terms.

* `test = "working-lrt"` scaled each model's working log-likelihood by its own
  estimated dispersion, so the statistic was not a likelihood ratio (for a
  Gaussian response it always equalled the number of tested parameters). Both
  models now use the dispersion of the larger model (Rotnitzky and Jewell,
  1990). Results change when the dispersion is estimated.

* The inverse-Gaussian quasi-log-likelihood had the wrong sign, affecting `QIC`,
  `QICu` and the working LRT for `inverse.gaussian` fits.

* P-values of the Wald, score and working tests were computed as `1 - pchisq()`
  and were exactly 0 for large statistics (above about 70 on one degree of
  freedom). They now use the upper tail.

* The bias-corrected covariance (Morel, Bokossa and Neerchal, 2003) is `NA`,
  with a warning, when the number of clusters does not exceed the number of
  parameters. Before, it could be negative or was silently capped.

* `opgee-jeffreys` and `hpgee-jeffreys` now estimate `alpha` and `phi` once, at
  the independence penalized estimate, instead of re-estimating them at the
  one-step estimate. Results change slightly.

* One-step estimators with `corstr = "fixed"` reset the supplied association
  parameters to zero in the second pass.

* The bias-reduced and penalized methods of `geewa_binary()` used a wrong
  variance derivative for clusters of size one.

* `anova()` dropped formula `offset()` terms from the null model, affecting
  score tests.

* Fits called as `geer::geewa(...)`, or with `phi_fixed` passed as a variable,
  were treated as odds-ratio fits or as having an estimated dispersion in tests,
  `geecriteria()`, `glance()` and `frechet_bounds_cor()`.

* A scalar `offset` failed with a `model.frame()` error and is now recycled.

* The `"unstructured"` and `"fixed"` structures failed when every cluster had a
  single observation.

* A negative test statistic (possible with the bias-corrected covariance in
  small samples) aborted `anova()`, `add1()`, `drop1()` and `step_p()`. It now
  gives `NA` with a warning.

* Refits in `add1()`, `drop1()`, `anova()` and `step_p()` reused a `beta_start`
  of the wrong length.

* `emmeans` support passed prior weights in the wrong row order for unsorted
  data.

* The design matrix of a one-column model lost its `assign` attribute and column
  name.

* `geecriteria()` used the whole deparsed fitted object as the row name of a
  model supplied through `do.call()`. Such models are now labelled `Model 1`,
  `Model 2` and so on.

* `anova()` with one model, `add1()`, `drop1()` and `step_p()` failed when the
  fit had omitted rows with missing values in a variable that a refit did not
  use, or when the data were not visible from the formula environment (for
  example a fit made inside a function). Refits now use the stored data and the
  rows of the original fit.

* `geewa()` did not check that the response lies in the support of the family.
  Negative responses for the Poisson-type families and non-positive responses
  for `Gamma` and `inverse.gaussian` gave a misleading starting-value error, or
  were fitted when `beta_start` was supplied. They now give the errors of
  `glm()`.

* The starting-value fit from `brglm2` emitted one "non-integer" warning per
  observation for Poisson responses and binomial proportions.

## Changes

* Default `cov_type` of `geecriteria()` is now `"robust"`.

* `cov_type` is validated in the same way by every function, requires an exact
  name, and its default is shown as a single value in the usage. Abbreviations
  such as `cov_type = "rob"`, which `vcov()`, `confint()`, `summary()`,
  `predict()` and `tidy()` accepted before, are now an error.

* Deviance residuals are now divided by the square root of the dispersion, as in
  `glmtoolbox`.

* `GESSC` counts only estimated association parameters, so its value changes for
  `corstr = "fixed"` and `orstr = "fixed"`.

* `print()` and `summary()` flag fits that did not converge.

* A converged binomial fit with fitted probabilities within 1e-12 of 0 or 1
  warns of possible separation.

* A numerical failure during the iterations no longer aborts the fit: the last
  accepted iterate is returned with a warning and `converged = FALSE`.
  Non-finite Newton steps are treated the same way.

* Refits in `anova()`, `add1()`, `drop1()` and `step_p()` stop when they do not
  converge.

* Nested comparisons stop when the fits differ in fitting function, method,
  dispersion handling or other settings. `anova()` with several models still
  omits, with a warning, a model with a different association structure.

* `geewa()` and `geewa_binary()` gain `subset` and `na.action`.

* `residuals()` gains `type = "response"`, an alias of `"working"`.

* `step_multiplier` in `geer_control()` may be any positive number, and values
  below 1 damp the step.

* Clearer errors for a factor, character or matrix response with non-binomial
  families, `beta_start` with `NA`, rank deficiency (naming the columns), an
  empty model matrix and numerical failures in the C++ code.

* Small dispersion estimates are no longer rejected, only zero or non-finite
  ones.

* Standard errors from a negative variance are `NA` with a warning.

* `step_p()`: the `Step` column is an integer, `direction = "forward"` without
  `scope` warns, and `direction = "both"` requires `p_enter <= p_remove` and
  stops instead of cycling.

* `geecriteria()` warns when models differ in family or response.

* `summary()` and `tidy()` stop on unused or misspelled arguments, as `vcov()`,
  `confint()` and the other methods already did.

* A character `repeated` variable gives a warning, because its values are
  ordered alphabetically (for example `"10"` before `"2"`).

* `glance()` computes only `QIC`, `QICu` and `CIC`.

* P-values below 0.0001 are printed as `<0.0001`.

* Faster C++ solvers with lower memory use. Estimates are unchanged.

## Documentation

* `geewa()` and `geewa_binary()` now state that `converged = TRUE` does not
  guarantee finite estimates under separation, that only `"brgee-naive"` and the
  Jeffreys-penalized methods protect against infinite estimates, and that the
  estimated dispersion of `geewa()` (the default) cancels this protection: use
  `phi_fixed = TRUE` or `geewa_binary()`.

* The documentation of the score tests states that `alpha` and `phi` are those
  of the larger model, not re-estimated under the null model.

* The AR(1) correlation estimator and the meaning of `Resid. Df` in `anova()`
  are documented.
