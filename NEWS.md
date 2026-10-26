# geer 0.1.1

## New features

* Added a `sigma()` method returning the square root of a fitted model's
  working dispersion.

* Added `cov_type = "jackknife"`, the leave-one-cluster-out jackknife covariance
  (each estimate a full refit with the association parameters held at their
  full-data values), available for all estimation methods and wherever
  `cov_type` is accepted. It is computed once and cached on the fit.

* Added nine small-sample sandwich estimators to `cov_type`, those of
  Kauermann and Carroll (2001), Mancl and DeRouen (2001), Ford and Westgate
  (2017), Fay and Graubard (2001), Fan, Zhang and Zhang (2013), Pan (2001),
  Gosho, Sato and Takeuchi (2014), Wang and Long (2011), and Westgate and
  Burchett (2016), for `geewa()` and `geewa_binary()` fits with any family and
  link. They are available wherever `cov_type` is accepted, including the score
  tests, where they are evaluated at the fit of the smaller model.

* Added `step_criterion()`, stepwise model selection that moves to the
  addition or deletion with the smallest QICu, QIC, QICC, AGPC, or SGPC. The
  candidate table of every step is returned, and printed with `trace = TRUE`.

* Added `hatvalues()`, `dfbeta()`, `dfbetas()`, and `cooks.distance()` methods,
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
  supplied, as `cluster_labels` and `repeated_labels`, and also
  `fit_function`, `use_p`, and `phi_fixed`. Cluster-level diagnostics are
  named by these identifiers, and observation-level diagnostics by
  `"id(repeated)"`.

* Added optional support for the `marginaleffects` package (predictions,
  comparisons, and slopes).

* Expanded `geecriteria()` with `QICHH`, `QICC`, `EQIC`, `GHYC`, `PAC`, `AGPC`,
  `SGPC`, and the eigenvalue-based `PT`, `WR`, and `RMR` of Jang (2011), and a
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

* `geewa()` and `geewa_binary()` gain `subset` and `na.action`.

* `residuals()` gains `type = "response"`, an alias of `"working"`.

## Bug fixes

* Fixed Gaussian `emmeans` integration and reference-grid recovery after
  missing observations or subsetting. Reference grids now use the fitting rows
  with their corresponding prior weights.

* Improved the Plackett covariance and derivative calculations for
  near-complementary marginal probabilities and extreme negative association
  by retaining the rounding residual in the discriminant.

* Used scale-relative symmetry checks and general linear solves for the
  nonsymmetric systems in cluster deletion diagnostics.

* Corrected cluster error context in association estimation and pooled
  covariance accumulation.

* Rejected overflowing integer controls, incomplete family specifications,
  nonnumeric fixed dispersion, and non-finite marginaleffects inputs.

* Applied the shared variance checks to MCAR coefficient tables and improved
  prediction and MCAR reconstruction input errors.

* Rejected unused print arguments while retaining supported formatting options.

* Corrected the plain-text alternative of the bias-correction equation in help.

* Preserved Plackett association near probability boundaries and avoided
  overflow for extreme finite odds ratios, including covariance derivatives.

* Corrected the local influence direction for rank-deficient perturbations.
  A zero-rank perturbation now gives an error for the local direction.

* Corrected Anderson-Darling counts for ties at the pooled maximum and
  rejected samples with only one distinct transformed value.

* Centered and scaled normal-theory MCAR EM calculations to stabilize
  estimation under translations and changes of measurement units.

* Required agreement of shared design spaces in nested-model comparisons.

* Assigned implicit repeated positions before subsetting and missing-value
  removal, preserving gaps in fitting and MCAR data reconstruction.

* The heading of the `step_p()` table split long formulas over several pieces
  joined by runs of spaces. Formulas are now printed on one line.

* `step_p()` with a scope that does not nest the initial model failed with the
  `stats::factor.scope()` message "upper scope has term ... not included in
  model". `step_p()` and `step_criterion()` now name the scope and the missing
  terms.

* `test = "working-lrt"` scaled each model's working log-likelihood by its own
  estimated dispersion, so the statistic was not a likelihood ratio (for a
  Gaussian response it always equaled the number of tested parameters). Both
  models now use the dispersion of the larger model (Rotnitzky and Jewell,
  1990). Results change when the dispersion is estimated.

* The inverse-Gaussian quasi-log-likelihood had the wrong sign, affecting `QIC`,
  `QICu`, and the working LRT for `inverse.gaussian` fits.

* P-values of the Wald, score, and working tests were computed as
  `1 - pchisq()` and were exactly 0 for large statistics (above about 70 on one
  degree of freedom). They now use the upper tail.

* The bias-corrected covariance (Morel, Bokossa and Neerchal, 2003) is
  undefined when the number of clusters does not exceed the number of
  parameters: the fit warns, and `vcov()`, the methods that use it, and the
  score tests give an error. Before, it could be negative or was silently
  capped.

* `opgee-jeffreys` and `hpgee-jeffreys` now estimate `alpha` and `phi` once, at
  the independence penalized estimate, instead of re-estimating them at the
  one-step estimate. Results change slightly.

* One-step estimators with `corstr = "fixed"` set the supplied association
  parameters to zero in the second pass. They now keep them.

* The bias-reduced and penalized methods of `geewa_binary()` used a wrong
  variance derivative for clusters of size one.

* `anova()` dropped formula `offset()` terms from the null model, affecting
  score tests.

* Fits called as `geer::geewa(...)`, or with `phi_fixed` passed as a variable,
  were treated as odds-ratio fits or as having an estimated dispersion in tests,
  `geecriteria()`, `glance()`, and `frechet_bounds_cor()`.

* A scalar `offset` failed with a `model.frame()` error and is now recycled.

* The `"unstructured"` and `"fixed"` structures failed when every cluster had a
  single observation.

* A negative test statistic, possible with the bias-corrected covariance in
  small samples, aborted `anova()`, `add1()`, `drop1()`, and `step_p()`. They
  now stop with an error that names the test and says that the covariance
  estimate is not positive semidefinite.

* Refits in `add1()`, `drop1()`, `anova()`, and `step_p()` reused a
  `beta_start` of the wrong length.

* `emmeans` support passed prior weights in the wrong row order for unsorted
  data.

* The design matrix of a one-column model lost its `assign` attribute and column
  name.

* `geecriteria()` used the whole deparsed fitted object as the row name of a
  model supplied through `do.call()`. Such models are now labeled `Model 1`,
  `Model 2`, and so on.

* `anova()` with one model, `add1()`, `drop1()`, and `step_p()` failed when the
  fit had omitted rows with missing values in a variable that a refit did not
  use, or when the data were not visible from the formula environment (for
  example a fit made inside a function). Refits now use the stored data and the
  rows of the original fit.

* `geewa()` did not check that the response lies in the support of the family.
  Negative responses for the Poisson-type families and nonpositive responses
  for `Gamma` and `inverse.gaussian` gave a misleading starting-value error, or
  were fitted when `beta_start` was supplied. They now give the errors of
  `glm()`.

* The starting-value fit from `brglm2` emitted one "non-integer" warning per
  observation for Poisson responses and binomial proportions.

* Fitted probabilities near 0 under the `cloglog` link lost relative accuracy
  through cancellation in the inverse link.

* `predict()` with `se.fit = TRUE` warned that the variance was negative or
  nonfinite for rows of `newdata` with missing covariates. Such rows now give
  `NA` without a warning.

* `mcar_logistic_test()` failed with "the binary GEE missingness model did not
  converge" when an occasion had no dropouts, as after monotone dropout that
  ends before the last occasion. Such occasions, and those with no observed
  responses, are now omitted from the risk set with a warning.

## Changes

* Strengthened scalar checks against rounded published values, independent
  model-selection references, subprocess isolation, and regression coverage.
  Removed redundant smoke checks and reduced repeated test fitting.

* Integration examples now print both results. Examples for test printers cover
  both result classes, and the runs plot shows how to omit cluster boundaries.

* Corrected the documented MCAR imputation policy and data precedence, and
  clarified numerical limitations of the Plackett covariance and information
  criteria. Updated cluster terminology, narrative citation punctuation, and
  documentation section order. These documentation changes do not change
  estimates, tests, or runtime behavior.

* Shared cluster data construction between the correlation and odds-ratio
  backends and normalized internal score vectors.

* Centralized fitting choice matching with exact character names and retained
  strict NULL handling. Method factors are no longer coerced to character.

* Documented the retained inactive-argument, case-insensitive criterion-name,
  and integration-argument exceptions.

* Standardized scalar vapply templates and negated membership syntax.

* Reconciled the scoring signs of the bias-corrected estimators with
  Touloumis (2026a), including the empirical sensitivity convention and the
  first contraction of its empirical adjustment. Numerical bias correction
  is unchanged.

* `add1()`, `drop1()`, `model.matrix()`, `glance()`, and the `print()` methods
  for fits and their summaries now reject unused arguments. `anova()` rejects
  named arguments and objects that are not `geer` fits instead of ignoring
  them.

* `geewa_binary()` checks `link` against the supported links before building
  the family, and `geewa()` checks that `beta_start`, `alpha_vector`, and
  `control_glm` have the right type before using them.

* The Kauermann-Carroll, Mancl-DeRouen, Ford-Westgate, Fan-Zhang-Zhang,
  Wang-Long, and Westgate-Burchett covariances now stop with an error that
  names the clusters with leverage 1, for which I - H is singular, instead of
  returning a covariance from a numerically singular inverse. The other
  small-sample covariances remain available for such fits. All small-sample
  covariances stop when the naive covariance is not positive definite.

* A working covariance or information matrix that is symmetric but not positive
  definite, such as an odds-ratio working covariance built from incompatible
  pairwise odds ratios, is still solved by LU decomposition, but now with a
  warning.

* Default `cov_type` of `geecriteria()` is now `"robust"`.

* `geecriteria()` returns the criteria unrounded, as an object of class
  `geer_criteria`, and `digits` sets the decimal places only when printing.

* `cov_type` is validated in the same way by every function, requires an exact
  name, and its default is shown as a single value in the usage. Abbreviations
  such as `cov_type = "rob"`, which `vcov()`, `confint()`, `summary()`,
  `predict()`, and `tidy()` accepted before, are now an error. The same holds
  for `vcov.method` in the `emmeans` support, which also accepted other
  capitalizations and separators (for example `"Bias.Corrected"`).

* In `geewa()` and `geewa_binary()`, `...` builds only `control`, and
  `control_glm` defaults to `list()`. Before, `...` also built `control_glm`:
  `brglm2` settings in `...` failed in `geer_control()`, and a partial name
  such as `maxit` silently set `maxiter`. Arguments of `geer_control()` must
  now be given by their full names.

* Deviance residuals are now divided by the square root of the dispersion, as in
  `glmtoolbox`.

* `GESSC` counts only estimated association parameters, so its value changes for
  `corstr = "fixed"` and `orstr = "fixed"`.

* `print()` and `summary()` flag fits that did not converge.

* A converged binomial fit with fitted probabilities within 1e-12 of 0 or 1
  warns of possible separation.

* A numerical failure during the iterations no longer aborts the fit: the last
  accepted iterate is returned with a warning and `converged = FALSE`.
  Nonfinite Newton steps are treated the same way.

* Refits in `anova()`, `add1()`, `drop1()`, and `step_p()` stop when they do
  not converge.

* Nested comparisons stop when the fits differ in fitting function, method,
  dispersion handling, or other settings. `anova()` with several models still
  omits, with a warning, a model with a different association structure.

* `step_multiplier` in `geer_control()` may be any positive number, and values
  below 1 damp the step.

* Clearer errors for a factor, character, or matrix response with nonbinomial
  families, `beta_start` with `NA`, rank deficiency (naming the columns), an
  empty model matrix, and numerical failures in the C++ code.

* Errors from numerical failures start with the name of the function that
  failed, for example `anova.geer failed:` or `vcov.geer failed:`, and errors
  in the arguments name the argument.

* `predict()` checks that `newdata` contains the variables of the offset
  before building the model frame, so a missing offset variable gives this
  error instead of an error from `model.frame()`.

* Small dispersion estimates are no longer rejected, only zero or nonfinite
  ones.

* `summary()`, `tidy()`, and `predict()` stop when a variance is negative or
  not finite, as `confint()` already did.

* `step_p()`: the `Step` column is an integer, `direction = "forward"` requires
  `scope`, and `direction = "both"` requires `p_enter <= p_remove` and stops
  instead of cycling.

* `geecriteria()` stops when the models differ in the number of
  observations, the family, or the response.

* `summary()` and `tidy()` stop on unused or misspelled arguments, as `vcov()`,
  `confint()`, and the other methods already did.

* A character `repeated` variable gives an error, because it has no natural
  order (sorted alphabetically, `"10"` would precede `"2"`).

* `glance()` computes only `QIC`, `QICu`, and `CIC`.

* `anova()`, `add1()`, `drop1()`, and the `anova` component of `step_p()`
  return objects of class `geer_anova`, which inherits from `anova`. P-values
  below 0.0001 are printed as `<0.0001`.

* Faster C++ solvers with lower memory use. Estimates are unchanged.

## Documentation

* `geewa()` and `geewa_binary()` now state that `converged = TRUE` does not
  guarantee finite estimates under separation. The Jeffreys-penalized methods
  give finite estimates for binary responses under working independence
  (Touloumis, 2026b), with `phi_fixed = TRUE` in `geewa()` or with
  `geewa_binary()`, because the estimated dispersion of `geewa()` (the
  default) removes this protection.

* The documentation of the score tests states that `alpha` and `phi` are those
  of the larger model, not re-estimated under the null model.

* The AR(1) correlation estimator and the meaning of `Resid. Df` in `anova()`
  are documented.

* `?geewa` states that the single-stage methods update `alpha` and `phi`
  jointly with the regression parameters, and that its Jeffreys-penalized
  methods extend Touloumis (2026b), which covers binary responses with
  odds-ratio structures.
