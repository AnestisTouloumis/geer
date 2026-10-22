fit_null <- geewa(
  formula = seizures ~ 1,
  data = test_data$epilepsy,
  id = id,
  family = poisson(link = "log"),
  corstr = "independence",
  method = "gee"
)

fit_trt <- geewa(
  formula = seizures ~ treatment,
  data = test_data$epilepsy,
  id = id,
  family = poisson(link = "log"),
  corstr = "independence",
  method = "gee"
)

test_that("add1 requires a non-empty scope", {
  expect_error(add1(fit_geewa_bin_exch), "'scope' contains no terms to add")
  expect_error(add1(fit_geewa_bin_exch, scope = NULL), "'scope' contains no terms to add")
})

test_that("add1 returns an anova-like table with the expected structure", {
  out <- add1(
    fit_geewa_bin_exch,
    scope = . ~ . + treatment:period,
    test = "score",
    cov_type = "robust"
  )
  expect_anova_table(
    out,
    candidate_pattern = "period:treatment"
  )
  expect_gte(nrow(out), 2L)
})

test_that("drop1 returns an anova-like table for explicit and default scope", {
  out_explicit <- drop1(
    fit_geewa_bin_exch,
    scope = "treatment",
    test = "score",
    cov_type = "robust"
  )
  expect_anova_table(
    out_explicit,
    candidate_pattern = "^treatment$",
    n_rows = 2L
  )
  out_default <- drop1(
    fit_geewa_bin_exch,
    test = "score",
    cov_type = "robust"
  )
  expect_anova_table(out_default)
  expect_gte(nrow(out_default), 2L)
})

test_that("add1 and drop1 require independence for working-lrt", {
  expect_error(
    add1(
      fit_geewa_bin_exch,
      scope = . ~ . + treatment:period,
      test = "working-lrt"
    ),
    "independence"
  )
  expect_error(
    drop1(
      fit_geewa_bin_exch,
      scope = "treatment",
      test = "working-lrt"
    ),
    "independence"
  )
})

test_that("anova single-model returns a sequential table with one row per term", {
  out <- anova(fit_geewa_pois_indep, test = "wald", cov_type = "robust")
  terms <- attr(fit_geewa_pois_indep$terms, "term.labels")
  expect_s3_class(out, "anova")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Df", "Resid. Df", "Chi", "Pr(>Chi)") %in% names(out)))
  expect_equal(nrow(out), 1L + length(terms))
  expect_identical(rownames(out), c("NULL", terms))
  expect_true(is.na(out$Chi[1L]))
  expect_true(is.na(out$`Pr(>Chi)`[1L]))
  expect_true(all(is.finite(out$Chi[-1L])))
  expect_true(all(out$Chi[-1L] >= 0))
  expect_true(all(out$`Pr(>Chi)`[-1L] >= 0 & out$`Pr(>Chi)`[-1L] <= 1))
})

test_that("anova heading contains family and link information", {
  out <- anova(fit_geewa_pois_indep, test = "wald", cov_type = "robust")
  heading <- attr(out, "heading")
  expect_type(heading, "character")
  expect_true(grepl("poisson", heading, ignore.case = TRUE))
  expect_true(grepl("log", heading, ignore.case = TRUE))
})

test_that("anova multi-model returns one row per model with decreasing residual df", {
  out <- anova(fit_null, fit_trt, fit_geewa_pois_indep, test = "wald", cov_type = "robust")
  expect_s3_class(out, "anova")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Df", "Resid. Df", "Chi", "Pr(>Chi)") %in% names(out)))
  expect_equal(nrow(out), 3L)
  expect_true(is.na(out$Chi[1L]))
  expect_true(is.na(out$`Pr(>Chi)`[1L]))
  expect_true(all(is.finite(out$Chi[-1L])))
  expect_true(all(out$Chi[-1L] >= 0))
  expect_true(all(out$`Pr(>Chi)`[-1L] >= 0 & out$`Pr(>Chi)`[-1L] <= 1))
  expect_true(all(diff(out$`Resid. Df`) < 0))
})

test_that("anova rejects working-lrt for non-independence models", {
  expect_error(
    anova(fit_trt, fit_geewa_pois_exch, test = "working-lrt", cov_type = "robust"),
    "independence"
  )
})

test_that("anova rejects non-geer and named objects in dots", {
  expect_error(
    anova(fit_null, fit_trt, "not_a_model", test = "wald", cov_type = "robust"),
    "anova.geer does not use the argument\\(s\\): <unnamed>"
  )
  expect_error(
    anova(fit_geewa_pois_indep, foo = "bar", test = "wald", cov_type = "robust"),
    "anova.geer does not use the argument\\(s\\): 'foo'"
  )
})

test_that("add1 and drop1 reject unused arguments", {
  expect_error(
    add1(fit_null, scope = ~ . + treatment, cov.type = "robust"),
    "add1.geer does not use the argument\\(s\\): 'cov.type'"
  )
  expect_error(
    drop1(fit_trt, tset = "score"),
    "drop1.geer does not use the argument\\(s\\): 'tset'"
  )
})

test_that("add1 accepts a scope when the formula was passed as a symbol", {
  user_formula <- seizures ~ treatment
  fit <- geewa(
    formula = user_formula,
    data = test_data$epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  out <- add1(fit, scope = ~ . + lnbaseline)
  expect_s3_class(out, "data.frame")
  expect_true("lnbaseline" %in% rownames(out))
})

test_that("anova and drop1 stop on a negative test statistic", {
  negative <- fit_geewa_pois_indep
  negative$robust_covariance <- -negative$robust_covariance
  smaller <- update(fit_geewa_pois_indep, . ~ . - lnage)
  expect_error(
    anova(smaller, negative, test = "wald", cov_type = "robust"),
    "Wald test failed: the test statistic is negative"
  )
  expect_error(
    drop1(negative, scope = "lnage", test = "wald", cov_type = "robust"),
    "Wald test failed: the test statistic is negative"
  )
})

test_that("add1 and drop1 pass pmethod to the working tests", {
  larger <- geewa(seizures ~ lnbaseline + factor(visit),
    data = test_data$epilepsy, id = id,
    repeated = visit, family = poisson(), corstr = "exchangeable"
  )
  smaller <- update(larger, . ~ . - factor(visit))
  p_values <- numeric(0)
  for (pmethod in c("rao-scott", "satterthwaite")) {
    expected <- compute_working_wald_test(smaller, larger,
      cov_type = "robust", pmethod = pmethod
    )
    dropped <- drop1(larger,
      scope = "factor(visit)", test = "working-wald",
      cov_type = "robust", pmethod = pmethod
    )
    added <- add1(smaller,
      scope = . ~ . + factor(visit), test = "working-wald",
      cov_type = "robust", pmethod = pmethod
    )
    expect_equal(dropped[["Pr(>Chi)"]][2], expected$test_p, tolerance = 1e-10)
    expect_equal(added[["Pr(>Chi)"]][2], expected$test_p, tolerance = 1e-10)
    p_values[pmethod] <- expected$test_p
  }
  expect_false(isTRUE(all.equal(p_values[["rao-scott"]], p_values[["satterthwaite"]])))
})

test_that("refits keep the rows of a fit with a missing covariate", {
  dat <- test_data$epilepsy
  dat$lnage[1:8] <- NA
  dat <- dat[c(seq(2, nrow(dat), by = 2), seq(1, nrow(dat), by = 2)), ]
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    family = poisson,
    data = dat,
    id = id,
    repeated = visit,
    corstr = "exchangeable"
  )
  complete <- dat[!is.na(dat$lnage), ]
  fit_reduced <- geewa(
    seizures ~ treatment + lnbaseline,
    family = poisson,
    data = complete,
    id = id,
    repeated = visit,
    corstr = "exchangeable"
  )
  for (test in c("wald", "score", "working-wald")) {
    expect_equal(
      drop1(fit, test = test)["lnage", "Chi"],
      anova(fit_reduced, fit, test = test)$Chi[2L],
      tolerance = 1e-10
    )
  }
  sequential <- anova(fit, test = "score")
  expect_equal(
    sequential["lnage", "Chi"],
    anova(fit_reduced, fit, test = "score")$Chi[2L],
    tolerance = 1e-10
  )
  expect_equal(
    add1(fit_reduced, scope = ~ . + lnage)["lnage", "Chi"],
    drop1(fit)["lnage", "Chi"],
    tolerance = 1e-10
  )
})

test_that("refits work for a fit made inside a function", {
  model_formula <- seizures ~ treatment + lnbaseline + lnage
  fit_in_function <- function(form, dat) {
    geewa(form, family = poisson, data = dat, id = id, repeated = visit)
  }
  early_visits <- test_data$epilepsy[test_data$epilepsy$visit < 4, ]
  fit <- fit_in_function(model_formula, early_visits)
  fit_reduced <- geewa(
    seizures ~ treatment + lnbaseline,
    family = poisson,
    data = early_visits,
    id = id,
    repeated = visit
  )
  expect_equal(
    drop1(fit)["lnage", "Chi"],
    anova(fit_reduced, fit)$Chi[2L],
    tolerance = 1e-10
  )
  expect_no_error(anova(fit))
})

test_that("refits keep the user-facing call and the original subset rows", {
  weights <- rep(c(1, 2), length.out = nrow(test_data$epilepsy))
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    family = poisson,
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    weights = weights,
    subset = visit > 1
  )
  refit <- refit_geer(fit, . ~ . - lnage)
  expect_identical(refit$call$data, quote(test_data$epilepsy))
  expect_identical(refit$call$subset, quote(visit > 1))
  expect_identical(refit$obs_no, fit$obs_no)
  expect_identical(refit$frame_rows, fit$frame_rows)
  reference <- geewa(
    seizures ~ treatment + lnbaseline,
    family = poisson,
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    weights = weights,
    subset = visit > 1
  )
  expect_equal(coef(refit), coef(reference), tolerance = 1e-10)
})

test_that("anova of a model without intercept tests the terms sequentially", {
  fit <- geewa(
    formula = seizures ~ 0 + treatment + lnbaseline,
    family = poisson(),
    data = test_data$epilepsy,
    id = id,
    repeated = visit
  )
  tab <- anova(fit)
  expect_identical(rownames(tab), c("treatment", "lnbaseline"))
  smaller <- update(fit, formula = . ~ . - lnbaseline)
  expect_equal(
    tab["lnbaseline", "Chi"],
    anova(smaller, fit)$Chi[2L],
    tolerance = 1e-10
  )
})

test_that("drop1 accepts a formula scope", {
  fit <- geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    family = poisson(),
    data = test_data$epilepsy,
    id = id,
    repeated = visit
  )
  expect_equal(drop1(fit, scope = ~lnage), drop1(fit, scope = "lnage"))
  expect_error(
    drop1(fit, scope = ~visit),
    "'scope' must contain only terms of the model"
  )
})
