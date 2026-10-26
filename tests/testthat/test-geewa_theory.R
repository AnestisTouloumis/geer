list2env(get_cached_test_fixtures("geewa_theory"), environment())

test_that("geewa under independence matches the corresponding GLM", {
  skip_on_cran()
  expect_equal(coef(fit_gee_probit), coef(fit_glm_probit), tolerance = 1e-8)
  expect_equal(
    vcov(fit_gee_probit, cov_type = "naive"),
    vcov(fit_glm_probit),
    tolerance = 1e-8
  )
  expect_equal(coef(fit_gee_logit), coef(fit_glm_logit), tolerance = 1e-8)
  expect_equal(
    vcov(fit_gee_logit, cov_type = "naive"),
    vcov(fit_glm_logit),
    tolerance = 1e-8
  )
})

test_that("brgee-naive matches brglm2 AS_mean for representative links", {
  skip_on_cran()
  fit_glm_as_mean_probit <- fit_brglm(fit_glm_probit, "AS_mean")
  fit_glm_as_mean_logit <- fit_brglm(fit_glm_logit, "AS_mean")
  expect_equal(
    coef(fit_brgee_naive_probit),
    coef(fit_glm_as_mean_probit),
    tolerance = 1e-8
  )
  expect_equal(
    vcov(fit_brgee_naive_probit, cov_type = "naive"),
    vcov(fit_glm_as_mean_probit),
    tolerance = 1e-8
  )
  expect_equal(
    coef(fit_brgee_naive_logit),
    coef(fit_glm_as_mean_logit),
    tolerance = 1e-8
  )
  expect_equal(
    vcov(fit_brgee_naive_logit, cov_type = "naive"),
    vcov(fit_glm_as_mean_logit),
    tolerance = 1e-8
  )
})

test_that("pgee-jeffreys matches brgee-naive for logit but not for probit", {
  skip_on_cran()
  expect_equal(
    coef(fit_brgee_naive_logit),
    coef(fit_pgee_jeffreys_logit),
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(
    coef(fit_brgee_naive_probit),
    coef(fit_pgee_jeffreys_probit)
  )))
})

test_that("use_p = FALSE reproduces gee::gee 4.13.26 (ar1, unstructured)", {
  fit_ar1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    family = poisson(),
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    corstr = "ar1",
    use_p = FALSE,
    control = geer_control(tolerance = 1e-10)
  )
  fit_unstr <- update(fit_ar1, corstr = "unstructured")
  expect_equal(
    unname(coef(fit_ar1)),
    c(-2.6467004621238, -0.0234155917577, 1.2452149614426, 0.6446043492813),
    tolerance = 1e-6
  )
  expect_equal(unname(fit_ar1$alpha), 0.481015591440, tolerance = 1e-6)
  expect_equal(
    unname(coef(fit_unstr)),
    c(-2.5133346699774, -0.0270451301904, 1.2336125872299, 0.6120666013723),
    tolerance = 1e-6
  )
  expect_equal(
    unname(fit_unstr$alpha),
    c(
      0.470858677472, 0.404341082539, 0.222763170694,
      0.577551171276, 0.292495071726, 0.392603234262
    ),
    tolerance = 1e-6
  )
})

test_that("pgee-jeffreys matches brglm2 Jeffreys fits for representative links", {
  skip_on_cran()
  fit_glm_jeffreys_logit <- fit_brglm(fit_glm_logit, "MPL_Jeffreys")
  fit_glm_jeffreys_probit <- fit_brglm(fit_glm_probit, "MPL_Jeffreys")
  expect_equal(
    coef(fit_pgee_jeffreys_logit),
    coef(fit_glm_jeffreys_logit),
    tolerance = 1e-8
  )
  expect_equal(
    vcov(fit_pgee_jeffreys_logit, cov_type = "naive"),
    vcov(fit_glm_jeffreys_logit),
    tolerance = 1e-8
  )
  expect_equal(
    coef(fit_pgee_jeffreys_probit),
    coef(fit_glm_jeffreys_probit),
    tolerance = 1e-8
  )
  expect_equal(
    vcov(fit_pgee_jeffreys_probit, cov_type = "naive"),
    vcov(fit_glm_jeffreys_probit),
    tolerance = 1e-8
  )
})
