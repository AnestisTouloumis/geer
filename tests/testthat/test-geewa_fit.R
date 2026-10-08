testthat::local_edition(3)

fit_geewa_gauss_indep <- geewa(
  formula = score ~ treatment + baseline + time,
  data = test_data$rinse,
  id = id,
  family = gaussian(link = "identity"),
  corstr = "independence",
  method = "gee"
)


test_that("geewa fits a coherent Gaussian independence model", {
  fit <- fit_geewa_gauss_indep
  expect_s3_class(fit, "geer")
  expect_true(fit$converged)
  expect_equal(fit$family$family, "gaussian")
  expect_equal(fit$family$link, "identity")
  expect_equal(names(coef(fit)), colnames(fit$x))
  expect_true(all(is.finite(fitted(fit))))
  expect_equal(length(fitted(fit)), fit$obs_no)
  expect_equal(
    residuals(fit, type = "working"),
    fit$y - fitted(fit),
    tolerance = 1e-10
  )
  expect_gt(fit$phi, 0)
})


test_that("geewa fits coherent Poisson independence and exchangeable models", {
  expect_true(fit_geewa_pois_indep$converged)
  expect_true(all(fitted(fit_geewa_pois_indep) > 0))
  expect_true(fit_geewa_pois_exch$converged)
  expect_true(all(fitted(fit_geewa_pois_exch) > 0))
  expect_length(fit_geewa_pois_exch$alpha, 1L)
  expect_gt(fit_geewa_pois_exch$alpha, -1)
  expect_lt(fit_geewa_pois_exch$alpha, 1)
  expect_equal(
    fit_geewa_pois_indep$clusters_no,
    length(unique(test_data$epilepsy$id))
  )
  expect_equal(
    fit_geewa_pois_indep$obs_no,
    nrow(test_data$epilepsy)
  )
})


test_that("geewa fits representative non-independence correlation structures", {
  fit_ar1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "ar1"
  )
  fit_unstr <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "unstructured"
  )
  fit_mdep <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "m-dependent",
    Mv = 1
  )
  fit_fixed <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "fixed",
    alpha_vector = rep(0.2, choose(4, 2))
  )
  expect_s3_class(fit_ar1, "geer")
  expect_true(fit_ar1$converged)
  expect_length(fit_ar1$alpha, 1L)
  expect_s3_class(fit_unstr, "geer")
  expect_true(fit_unstr$converged)
  expect_length(fit_unstr$alpha, choose(4, 2))
  expect_s3_class(fit_mdep, "geer")
  expect_true(fit_mdep$converged)
  expect_s3_class(fit_fixed, "geer")
  expect_true(fit_fixed$converged)
})


test_that("geewa respects phi_fixed and phi_value", {
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence",
    phi_fixed = TRUE,
    phi_value = 2
  )
  expect_s3_class(fit, "geer")
  expect_equal(fit$phi, 2, tolerance = 1e-10)
})


test_that("geewa converges for representative alternative methods", {
  methods <- c("brgee-robust", "bcgee-robust", "pgee-jeffreys")
  for (method_name in methods) {
    fit <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "exchangeable",
      method = method_name
    )
    expect_s3_class(fit, "geer")
    expect_true(fit$converged)
    expect_true(all(is.finite(coef(fit))))
  }
})


test_that("bias-reduced estimates differ from plain GEE estimates", {
  fit_br <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "exchangeable",
    method = "brgee-robust"
  )
  expect_false(isTRUE(all.equal(
    coef(fit_geewa_pois_exch),
    coef(fit_br),
    tolerance = 1e-12
  )))
  ## Bias reduction should move each estimate by only a small fraction of its
  ## standard error with 59 clusters; a relative tolerance on the coefficients
  ## would say nothing about the size of the shift.
  shift_in_se <- abs(coef(fit_br) - coef(fit_geewa_pois_exch)) /
    sqrt(diag(vcov(fit_geewa_pois_exch, cov_type = "robust")))
  expect_true(all(shift_in_se < 0.5))
})


test_that("geewa is invariant to row order via internal sorting", {
  set.seed(1)
  fit_1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  epilepsy_shuffled <- test_data$epilepsy[
    sample.int(nrow(test_data$epilepsy)),
    ,
    drop = FALSE
  ]
  fit_2 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy_shuffled,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_equal(
    unname(fit_1$coefficients),
    unname(fit_2$coefficients),
    tolerance = 1e-8
  )
  expect_equal(names(fit_1$coefficients), names(fit_2$coefficients))
})


testthat::test_that("grouped binomial fit does not fail from an n - p^2 phi denominator", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1),
    id = 1:4
  )
  testthat::expect_no_error(
    fit <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      use_p = TRUE
    )
  )
  testthat::expect_s3_class(fit, "geer")
})


testthat::test_that("use_p = TRUE and use_p = FALSE both fit for a small grouped binomial model", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1),
    id = 1:4
  )
  testthat::expect_no_error(
    fit_true <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      use_p = TRUE
    )
  )
  testthat::expect_no_error(
    fit_false <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      use_p = FALSE
    )
  )
  testthat::expect_s3_class(fit_true, "geer")
  testthat::expect_s3_class(fit_false, "geer")
})


test_that("a fixed working correlation works with the independence-start one-step methods", {
  alpha_fixed <- rep(0.2, choose(4, 2))
  for (method in c("opgee-jeffreys", "hpgee-jeffreys")) {
    fit <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "fixed",
      alpha_vector = alpha_fixed,
      method = method
    )
    expect_s3_class(fit, "geer")
    expect_true(all(is.finite(coef(fit))))
    expect_equal(as.numeric(fit$alpha), alpha_fixed)
  }
})


test_that("geewa honours 'subset' and 'na.action' and rejects stray arguments", {
  epi <- test_data$epilepsy
  epi_sub <- epi[as.numeric(epi$id) <= 30, , drop = FALSE]
  fit_subset <- geewa(
    seizures ~ treatment + lnbaseline,
    data = epi,
    id = id,
    family = poisson("log"),
    subset = as.numeric(id) <= 30
  )
  fit_manual <- geewa(
    seizures ~ treatment + lnbaseline,
    data = epi_sub,
    id = id,
    family = poisson("log")
  )
  expect_equal(fit_subset$obs_no, nrow(epi_sub))
  expect_equal(coef(fit_subset), coef(fit_manual), tolerance = 1e-8)
  ## 'subset' must not be ignored when 'control' and 'control_glm' are given
  fit_explicit <- geewa(
    seizures ~ treatment + lnbaseline,
    data = epi,
    id = id,
    family = poisson("log"),
    control = geer_control(),
    control_glm = list(),
    subset = as.numeric(id) <= 30
  )
  expect_equal(fit_explicit$obs_no, nrow(epi_sub))
  epi_na <- epi
  epi_na$lnbaseline[1L] <- NA
  fit_na <- geewa(
    seizures ~ treatment + lnbaseline,
    data = epi_na,
    id = id,
    family = poisson("log")
  )
  expect_equal(fit_na$obs_no, nrow(epi) - 1L)
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline,
      data = epi_na,
      id = id,
      family = poisson("log"),
      na.action = stats::na.fail
    )
  )
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline,
      data = epi,
      id = id,
      family = poisson("log"),
      control = geer_control(),
      control_glm = list(),
      not_an_argument = TRUE
    ),
    "does not use the argument"
  )
})


test_that("one-step penalized methods hold phi at the independence penalized fit", {
  epi <- test_data$epilepsy
  fit_ind <- geewa(
    seizures ~ treatment + lnbaseline,
    data = epi,
    id = id,
    family = poisson("log"),
    corstr = "independence",
    method = "pgee-jeffreys"
  )
  for (method in c("opgee-jeffreys", "hpgee-jeffreys")) {
    fit_one <- geewa(
      seizures ~ treatment + lnbaseline,
      data = epi,
      id = id,
      family = poisson("log"),
      corstr = "exchangeable",
      method = method
    )
    expect_true(fit_one$converged)
    expect_equal(fit_one$phi, fit_ind$phi, tolerance = 1e-8)
  }
})


test_that("the fit stores the row order and residuals accept type = 'response'", {
  epi <- test_data$epilepsy
  set.seed(2026)
  epi_shuffled <- epi[sample(nrow(epi)), , drop = FALSE]
  fit <- geewa(
    seizures ~ treatment + lnbaseline,
    data = epi_shuffled,
    id = id,
    family = poisson("log")
  )
  expect_equal(sort(fit$row_order), seq_len(nrow(epi_shuffled)))
  expect_equal(
    fit$id[order(fit$row_order)],
    as.numeric(factor(epi_shuffled$id))
  )
  expect_equal(
    residuals(fit, type = "response"),
    residuals(fit, type = "working")
  )
})


test_that("geewa fits a response on a very small scale (dispersion below machine epsilon)", {
  set.seed(101)
  n_id <- 30L
  dat <- data.frame(id = rep(seq_len(n_id), each = 3L))
  dat$x <- rnorm(nrow(dat))
  dat$y_unit <- 1 + dat$x + rnorm(nrow(dat))
  dat$y_small <- 1e-9 * dat$y_unit
  for (corstr in c("independence", "exchangeable")) {
    fit_unit <- geewa(
      y_unit ~ x, data = dat, id = id,
      family = gaussian("identity"), corstr = corstr
    )
    expect_no_error(
      fit_small <- geewa(
        y_small ~ x, data = dat, id = id,
        family = gaussian("identity"), corstr = corstr
      )
    )
    expect_equal(unname(coef(fit_small)), unname(1e-9 * coef(fit_unit)),
                 tolerance = 1e-6)
    expect_equal(fit_small$phi, 1e-18 * fit_unit$phi, tolerance = 1e-6)
    expect_equal(fit_small$alpha, fit_unit$alpha, tolerance = 1e-6)
  }
})


test_that("an unstructured correlation structure works when every cluster has one observation", {
  set.seed(102)
  dat <- data.frame(id = seq_len(40L), x = rnorm(40L))
  dat$y <- 1 + dat$x + rnorm(40L)
  fit_ind <- geewa(
    y ~ x, data = dat, id = id,
    family = gaussian("identity"), corstr = "independence"
  )
  expect_no_error(
    fit_unstr <- geewa(
      y ~ x, data = dat, id = id,
      family = gaussian("identity"), corstr = "unstructured"
    )
  )
  expect_length(fit_unstr$alpha, 0L)
  expect_equal(coef(fit_unstr), coef(fit_ind), tolerance = 1e-6)
})


test_that("the standard-correlation solver stops on a non-finite starting step", {
  ## A NaN response makes the first Newton step non-finite; the solver must
  ## signal an error instead of returning NaN estimates.
  model_matrix <- cbind(1, c(-1, 0, 1, 2))
  expect_error(
    fit_geesolver_cc(
      y_vector = c(1, NaN, 2, 3), model_matrix = model_matrix,
      id_vector = c(1, 1, 2, 2), repeated_vector = c(1, 2, 1, 2),
      weights_vector = rep(1, 4), link = "identity", family = "gaussian",
      beta_vector = c(0, 0), offset = rep(0, 4), maxiter = 5L,
      tolerance = 1e-6, step_maxiter = 3L, step_multiplier = 1,
      jeffreys_power = 0.5, method = "gee", use_params = 0L,
      alpha_vector = 0, alpha_fixed = 1L,
      correlation_structure = "independence", mdependence = 1L,
      phi = 1, phi_fixed = 1L, hold_nuisance = 0L
    )
  )
})

