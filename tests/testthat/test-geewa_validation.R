epilepsy <- test_data$epilepsy

test_that("geewa rejects unsupported correlation structures and methods", {
  expect_error(
    geewa(
      seizures ~ treatment,
      data = epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "compound"
    ),
    "'corstr' must be one of"
  )
  expect_error(
    geewa(
      seizures ~ treatment,
      data = epilepsy,
      id = id,
      family = poisson("log"),
      method = "em-gee"
    ),
    "'method' must be one of"
  )
})

test_that("geewa validates m-dependent and fixed-correlation arguments", {
  expect_error(
    geewa(
      seizures ~ treatment,
      data = epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "m-dependent",
      Mv = 0
    ),
    "'Mv' must be a positive integer"
  )
  expect_error(
    geewa(
      seizures ~ treatment,
      data = epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "fixed"
    ),
    "'alpha_vector' must be provided"
  )
})

test_that("geewa validates beta_start", {
  expect_error(
    geewa(
      seizures ~ treatment,
      data = epilepsy,
      id = id,
      family = poisson("log"),
      beta_start = c(0, 0, 0)
    ),
    "'beta_start' must be a numeric vector of length"
  )
})

test_that("geewa rejects unsupported quasi variance functions early", {
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    x = c(0, 1, 0, 1),
    id = 1:4
  )
  expect_error(
    geewa(
      y ~ x,
      data = dat,
      id = id,
      family = quasi(link = "identity", variance = "mu^4"),
      corstr = "independence"
    ),
    "'variance'"
  )
})

test_that("geewa rejects duplicated repeated values within cluster", {
  epilepsy_bad <- epilepsy
  epilepsy_bad$rep2 <- epilepsy_bad$visit
  idx <- which(epilepsy_bad$id == epilepsy_bad$id[1])[1:2]
  epilepsy_bad$rep2[idx] <- 1
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy_bad,
      id = id,
      repeated = rep2,
      family = poisson(link = "log"),
      corstr = "independence"
    ),
    "distinct",
    ignore.case = TRUE
  )
})

test_that("geewa accepts row-level prior weights for two-column binomial responses", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1),
    id = 1:4,
    w = c(1, 2, 1, 3)
  )
  expect_no_error(
    fit <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = binomial(link = "logit"),
      corstr = "independence",
      weights = w,
      phi_fixed = TRUE,
      phi_value = 1
    )
  )
  expect_length(fit$prior.weights, nrow(dat))
  expect_equal(
    unname(fit$prior.weights),
    unname(dat$w * (dat$success + dat$failure))
  )
})

test_that("extract_geer_response_weights treats quasibinomial grouped responses like binomial", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1)
  )
  mf <- model.frame(
    as.formula(cbind(success, failure) ~ x),
    data = dat
  )
  out_binomial <- extract_geer_response_weights(
    mf,
    binomial(link = "logit")
  )
  out_quasibinomial <- extract_geer_response_weights(
    mf,
    quasibinomial(link = "logit")
  )
  expect_equal(out_quasibinomial$y, out_binomial$y)
  expect_equal(out_quasibinomial$weights, out_binomial$weights)
})

test_that("extract_geer_response_weights treats quasibinomial factor responses like binomial", {
  dat <- data.frame(
    y = factor(c("no", "yes", "no", "yes")),
    x = c(0, 1, 0, 1)
  )
  mf <- model.frame(
    as.formula(y ~ x),
    data = dat
  )
  out_binomial <- extract_geer_response_weights(
    mf,
    binomial(link = "logit")
  )
  out_quasibinomial <- extract_geer_response_weights(
    mf,
    quasibinomial(link = "logit")
  )
  expect_equal(out_quasibinomial$y, out_binomial$y)
  expect_equal(out_quasibinomial$weights, out_binomial$weights)
})

test_that("geewa rejects non-finite beta_start and names aliased columns", {
  epi <- epilepsy
  expect_error(
    geewa(seizures ~ treatment + lnbaseline,
      data = epi, id = id,
      family = poisson("log"), beta_start = c(NA, 0, 0)
    ),
    "'beta_start' must contain only finite values"
  )
  expect_error(
    geewa(seizures ~ treatment + lnbaseline,
      data = epi, id = id,
      family = poisson("log"), beta_start = c(0, Inf, 0)
    ),
    "'beta_start' must contain only finite values"
  )
  epi$lnbaseline2 <- 2 * epi$lnbaseline
  expect_error(
    geewa(seizures ~ lnbaseline + lnbaseline2,
      data = epi, id = id,
      family = poisson("log")
    ),
    "rank-deficient model matrix: the column\\(s\\) 'lnbaseline2'"
  )
  expect_error(
    geewa(seizures ~ 0, data = epi, id = id, family = poisson("log")),
    "no columns"
  )
})

test_that("geewa rejects factor, character and multi-column responses for non-binomial families", {
  epi <- epilepsy
  epi$seizures_factor <- factor(epi$seizures)
  epi$seizures_character <- as.character(epi$seizures)
  expect_error(
    geewa(seizures_factor ~ treatment,
      data = epi, id = id,
      family = poisson("log")
    ),
    "unless the family is binomial-type"
  )
  expect_error(
    geewa(seizures_character ~ treatment,
      data = epi, id = id,
      family = gaussian()
    ),
    "unless the family is binomial-type"
  )
  expect_error(
    geewa(cbind(seizures, lnbaseline) ~ treatment,
      data = epi, id = id,
      family = poisson("log")
    ),
    "a matrix response must have one column"
  )
})

test_that("character repeated values are sorted alphabetically with a warning", {
  dat <- epilepsy
  visit_code <- c(1, 2, 10, 11)[as.integer(factor(dat$visit))]
  dat$visit_numeric <- visit_code
  dat$visit_character <- as.character(visit_code)
  dat$visit_factor <- factor(visit_code)
  fit_numeric <- geewa(
    seizures ~ treatment,
    family = poisson,
    data = dat,
    id = id,
    repeated = visit_numeric,
    corstr = "ar1"
  )
  expect_warning(
    fit_character <- geewa(
      seizures ~ treatment,
      family = poisson,
      data = dat,
      id = id,
      repeated = visit_character,
      corstr = "ar1"
    ),
    "'repeated' is a character vector",
    fixed = TRUE
  )
  expect_no_warning(
    fit_factor <- geewa(
      seizures ~ treatment,
      family = poisson,
      data = dat,
      id = id,
      repeated = visit_factor,
      corstr = "ar1"
    )
  )
  expect_equal(fit_factor$alpha, fit_numeric$alpha, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fit_character$alpha, fit_numeric$alpha)))
  expect_no_warning(runs_test(fit_character))
})

test_that("geewa rejects responses outside the support of the family", {
  dat <- epilepsy
  dat$shifted <- dat$seizures - 5
  expect_error(
    geewa(shifted ~ treatment, family = poisson, data = dat, id = id),
    "negative values not allowed for the 'Poisson' family",
    fixed = TRUE
  )
  expect_error(
    geewa(
      shifted ~ treatment,
      family = poisson,
      data = dat,
      id = id,
      beta_start = c(1, 0)
    ),
    "negative values not allowed for the 'Poisson' family",
    fixed = TRUE
  )
  expect_error(
    geewa(shifted ~ treatment, family = quasipoisson, data = dat, id = id),
    "negative values not allowed for the 'quasiPoisson' family",
    fixed = TRUE
  )
  expect_error(
    geewa(seizures ~ treatment, family = Gamma("log"), data = dat, id = id),
    "non-positive values not allowed for the 'Gamma' family",
    fixed = TRUE
  )
  expect_error(
    geewa(seizures ~ treatment, family = inverse.gaussian, data = dat, id = id),
    "positive values only are allowed for the 'inverse.gaussian' family",
    fixed = TRUE
  )
  expect_error(
    geewa(
      shifted ~ treatment,
      family = quasi(link = "log", variance = "mu"),
      data = dat,
      id = id
    ),
    "negative values not allowed for quasi(variance = \"mu\")",
    fixed = TRUE
  )
  expect_no_error(
    geewa(shifted ~ treatment, family = gaussian, data = dat, id = id)
  )
  expect_no_error(
    geewa(
      I(seizures + 1) ~ treatment,
      family = Gamma("log"),
      data = dat,
      id = id
    )
  )
})
