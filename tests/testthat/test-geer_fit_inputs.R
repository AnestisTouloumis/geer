fit_inputs_data <- data.frame(
  id = rep(1:3, each = 2),
  t = rep(1:2, times = 3),
  y = c(1.2, 0.4, 2.2, 1.9, 0.7, 3.1),
  x = c(0.1, 0.5, 0.2, 0.9, 0.4, 0.8)
)

fit_or_offset_formula <- geewa_binary(
  formula = ecg ~ treatment + offset(period),
  link = "logit",
  id = id,
  data = test_data$cerebrovascular,
  orstr = "independence",
  method = "gee"
)

fit_or_offset_argument <- geewa_binary(
  formula = ecg ~ treatment,
  link = "logit",
  id = id,
  offset = period,
  data = test_data$cerebrovascular,
  orstr = "independence",
  method = "gee"
)

fit_cc_offset_formula <- geewa(
  formula = ecg ~ treatment + offset(period),
  family = binomial(link = "logit"),
  id = id,
  data = test_data$cerebrovascular,
  corstr = "independence",
  method = "gee",
  phi_fixed = TRUE,
  phi_value = 1
)

test_that("extract_geer_offset returns zeros when no offset is supplied", {
  mf <- model.frame(y ~ x, data = fit_inputs_data)

  expect_identical(extract_geer_offset(mf, 6L), rep(0, 6))
})

test_that("extract_geer_offset returns the supplied offset as a double vector", {
  offset <- seq(0.1, 0.6, length.out = 6)
  mf <- model.frame(y ~ x, data = fit_inputs_data, offset = offset)

  out <- extract_geer_offset(mf, 6L)

  expect_type(out, "double")
  expect_equal(out, offset)
})

test_that("extract_geer_offset recycles a scalar offset", {
  mf <- list("(offset)" = 2)

  expect_equal(extract_geer_offset(mf, 4L), rep(2, 4))
})

test_that("extract_geer_offset validates the offset", {
  expect_error(
    extract_geer_offset(list("(offset)" = c(1, 2, 3)), 6L),
    "'offset' must have one value per observation"
  )

  mf_factor <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    offset = factor(rep(c("a", "b"), 3))
  )
  expect_error(
    extract_geer_offset(mf_factor, 6L),
    "'offset' must be numeric"
  )

  mf_na <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    offset = c(0, NA, 0, 0, 0, 0),
    na.action = na.pass
  )
  expect_error(
    extract_geer_offset(mf_na, 6L),
    "'offset' must be finite"
  )

  mf_inf <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    offset = c(0, Inf, 0, 0, 0, 0)
  )
  expect_error(
    extract_geer_offset(mf_inf, 6L),
    "'offset' must be finite"
  )
})

test_that("build_geer_design_matrix returns the design matrix and its QR", {
  mf <- model.frame(y ~ x, data = fit_inputs_data)

  out <- build_geer_design_matrix(mf)

  expect_identical(out$xnames, c("(Intercept)", "x"))
  expect_identical(dim(out$x), c(6L, 2L))
  expect_identical(out$qr$rank, 2L)
  expect_identical(out$assign, c(0L, 1L))
})

test_that("build_geer_design_matrix keeps a one-column design as a matrix", {
  mf <- model.frame(y ~ 1, data = fit_inputs_data)

  out <- build_geer_design_matrix(mf)

  expect_true(is.matrix(out$x))
  expect_identical(dim(out$x), c(6L, 1L))
})

test_that("a one-column design keeps its column name and assign attribute", {
  fit_no_intercept <- geewa(
    seizures ~ 0 + lnbaseline,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence"
  )
  expect_identical(colnames(fit_no_intercept$x), "lnbaseline")
  expect_identical(attr(fit_no_intercept$x, "assign"), 1L)

  fit_intercept_only <- geewa(
    seizures ~ 1,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence"
  )
  expect_identical(colnames(fit_intercept_only$x), "(Intercept)")
  expect_identical(attr(fit_intercept_only$x, "assign"), 0L)
})

test_that("build_geer_design_matrix rejects a rank-deficient design", {
  data <- fit_inputs_data
  data$x2 <- 2 * data$x
  mf <- model.frame(y ~ x + x2, data = data)

  expect_error(
    build_geer_design_matrix(mf),
    "rank-deficient model matrix"
  )
})

test_that("normalize_family accepts names, functions and family objects", {
  from_name <- normalize_family("poisson")
  from_function <- normalize_family(poisson)
  from_object <- normalize_family(binomial(link = "probit"))

  expect_identical(from_name$family, "poisson")
  expect_identical(from_name$link, "log")
  expect_identical(from_function$family, "poisson")
  expect_identical(from_object$family, "binomial")
  expect_identical(from_object$link, "probit")
})

test_that("normalize_family rejects invalid specifications", {
  expect_error(
    normalize_family("no_such_family_function"),
    "'family' must be a family object, a family function, or the name of one"
  )
  expect_error(
    normalize_family(list(a = 1)),
    "'family' must be a valid family object"
  )
  expect_error(
    normalize_family(42),
    "'family' must be a valid family object"
  )
})

test_that("extract_geer_response_weights passes a numeric response through", {
  mf <- model.frame(y ~ x, data = fit_inputs_data)

  out <- extract_geer_response_weights(mf, gaussian())

  expect_identical(out$y, fit_inputs_data$y)
  expect_identical(out$weights, rep(1, 6))
})

test_that("extract_geer_response_weights reports a missing response", {
  mf <- model.frame(~x, data = fit_inputs_data)

  expect_error(
    extract_geer_response_weights(mf, gaussian()),
    "the formula must have a response"
  )
})

test_that("extract_geer_response_weights rejects non-finite responses", {
  data <- fit_inputs_data
  data$y[3] <- NA_real_
  mf <- model.frame(y ~ x, data = data, na.action = na.pass)

  expect_error(
    extract_geer_response_weights(mf, gaussian()),
    "the response must contain only finite values"
  )
})

test_that("binomial responses may be coded as two-level factors or characters", {
  data <- fit_inputs_data
  data$yes_no <- rep(c("no", "yes", "yes"), 2)
  data$yes_no_factor <- factor(data$yes_no)

  out_chr <- extract_geer_response_weights(
    model.frame(yes_no ~ x, data = data),
    binomial()
  )
  out_fct <- extract_geer_response_weights(
    model.frame(yes_no_factor ~ x, data = data),
    binomial()
  )

  expected <- as.numeric(data$yes_no == "yes")
  expect_identical(out_chr$y, expected)
  expect_identical(out_fct$y, expected)
})

test_that("binomial responses with more than two categories are rejected", {
  data <- fit_inputs_data
  data$three <- rep(c("a", "b", "c"), 2)
  data$three_factor <- factor(data$three)

  expect_error(
    extract_geer_response_weights(
      model.frame(three_factor ~ x, data = data),
      binomial()
    ),
    "a factor response must have exactly two levels"
  )
  expect_error(
    extract_geer_response_weights(
      model.frame(three ~ x, data = data),
      binomial()
    ),
    "exactly two (levels|distinct values)"
  )
})

test_that("binomial responses outside [0, 1] are rejected", {
  data <- fit_inputs_data
  data$bad <- c(0, 1, 2, 0, 1, 0)
  mf <- model.frame(bad ~ x, data = data)

  expect_error(
    extract_geer_response_weights(mf, binomial()),
    "the response must be coded as 0/1, proportions in \\[0, 1\\]"
  )
  expect_error(
    extract_geer_response_weights(
      mf,
      quasi(link = "logit", variance = "mu(1-mu)")
    ),
    "the response must be coded as 0/1, proportions in \\[0, 1\\]"
  )
})

test_that("two-column binomial responses become proportions with trial weights", {
  data <- fit_inputs_data
  data$s <- c(2, 1, 3, 0, 4, 1)
  data$f <- c(1, 2, 0, 3, 1, 2)
  trials <- data$s + data$f

  out <- extract_geer_response_weights(
    model.frame(cbind(s, f) ~ x, data = data),
    binomial()
  )
  expect_equal(out$y, data$s / trials)
  expect_equal(out$weights, trials, ignore_attr = TRUE)

  prior <- c(1, 2, 1, 2, 1, 2)
  out_weighted <- extract_geer_response_weights(
    model.frame(cbind(s, f) ~ x, data = data, weights = prior),
    binomial()
  )
  expect_equal(out_weighted$y, data$s / trials)
  expect_equal(out_weighted$weights, prior * trials, ignore_attr = TRUE)
})

test_that("binomial matrix responses are validated", {
  data <- fit_inputs_data
  data$s <- c(2, 1, 3, 0, 4, 1)
  data$f <- c(1, 2, 0, 3, 1, 2)
  data$g <- c(1, 1, 1, 1, 1, 1)

  expect_error(
    extract_geer_response_weights(
      model.frame(cbind(s, f, g) ~ x, data = data),
      binomial()
    ),
    "a matrix response must have exactly two columns"
  )

  data_na <- data
  data_na$s[2] <- NA_real_
  expect_error(
    extract_geer_response_weights(
      model.frame(
        cbind(s, f) ~ x,
        data = data_na,
        na.action = na.pass
      ),
      binomial()
    ),
    "all entries must be finite"
  )

  data_negative <- data
  data_negative$s[2] <- -1
  expect_error(
    extract_geer_response_weights(
      model.frame(cbind(s, f) ~ x, data = data_negative),
      binomial()
    ),
    "all counts must be nonnegative"
  )

  data_zero <- data
  data_zero$s[2] <- 0
  data_zero$f[2] <- 0
  expect_error(
    extract_geer_response_weights(
      model.frame(cbind(s, f) ~ x, data = data_zero),
      binomial()
    ),
    "row sums \\(trials\\) must be positive"
  )
})

test_that("prior weights are validated", {
  mf_chr <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = rep(c("a", "b"), 3)
  )
  expect_error(
    extract_geer_response_weights(mf_chr, gaussian()),
    "'weights' must be a numeric vector"
  )

  mf_na <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = c(1, NA, 1, 1, 1, 1),
    na.action = na.pass
  )
  expect_error(
    extract_geer_response_weights(mf_na, gaussian()),
    "'weights' must be finite"
  )

  mf_zero <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = c(1, 0, 1, 1, 1, 1)
  )
  expect_error(
    extract_geer_response_weights(mf_zero, gaussian()),
    "'weights' must be strictly positive"
  )

  mf_ok <- model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = c(1, 2, 3, 1, 2, 3)
  )
  expect_identical(
    extract_geer_response_weights(mf_ok, gaussian())$weights,
    c(1, 2, 3, 1, 2, 3)
  )
})

test_that("extract_geer_id_repeated numbers clusters and occasions", {
  mf <- model.frame(y ~ x, data = fit_inputs_data, id = id)

  out <- extract_geer_id_repeated(mf, 6L)

  expect_equal(out$id, rep(1:3, each = 2))
  expect_equal(out$repeated, rep(1:2, times = 3))
})

test_that("extract_geer_id_repeated uses the supplied repeated variable", {
  data <- fit_inputs_data
  data$when <- rep(c("second", "first"), times = 3)
  mf <- model.frame(y ~ x, data = data, id = id, repeated = when)

  expect_error(
    extract_geer_id_repeated(mf, 6L),
    paste0(
      "'repeated' must be numeric or a factor; a character vector has no ",
      "natural order"
    ),
    fixed = TRUE
  )
  data$when <- factor(data$when, levels = c("first", "second"))
  mf <- model.frame(y ~ x, data = data, id = id, repeated = when)
  out <- extract_geer_id_repeated(mf, 6L)
  expect_equal(out$repeated, rep(c(2, 1), times = 3))
})

test_that("extract_geer_id_repeated validates id and repeated", {
  mf_no_id <- model.frame(y ~ x, data = fit_inputs_data)
  expect_error(
    extract_geer_id_repeated(mf_no_id, 6L),
    "'id' must be supplied"
  )

  data_na_id <- fit_inputs_data
  data_na_id$id[2] <- NA
  mf_na_id <- model.frame(
    y ~ x,
    data = data_na_id,
    id = id,
    na.action = na.pass
  )
  expect_error(
    extract_geer_id_repeated(mf_na_id, 6L),
    "'id' must not contain missing values"
  )

  mf <- model.frame(y ~ x, data = fit_inputs_data, id = id)
  expect_error(
    extract_geer_id_repeated(mf, 5L),
    "'id' must have one value per observation"
  )

  data_na_rep <- fit_inputs_data
  data_na_rep$t[3] <- NA
  mf_na_rep <- model.frame(
    y ~ x,
    data = data_na_rep,
    id = id,
    repeated = t,
    na.action = na.pass
  )
  expect_error(
    extract_geer_id_repeated(mf_na_rep, 6L),
    "'repeated' must not contain missing values"
  )

  data_dup <- fit_inputs_data
  data_dup$t <- 1L
  mf_dup <- model.frame(
    y ~ x,
    data = data_dup,
    id = id,
    repeated = t
  )
  expect_error(
    extract_geer_id_repeated(mf_dup, 6L),
    "'repeated' must have distinct values within each cluster"
  )
})

test_that("normalize_phi resolves fixed and estimated dispersion", {
  expect_identical(
    normalize_phi(FALSE, 3),
    list(phi_fixed = FALSE, phi_value = 1)
  )
  expect_identical(
    normalize_phi(TRUE, 2),
    list(phi_fixed = TRUE, phi_value = 2)
  )
})

test_that("normalize_phi validates its arguments", {
  for (bad in list(NA, "yes", c(TRUE, FALSE), NULL, 1)) {
    expect_error(
      normalize_phi(bad, 1),
      "'phi_fixed' must be TRUE or FALSE"
    )
  }
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), numeric(0))) {
    expect_error(
      normalize_phi(TRUE, bad),
      "'phi_value' must be a positive finite number"
    )
  }
})

test_that("normalize_use_p returns a validated logical flag", {
  expect_true(normalize_use_p(TRUE))
  expect_false(normalize_use_p(FALSE))
  for (bad in list(NA, "TRUE", c(TRUE, FALSE), NULL, 1)) {
    expect_error(
      normalize_use_p(bad),
      "'use_p' must be TRUE or FALSE"
    )
  }
})

test_that("geewa_binary treats formula and argument offsets equivalently", {
  expect_s3_class(fit_or_offset_formula, "geer")
  expect_s3_class(fit_or_offset_argument, "geer")
  expect_equal(coef(fit_or_offset_formula), coef(fit_or_offset_argument))
})

test_that("geewa and geewa_binary agree for the same offset specification under independence", {
  expect_s3_class(fit_cc_offset_formula, "geer")
  expect_equal(coef(fit_or_offset_formula), coef(fit_cc_offset_formula))
})

test_that("predict uses offset supplied via the offset argument with newdata", {
  nd <- test_data$cerebrovascular[1:8, , drop = FALSE]
  pred_formula <- predict(
    fit_or_offset_formula,
    newdata = nd,
    type = "link"
  )
  pred_argument <- predict(
    fit_or_offset_argument,
    newdata = nd,
    type = "link"
  )
  expect_equal(pred_formula, pred_argument)
})

test_that("formula offsets are reported by extract_formula_offset_terms", {
  expect_identical(
    extract_formula_offset_terms(fit_or_offset_formula),
    "offset(period)"
  )
  expect_identical(
    extract_formula_offset_terms(fit_or_offset_argument),
    character(0)
  )
})

test_that("sequential anova keeps a formula offset in the null model", {
  fit_null <- geewa(
    formula = ecg ~ 1 + offset(period),
    family = binomial(link = "logit"),
    id = id,
    data = test_data$cerebrovascular,
    corstr = "independence",
    method = "gee",
    phi_fixed = TRUE,
    phi_value = 1
  )
  tab <- anova(fit_cc_offset_formula, test = "score", cov_type = "robust")
  expected <- compute_score_test(fit_null, fit_cc_offset_formula, "robust")
  expect_equal(tab["treatment", "Chi"], expected$test_stat)
})

test_that("a scalar offset is recycled and equals the equivalent vector offset", {
  fit_scalar <- geewa(
    ecg ~ treatment,
    family = binomial(link = "logit"),
    id = id,
    data = test_data$cerebrovascular,
    offset = 0.5
  )
  fit_vector <- geewa(
    ecg ~ treatment,
    family = binomial(link = "logit"),
    id = id,
    data = test_data$cerebrovascular,
    offset = rep(0.5, nrow(test_data$cerebrovascular))
  )
  expect_equal(coef(fit_scalar), coef(fit_vector), tolerance = 1e-8)
  expect_true(all(fit_scalar$call$offset == 0.5))
  fit_subset <- geewa(
    ecg ~ treatment,
    family = binomial(link = "logit"),
    id = id,
    data = test_data$cerebrovascular,
    offset = 0.5,
    subset = as.numeric(id) <= 40
  )
  expect_s3_class(fit_subset, "geer")
})

test_that("refits for anova, add1, drop1 and step_p keep a formula offset", {
  refit_drop <- refit_geer(
    fit_cc_offset_formula,
    . ~ . - treatment,
    "anova.geer"
  )
  expect_identical(extract_formula_offset_terms(refit_drop), "offset(period)")
  refit_add <- refit_geer(
    fit_cc_offset_formula,
    . ~ . + factor(period),
    "anova.geer"
  )
  expect_identical(extract_formula_offset_terms(refit_add), "offset(period)")
  out_drop <- drop1(fit_cc_offset_formula, test = "wald", cov_type = "robust")
  expect_s3_class(out_drop, "anova")
  expect_true(all(is.finite(out_drop[["CIC"]])))
})
