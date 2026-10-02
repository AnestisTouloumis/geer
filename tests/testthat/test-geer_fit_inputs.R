testthat::local_edition(3)


fit_inputs_data <- data.frame(
  id = rep(1:3, each = 2),
  t = rep(1:2, times = 3),
  y = c(1.2, 0.4, 2.2, 1.9, 0.7, 3.1),
  x = c(0.1, 0.5, 0.2, 0.9, 0.4, 0.8)
)


## ---------------------------------------------------------------- offset ----

test_that("extract_geer_offset returns zeros when no offset is supplied", {
  mf <- stats::model.frame(y ~ x, data = fit_inputs_data)

  expect_identical(geer:::extract_geer_offset(mf, 6L), rep(0, 6))
})


test_that("extract_geer_offset returns the supplied offset as a double vector", {
  offset <- seq(0.1, 0.6, length.out = 6)
  mf <- stats::model.frame(y ~ x, data = fit_inputs_data, offset = offset)

  out <- geer:::extract_geer_offset(mf, 6L)

  expect_type(out, "double")
  expect_equal(out, offset)
})


test_that("extract_geer_offset recycles a scalar offset", {
  ## A length-one offset cannot be stored in a model frame, so use the minimal
  ## list structure that stats::model.offset() reads.
  mf <- list("(offset)" = 2)

  expect_equal(geer:::extract_geer_offset(mf, 4L), rep(2, 4))
})


test_that("extract_geer_offset validates the offset", {
  expect_error(
    geer:::extract_geer_offset(list("(offset)" = c(1, 2, 3)), 6L),
    "response variable and 'offset' are not of same length"
  )

  mf_factor <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    offset = factor(rep(c("a", "b"), 3))
  )
  ## stats::model.offset() rejects a non-numeric offset before geer's own
  ## check is reached; either way the offset must not be accepted.
  expect_error(
    geer:::extract_geer_offset(mf_factor, 6L),
    "'offset' must be numeric"
  )

  mf_na <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    offset = c(0, NA, 0, 0, 0, 0),
    na.action = stats::na.pass
  )
  expect_error(
    geer:::extract_geer_offset(mf_na, 6L),
    "'offset' must be finite"
  )

  mf_inf <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    offset = c(0, Inf, 0, 0, 0, 0)
  )
  expect_error(
    geer:::extract_geer_offset(mf_inf, 6L),
    "'offset' must be finite"
  )
})


## --------------------------------------------------------- design matrix ----

test_that("build_geer_design_matrix returns the design matrix and its QR", {
  mf <- stats::model.frame(y ~ x, data = fit_inputs_data)

  out <- geer:::build_geer_design_matrix(mf)

  expect_identical(out$xnames, c("(Intercept)", "x"))
  expect_identical(dim(out$x), c(6L, 2L))
  expect_identical(out$qr$rank, 2L)
  expect_identical(out$assign, c(0L, 1L))
})


test_that("build_geer_design_matrix keeps a one-column design as a matrix", {
  mf <- stats::model.frame(y ~ 1, data = fit_inputs_data)

  out <- geer:::build_geer_design_matrix(mf)

  expect_true(is.matrix(out$x))
  expect_identical(dim(out$x), c(6L, 1L))
})


test_that("build_geer_design_matrix rejects a rank-deficient design", {
  data <- fit_inputs_data
  data$x2 <- 2 * data$x
  mf <- stats::model.frame(y ~ x + x2, data = data)

  expect_error(
    geer:::build_geer_design_matrix(mf),
    "rank-deficient model matrix"
  )
})


## -------------------------------------------------------------- control ----

test_that("normalize_geer_control fills in defaults and validates its input", {
  expect_identical(geer:::normalize_geer_control(NULL), geer_control())

  out <- geer:::normalize_geer_control(list(maxiter = 7, tolerance = 1e-4))
  expect_equal(out$maxiter, 7)
  expect_equal(out$tolerance, 1e-4)
  expect_identical(out$or_adding, geer_control()$or_adding)

  expect_error(
    geer:::normalize_geer_control("tight"),
    "'control' must be NULL or a list"
  )
  expect_error(
    geer:::normalize_geer_control(list(maxiter = -1)),
    "'maxiter' must be a positive integer"
  )
})


## --------------------------------------------------------------- family ----

test_that("normalize_family accepts names, functions and family objects", {
  from_name <- geer:::normalize_family("poisson")
  from_function <- geer:::normalize_family(poisson)
  from_object <- geer:::normalize_family(binomial(link = "probit"))

  expect_identical(from_name$family, "poisson")
  expect_identical(from_name$link, "log")
  expect_identical(from_function$family, "poisson")
  expect_identical(from_object$family, "binomial")
  expect_identical(from_object$link, "probit")
})


test_that("normalize_family rejects invalid specifications", {
  expect_error(
    geer:::normalize_family("no_such_family_function"),
    "'family' must be a family object, a family function, or the name of one"
  )
  expect_error(
    geer:::normalize_family(list(a = 1)),
    "'family' must be a valid family object"
  )
  expect_error(
    geer:::normalize_family(42),
    "'family' must be a valid family object"
  )
})


## -------------------------------------------------- response and weights ----

test_that("extract_geer_response_weights passes a numeric response through", {
  mf <- stats::model.frame(y ~ x, data = fit_inputs_data)

  out <- geer:::extract_geer_response_weights(mf, stats::gaussian())

  expect_identical(out$y, fit_inputs_data$y)
  expect_identical(out$weights, rep(1, 6))
})


test_that("extract_geer_response_weights reports a missing response", {
  mf <- stats::model.frame(~ x, data = fit_inputs_data)

  expect_error(
    geer:::extract_geer_response_weights(mf, stats::gaussian()),
    "response variable not found"
  )
})


test_that("extract_geer_response_weights rejects non-finite responses", {
  data <- fit_inputs_data
  data$y[3] <- NA_real_
  mf <- stats::model.frame(y ~ x, data = data, na.action = stats::na.pass)

  expect_error(
    geer:::extract_geer_response_weights(mf, stats::gaussian()),
    "response variable contains non-finite values"
  )
})


test_that("binomial responses may be coded as two-level factors or characters", {
  data <- fit_inputs_data
  data$yes_no <- rep(c("no", "yes", "yes"), 2)
  data$yes_no_factor <- factor(data$yes_no)

  out_chr <- geer:::extract_geer_response_weights(
    stats::model.frame(yes_no ~ x, data = data),
    stats::binomial()
  )
  out_fct <- geer:::extract_geer_response_weights(
    stats::model.frame(yes_no_factor ~ x, data = data),
    stats::binomial()
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
    geer:::extract_geer_response_weights(
      stats::model.frame(three_factor ~ x, data = data),
      stats::binomial()
    ),
    "a factor response must have exactly two levels"
  )
  expect_error(
    geer:::extract_geer_response_weights(
      stats::model.frame(three ~ x, data = data),
      stats::binomial()
    ),
    "exactly two (levels|distinct values)"
  )
})


test_that("binomial responses outside [0, 1] are rejected", {
  data <- fit_inputs_data
  data$bad <- c(0, 1, 2, 0, 1, 0)
  mf <- stats::model.frame(bad ~ x, data = data)

  expect_error(
    geer:::extract_geer_response_weights(mf, stats::binomial()),
    "the response must be coded as 0/1, proportions in \\[0, 1\\]"
  )
  expect_error(
    geer:::extract_geer_response_weights(
      mf,
      stats::quasi(link = "logit", variance = "mu(1-mu)")
    ),
    "the response must be coded as 0/1, proportions in \\[0, 1\\]"
  )
})


test_that("two-column binomial responses become proportions with trial weights", {
  data <- fit_inputs_data
  data$s <- c(2, 1, 3, 0, 4, 1)
  data$f <- c(1, 2, 0, 3, 1, 2)
  trials <- data$s + data$f

  out <- geer:::extract_geer_response_weights(
    stats::model.frame(cbind(s, f) ~ x, data = data),
    stats::binomial()
  )
  expect_equal(out$y, data$s / trials)
  expect_equal(out$weights, trials, ignore_attr = TRUE)

  prior <- c(1, 2, 1, 2, 1, 2)
  out_weighted <- geer:::extract_geer_response_weights(
    stats::model.frame(cbind(s, f) ~ x, data = data, weights = prior),
    stats::binomial()
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
    geer:::extract_geer_response_weights(
      stats::model.frame(cbind(s, f, g) ~ x, data = data),
      stats::binomial()
    ),
    "a matrix response must have exactly two columns"
  )

  data_na <- data
  data_na$s[2] <- NA_real_
  expect_error(
    geer:::extract_geer_response_weights(
      stats::model.frame(
        cbind(s, f) ~ x,
        data = data_na,
        na.action = stats::na.pass
      ),
      stats::binomial()
    ),
    "all entries must be finite"
  )

  data_negative <- data
  data_negative$s[2] <- -1
  expect_error(
    geer:::extract_geer_response_weights(
      stats::model.frame(cbind(s, f) ~ x, data = data_negative),
      stats::binomial()
    ),
    "all counts must be nonnegative"
  )

  data_zero <- data
  data_zero$s[2] <- 0
  data_zero$f[2] <- 0
  expect_error(
    geer:::extract_geer_response_weights(
      stats::model.frame(cbind(s, f) ~ x, data = data_zero),
      stats::binomial()
    ),
    "row sums \\(trials\\) must be positive"
  )
})


test_that("prior weights are validated", {
  mf_chr <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = rep(c("a", "b"), 3)
  )
  expect_error(
    geer:::extract_geer_response_weights(mf_chr, stats::gaussian()),
    "'weights' must be a numeric vector"
  )

  mf_na <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = c(1, NA, 1, 1, 1, 1),
    na.action = stats::na.pass
  )
  expect_error(
    geer:::extract_geer_response_weights(mf_na, stats::gaussian()),
    "'weights' must be finite"
  )

  mf_zero <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = c(1, 0, 1, 1, 1, 1)
  )
  expect_error(
    geer:::extract_geer_response_weights(mf_zero, stats::gaussian()),
    "'weights' must be strictly positive"
  )

  mf_ok <- stats::model.frame(
    y ~ x,
    data = fit_inputs_data,
    weights = c(1, 2, 3, 1, 2, 3)
  )
  expect_identical(
    geer:::extract_geer_response_weights(mf_ok, stats::gaussian())$weights,
    c(1, 2, 3, 1, 2, 3)
  )
})


## ----------------------------------------------------------- id/repeated ----

test_that("extract_geer_id_repeated numbers clusters and occasions", {
  mf <- stats::model.frame(y ~ x, data = fit_inputs_data, id = id)

  out <- geer:::extract_geer_id_repeated(mf, 6L)

  expect_equal(out$id, rep(1:3, each = 2))
  expect_equal(out$repeated, rep(1:2, times = 3))
})


test_that("extract_geer_id_repeated uses the supplied repeated variable", {
  data <- fit_inputs_data
  data$when <- rep(c("second", "first"), times = 3)
  mf <- stats::model.frame(y ~ x, data = data, id = id, repeated = when)

  out <- geer:::extract_geer_id_repeated(mf, 6L)

  ## Labels are mapped through factor(), so "first" -> 1 and "second" -> 2.
  expect_equal(out$repeated, rep(c(2, 1), times = 3))
})


test_that("extract_geer_id_repeated validates id and repeated", {
  mf_no_id <- stats::model.frame(y ~ x, data = fit_inputs_data)
  expect_error(
    geer:::extract_geer_id_repeated(mf_no_id, 6L),
    "'id' not found"
  )

  data_na_id <- fit_inputs_data
  data_na_id$id[2] <- NA
  mf_na_id <- stats::model.frame(
    y ~ x,
    data = data_na_id,
    id = id,
    na.action = stats::na.pass
  )
  expect_error(
    geer:::extract_geer_id_repeated(mf_na_id, 6L),
    "'id' cannot contain missing values"
  )

  mf <- stats::model.frame(y ~ x, data = fit_inputs_data, id = id)
  expect_error(
    geer:::extract_geer_id_repeated(mf, 5L),
    "response variable and 'id' are not of same length"
  )

  data_na_rep <- fit_inputs_data
  data_na_rep$t[3] <- NA
  mf_na_rep <- stats::model.frame(
    y ~ x,
    data = data_na_rep,
    id = id,
    repeated = t,
    na.action = stats::na.pass
  )
  expect_error(
    geer:::extract_geer_id_repeated(mf_na_rep, 6L),
    "'repeated' cannot contain missing values"
  )

  data_dup <- fit_inputs_data
  data_dup$t <- 1L
  mf_dup <- stats::model.frame(
    y ~ x,
    data = data_dup,
    id = id,
    repeated = t
  )
  expect_error(
    geer:::extract_geer_id_repeated(mf_dup, 6L),
    "'repeated' does not have unique values per 'id'"
  )
})


## ------------------------------------------------------------ phi, use_p ----

test_that("normalize_phi resolves fixed and estimated dispersion", {
  expect_identical(
    geer:::normalize_phi(FALSE, 3),
    list(phi_fixed = FALSE, phi_value = 1)
  )
  expect_identical(
    geer:::normalize_phi(TRUE, 2),
    list(phi_fixed = TRUE, phi_value = 2)
  )
})


test_that("normalize_phi validates its arguments", {
  for (bad in list(NA, "yes", c(TRUE, FALSE), NULL, 1)) {
    expect_error(
      geer:::normalize_phi(bad, 1),
      "'phi_fixed' must be a single non-missing logical value"
    )
  }
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), numeric(0))) {
    expect_error(
      geer:::normalize_phi(TRUE, bad),
      "'phi_value' must be a single positive number"
    )
  }
})


test_that("normalize_use_p returns a validated logical flag", {
  expect_true(geer:::normalize_use_p(TRUE))
  expect_false(geer:::normalize_use_p(FALSE))
  for (bad in list(NA, "TRUE", c(TRUE, FALSE), NULL, 1)) {
    expect_error(
      geer:::normalize_use_p(bad),
      "'use_p' must be a single non-missing logical value"
    )
  }
})
