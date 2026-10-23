make_mcar_dropout_data <- function() {
  withr::local_seed(7401)
  subjects <- 140L
  occasions <- 5L
  id <- rep(seq_len(subjects), each = occasions)
  visit <- rep(seq_len(occasions), times = subjects)
  treatment_subject <- rep(c(0, 1), length.out = subjects)
  treatment <- rep(treatment_subject, each = occasions)
  age_subject <- seq(35, 70, length.out = subjects)
  age <- rep(age_subject, each = occasions)
  y_complete <-
    2.5 + 0.5 * treatment + 0.025 * age + 0.18 * visit +
    rnorm(length(id), sd = 0.8)
  y <- y_complete

  for (i in seq_len(subjects)) {
    idx <- which(id == i)
    for (j in 2:occasions) {
      eta <-
        -3.2 + 0.55 * treatment[idx[j]] +
        0.45 * y_complete[idx[j - 1L]] + 0.12 * (j - 2L)
      if (rbinom(1L, 1L, plogis(eta)) == 1L) {
        y[idx[j:occasions]] <- NA_real_
        break
      }
    }
  }

  data.frame(id, visit, treatment, age, y, y_complete)
}

make_direct_dropout_risk_set <- function(dat) {
  rows <- vector("list", 0L)
  counter <- 0L
  for (i in unique(dat$id)) {
    idx <- which(dat$id == i)
    idx <- idx[order(dat$visit[idx])]
    if (length(idx) < 2L) next
    for (j in 2:length(idx)) {
      if (!is.na(dat$y[idx[j - 1L]])) {
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          missing = as.integer(is.na(dat$y[idx[j]])),
          id = dat$id[idx[j]],
          visit = dat$visit[idx[j]],
          treatment = dat$treatment[idx[j]],
          age = dat$age[idx[j]],
          lag_response = dat$y[idx[j - 1L]]
        )
      }
    }
  }
  do.call(rbind, rows)
}

test_that("mcar_logistic_test agrees with a direct geewa_binary risk-set fit", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  out <- mcar_logistic_test(
    fit,
    formula = ~ treatment + age + visit,
    orstr = "independence"
  )

  risk <- make_direct_dropout_risk_set(dat)
  direct <- geewa_binary(
    missing ~ factor(visit) + treatment + age + lag_response,
    data = risk,
    id = id,
    repeated = visit,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )
  beta <- coef(direct)
  covariance <- vcov(direct, cov_type = "bias-corrected")

  overall_names <- c("treatment", "age", "lag_response")
  overall_index <- match(overall_names, names(beta))
  expected_overall <- as.numeric(crossprod(
    beta[overall_index],
    solve(
      covariance[overall_index, overall_index, drop = FALSE],
      beta[overall_index]
    )
  ))
  lag_index <- match("lag_response", names(beta))
  expected_history <- as.numeric(beta[lag_index]^2 / covariance[lag_index, lag_index])

  expect_s3_class(out, "htest")
  expect_s3_class(out$model, "geer")
  expect_match(deparse(out$model$call[[1L]]), "geewa_binary")
  expect_equal(unname(out$statistic), expected_history, tolerance = 1e-8)
  expect_equal(
    out$tests$statistic[out$tests$test == "response_history"],
    expected_history,
    tolerance = 1e-8
  )
  expect_equal(
    out$tests$statistic[out$tests$test == "overall"],
    expected_overall,
    tolerance = 1e-8
  )
  expect_equal(unname(out$parameter), 1)
  expect_identical(out$orstr, "independence")
  expect_identical(out$test, "wald")
  expect_identical(out$cov_type, "bias-corrected")
  expect_null(out$pmethod)
  expect_true(all(out$tests$procedure == "wald"))
  expect_false(out$intermittent)
  expect_true("visit" %in% out$dropped_covariates)
  expect_identical(
    out$coefficients$term[grepl("^(treatment|age|lag_response)$", out$coefficients$term)],
    c("treatment", "age", "lag_response")
  )
})

test_that("mcar_logistic_test separates response-history and covariate tests", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  out <- mcar_logistic_test(fit, formula = ~ treatment + age)

  expect_identical(
    out$tests$test,
    c("response_history", "covariates", "overall")
  )
  expect_equal(out$tests$df, c(1, 2, 3))
  expect_true(all(is.finite(out$tests$statistic)))
  expect_true(all(out$tests$p.value >= 0 & out$tests$p.value <= 1))
  expect_equal(
    out$p.value,
    out$tests$p.value[out$tests$test == "response_history"],
    tolerance = 1e-12
  )
})

test_that("mcar_logistic_test uses fitted predictors and absorbs occasion main effects", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  default <- mcar_logistic_test(fit)
  explicit <- mcar_logistic_test(fit, formula = ~ treatment + age + visit)

  expect_equal(default$statistic, explicit$statistic, tolerance = 1e-8)
  expect_equal(default$tests$p.value, explicit$tests$p.value, tolerance = 1e-10)
  expect_true("visit" %in% default$dropped_covariates)
})

test_that("mcar_logistic_test supports every geer hypothesis-test procedure", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  risk <- make_direct_dropout_risk_set(dat)
  full <- geewa_binary(
    missing ~ factor(visit) + treatment + age + lag_response,
    data = risk,
    id = id,
    repeated = visit,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )
  null_history <- geewa_binary(
    missing ~ factor(visit) + treatment + age,
    data = risk,
    id = id,
    repeated = visit,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )

  procedures <- c(
    "wald", "score", "working-wald", "working-score", "working-lrt"
  )
  for (procedure in procedures) {
    out <- mcar_logistic_test(
      fit,
      formula = ~ treatment + age,
      orstr = "independence",
      test = procedure,
      cov_type = "robust",
      pmethod = "rao-scott"
    )
    expected <- anova(
      null_history,
      full,
      test = procedure,
      cov_type = "robust",
      pmethod = "rao-scott"
    )

    expect_identical(out$test, procedure)
    expect_true(all(out$tests$procedure == procedure))
    expect_equal(unname(out$statistic), expected$Chi[2L], tolerance = 1e-8)
    expect_equal(unname(out$parameter), expected$Df[2L], tolerance = 1e-8)
    expect_equal(out$p.value, expected$`Pr(>Chi)`[2L], tolerance = 1e-8)
  }

  out_satt <- mcar_logistic_test(
    fit,
    formula = ~ treatment + age,
    test = "working-score",
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expected_satt <- anova(
    null_history,
    full,
    test = "working-score",
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expect_identical(out_satt$pmethod, "satterthwaite")
  expect_equal(
    unname(out_satt$statistic),
    expected_satt$Chi[2L],
    tolerance = 1e-8
  )
  expect_equal(
    unname(out_satt$parameter),
    expected_satt$Df[2L],
    tolerance = 1e-8
  )
})

test_that("mcar_logistic_test supports transformed and factor covariates", {
  dat <- make_mcar_dropout_data()
  dat$group <- factor(ifelse(dat$age > 52, "older", "younger"))
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  out <- mcar_logistic_test(
    fit,
    formula = ~ treatment + log(age) + group
  )

  expect_s3_class(out$model, "geer")
  expect_true("log(age)" %in% out$coefficients$term)
  expect_true(any(grepl("^group", out$coefficients$term)))
  expect_true("lag_response" %in% out$coefficients$term)
})

test_that("mcar_logistic_test permits a response-history-only diagnostic", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  out <- mcar_logistic_test(fit, formula = ~1)

  expect_equal(out$tests$df, c(1, 0, 1))
  expect_equal(
    out$tests$statistic[out$tests$test == "response_history"],
    out$tests$statistic[out$tests$test == "overall"],
    tolerance = 1e-12
  )
})

test_that("mcar_logistic_test requires fully observed covariates", {
  dat <- make_mcar_dropout_data()
  dat$aux <- seq_len(nrow(dat))
  dat$aux[c(5, 18)] <- NA_real_
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  expect_error(
    mcar_logistic_test(fit, formula = ~ treatment + aux),
    "fully observed"
  )
})

test_that("mcar_logistic_test can recover data supplied explicitly", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )
  fit$data <- NULL

  out <- mcar_logistic_test(fit, data = dat, formula = ~ treatment + age)
  expect_s3_class(out, "htest")
})

test_that("mcar_logistic_test passes control to the binary GEE fits", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )
  expect_error(
    suppressWarnings(mcar_logistic_test(fit,
      formula = ~ treatment + age,
      orstr = "exchangeable", control = geer_control(maxiter = 1, tolerance = 1e-15)
    )),
    paste0(
      "mcar_logistic_test failed: the binary GEE fit of the missingness ",
      "model did not converge"
    ),
    fixed = TRUE
  )
})

test_that("mcar_logistic_test warns for intermittent response missingness", {
  dat <- make_mcar_dropout_data()
  subject <- dat$id == 1
  dat$y[subject] <- dat$y_complete[subject]
  subject_rows <- which(subject)
  dat$y[subject_rows[3L]] <- NA_real_
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  expect_warning(
    out <- mcar_logistic_test(fit, formula = ~ treatment + age),
    "intermittent missingness detected"
  )
  expect_true(out$intermittent)
})

test_that("mcar_logistic_test handles no response missingness", {
  id <- rep(seq_len(40), each = 4)
  visit <- rep(seq_len(4), times = 40)
  x <- rep(c(0, 1), length.out = length(id))
  y <- 1 + x + 0.2 * visit + sin(id / 4)
  dat <- data.frame(id, visit, x, y)
  fit <- geewa(
    y ~ x + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  out <- mcar_logistic_test(fit)

  expect_equal(unname(out$statistic), 0)
  expect_equal(unname(out$parameter), 0)
  expect_equal(out$p.value, 1)
  expect_identical(out$missing, 0L)
  expect_identical(out$observed, 120L)
  expect_identical(out$transitions, 120L)
  expect_null(out$model)
  expect_null(out$coefficients)
})

test_that("mcar_logistic_test validates its covariate formula and odds-ratio structure", {
  dat <- make_mcar_dropout_data()
  fit <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  expect_error(mcar_logistic_test(fit, formula = treatment ~ age), "one-sided")
  expect_error(mcar_logistic_test(fit, formula = ~ 0 + treatment), "include an intercept")
  expect_error(mcar_logistic_test(fit, formula = "~ treatment"), "one-sided formula")
  expect_error(
    mcar_logistic_test(
      fit,
      formula = ~ treatment + age,
      orstr = "exchangeable",
      test = "working-lrt"
    ),
    "requires an independence working structure"
  )
})

test_that("normalize_mcar_response codes two-level binary responses as 0/1", {
  binomial_fit <- list(family = binomial())
  response <- factor(c("no", "yes", NA, "yes"), levels = c("no", "yes"))
  expect_identical(normalize_mcar_response(response, binomial_fit), c(0, 1, NA, 1))
  expect_identical(
    normalize_mcar_response(c("a", "b", NA, "b"), binomial_fit),
    c(0, 1, NA, 1)
  )
  quasi_fit <- list(family = quasi(link = "logit", variance = "mu(1-mu)"))
  expect_identical(normalize_mcar_response(response, quasi_fit), c(0, 1, NA, 1))
  expect_error(
    normalize_mcar_response(factor(c("a", "b", "c")), binomial_fit),
    "factor response must have exactly two levels"
  )
  expect_error(
    normalize_mcar_response(c("a", "b", "c"), binomial_fit),
    "character response must have exactly two distinct values"
  )
  expect_error(
    normalize_mcar_response(factor(c("a", "b")), list(family = gaussian())),
    "must be numeric or binary with two levels"
  )
  expect_error(
    normalize_mcar_response(c(1, Inf), list(family = gaussian())),
    "observed response values must be finite"
  )
  expect_error(
    normalize_mcar_response(cbind(1:2, 2:1), binomial_fit),
    "require a univariate response"
  )
})

test_that("mcar_logistic_test gives the same result for 0/1 and factor responses", {
  dat <- make_mcar_dropout_data()
  dat$high <- as.numeric(dat$y > median(dat$y, na.rm = TRUE))
  dat$high_factor <- factor(dat$high, levels = c(0, 1), labels = c("low", "high"))
  fit_numeric <- geewa(
    high ~ treatment + age,
    data = dat,
    id = id,
    repeated = visit,
    family = binomial(),
    corstr = "independence"
  )
  fit_factor <- update(fit_numeric, high_factor ~ .)
  out_numeric <- mcar_logistic_test(fit_numeric, formula = ~ treatment + age)
  out_factor <- mcar_logistic_test(fit_factor, formula = ~ treatment + age)
  expect_equal(out_factor$statistic, out_numeric$statistic, tolerance = 1e-10)
  expect_equal(out_factor$parameter, out_numeric$parameter)
  expect_equal(out_factor$p.value, out_numeric$p.value, tolerance = 1e-10)
})

test_that("build_mcar_null_formula falls back to an intercept-only model", {
  expect_identical(
    deparse(build_mcar_null_formula(character(0))),
    "mcar_missing ~ 1"
  )
  expect_identical(
    deparse(build_mcar_null_formula(c("a", "b"))),
    "mcar_missing ~ a + b"
  )
})

test_that("mcar_logistic_test works on data with only two occasions", {
  dat <- make_mcar_dropout_data()
  dat2 <- dat[dat$visit <= 2L, , drop = FALSE]
  fit2 <- geewa(
    y ~ treatment + age,
    data = dat2,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )
  out2 <- mcar_logistic_test(
    fit2,
    formula = ~ treatment + age,
    orstr = "independence"
  )
  expect_s3_class(out2, "htest")
  expect_true(is.finite(out2$statistic))
  expect_gte(out2$p.value, 0)
  expect_lte(out2$p.value, 1)
})

test_that("mcar_logistic_test warns when dropout is recorded by absent rows", {
  dat <- make_mcar_dropout_data()
  observed <- dat[!is.na(dat$y), , drop = FALSE]
  fit <- geewa(
    y ~ treatment + age + visit,
    data = observed,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )
  expect_warning(
    out <- mcar_logistic_test(fit),
    "dropout recorded by omitted rows cannot be detected",
    fixed = TRUE
  )
  expect_identical(out$missing, 0L)
  expect_equal(out$p.value, 1)

  fit_full <- geewa(
    y ~ treatment + age + visit,
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )
  expect_no_warning(out_full <- mcar_logistic_test(fit_full))
  expect_gt(out_full$missing, 0L)
})

test_that("occasions without dropouts are omitted from the risk set", {
  dat <- make_mcar_dropout_data()
  observed_at_4 <- dat$id[dat$visit == 4L & !is.na(dat$y)]
  last <- dat$visit == 5L & dat$id %in% observed_at_4
  dat$y[last] <- dat$y_complete[last]
  fit <- geewa(
    formula = y ~ treatment + age,
    family = gaussian(),
    data = dat,
    id = id,
    repeated = visit
  )
  expect_warning(
    out <- mcar_logistic_test(fit, formula = ~ treatment + age),
    paste0(
      "occasion(s) 5 were omitted from the risk set because they have no ",
      "dropouts or no observed responses"
    ),
    fixed = TRUE
  )
  risk <- make_direct_dropout_risk_set(dat)
  expect_warning(
    reference <- stats::glm(
      missing ~ factor(visit) + treatment + age + lag_response,
      family = stats::binomial(),
      data = risk,
      control = stats::glm.control(epsilon = 1e-14, maxit = 100L)
    ),
    "fitted probabilities numerically 0 or 1 occurred",
    fixed = TRUE
  )
  terms <- c("treatment", "age", "lag_response")
  estimates <- out$coefficients$estimate[match(terms, out$coefficients$term)]
  expect_equal(estimates, unname(coef(reference)[terms]), tolerance = 1e-6)
  expect_identical(out$transitions, sum(risk$visit != 5L))
})

test_that("a risk set without an estimable occasion is an error", {
  dat <- data.frame(
    id = rep(seq_len(10L), each = 3L),
    visit = rep(seq_len(3L), times = 10L),
    x = rep(c(0, 1), each = 15L)
  )
  dat$y <- dat$visit + dat$x + rep(seq(-1, 1, length.out = 10L), each = 3L)
  dat$y[dat$visit == 3L] <- NA_real_
  fit <- geewa(
    formula = y ~ x,
    family = gaussian(),
    data = dat,
    id = id,
    repeated = visit
  )
  expect_error(
    mcar_logistic_test(fit, formula = ~ x),
    paste0(
      "every occasion has either no dropouts or no observed responses in ",
      "the risk set"
    ),
    fixed = TRUE
  )
})

