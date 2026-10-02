testthat::local_edition(3)


make_mcar_dropout_data <- function() {
  set.seed(7401)
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
    stats::rnorm(length(id), sd = 0.8)
  y <- y_complete

  for (i in seq_len(subjects)) {
    idx <- which(id == i)
    for (j in 2:occasions) {
      eta <-
        -3.2 + 0.55 * treatment[idx[j]] +
        0.45 * y_complete[idx[j - 1L]] + 0.12 * (j - 2L)
      if (stats::rbinom(1L, 1L, stats::plogis(eta)) == 1L) {
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

  out <- mcar_logistic_test(fit, formula = ~ 1)

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

  expect_error(mcar_logistic_test(fit), "original data are not available")
  out <- mcar_logistic_test(fit, data = dat, formula = ~ treatment + age)
  expect_s3_class(out, "htest")
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
    "intermittent response missingness"
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
  expect_error(mcar_logistic_test(fit, orstr = "ar1"), "should be one of")
  expect_error(mcar_logistic_test(fit, test = "lrt"), "must be one of")
  expect_error(
    mcar_logistic_test(
      fit,
      formula = ~ treatment + age,
      orstr = "exchangeable",
      test = "working-lrt"
    ),
    "independence working model"
  )
})
