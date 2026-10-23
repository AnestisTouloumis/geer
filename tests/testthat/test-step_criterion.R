fit_resp_main <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = test_data$respiratory2,
  link = "probit",
  orstr = "exchangeable",
  method = "pgee-jeffreys"
)

fit_resp_small <- update(fit_resp_main, formula = status ~ baseline)

qicu_trace <- capture.output(
  qicu_search <- step_criterion(fit_resp_main, criterion = "QICu", trace = TRUE)
)

test_that("step_criterion moves to the candidate with the smallest criterion", {
  expect_gt(nrow(qicu_search$anova), 1L)
  expect_greedy_path(qicu_search, fit_resp_main, "QICu", "backward")
  out <- step_criterion(
    fit_resp_small,
    scope = ~ baseline + treatment + gender + visit + age,
    direction = "forward",
    criterion = "AGPC"
  )
  expect_greedy_path(out, fit_resp_small, "AGPC", "forward")
  out <- step_criterion(
    fit_resp_main,
    scope = list(
      lower = ~baseline,
      upper = ~ baseline + treatment + gender + visit + age
    ),
    direction = "both",
    criterion = "QIC",
    cov_type = "naive"
  )
  expect_greedy_path(
    out,
    fit_resp_main,
    "QIC",
    "both",
    lower = "baseline",
    cov_type = "naive"
  )
})

test_that("step_criterion returns the selected fit with a criterion table", {
  out <- step_criterion(fit_resp_main, criterion = "SGPC")
  expect_s3_class(out, "geer")
  expect_s3_class(out$anova, "geer_anova")
  expect_named(out$anova, c("Step", "Df", "SGPC"))
  expect_identical(out$anova$Step[1L], NA_integer_)
  expect_identical(rownames(out$anova)[1L], "Initial model")
  expect_true(all(out$anova$Df[-1L] == 1))
  expect_true(all(diff(out$anova$SGPC) < 0))
  heading <- attr(out$anova, "heading")
  expect_true("SGPC" %in% heading)
  out <- step_criterion(
    fit_resp_main,
    criterion = "QICC",
    cov_type = "bias-corrected",
    steps = 0
  )
  expect_true(
    "QICC with the bias-corrected covariance matrix" %in%
      attr(out$anova, "heading")
  )
})

test_that("step_criterion records the candidate table of every step", {
  out <- qicu_search
  tables <- out$candidates
  expect_length(tables, nrow(out$anova))
  for (i in seq_along(tables)) {
    tab <- tables[[i]]
    expect_s3_class(tab, "geer_anova")
    expect_named(tab, c("Df", "QICu"))
    expect_true("<none>" %in% rownames(tab))
    expect_false(is.unsorted(tab$QICu, na.rm = TRUE))
    expect_equal(tab["<none>", "QICu"], out$anova$QICu[i])
  }
  for (i in seq_len(nrow(out$anova) - 1L)) {
    expect_identical(
      paste0("Step ", i, ": ", rownames(tables[[i]])[1L]),
      rownames(out$anova)[i + 1L]
    )
  }
  expect_identical(rownames(tables[[length(tables)]])[1L], "<none>")
  expect_true(any(grepl("Step 1 candidates", qicu_trace, fixed = TRUE)))
})

test_that("step_criterion prints long formulas on one line in the heading", {
  fit_long <- update(
    fit_resp_main,
    formula = status ~ (baseline + treatment + gender + visit + age)^2
  )
  out <- step_criterion(fit_long, criterion = "QIC", steps = 0)
  heading <- attr(out$anova, "heading")
  final <- heading[which(heading == "\nFinal Model:") + 1L]
  expect_identical(
    final,
    paste(deparse(formula(out), width.cutoff = 500L), collapse = " ")
  )
  expect_false(grepl("  ", final))
})

test_that("step_criterion does not refit candidate models fitted earlier", {
  cache <- new.env(parent = emptyenv())
  key <- build_step_criterion_key(fit_resp_main, "- age")
  assign(key, list(fit = fit_resp_main, value = -1), envir = cache)
  candidate <- fit_step_criterion_candidate(
    fitted_model = fit_resp_main,
    change = "- age",
    criterion = "QICu",
    cov_type = "robust",
    obs_no = fit_resp_main$obs_no,
    cache = cache
  )
  expect_identical(candidate$value, -1)
  fresh <- fit_step_criterion_candidate(
    fitted_model = fit_resp_main,
    change = "- gender",
    criterion = "QICu",
    cov_type = "robust",
    obs_no = fit_resp_main$obs_no,
    cache = cache
  )
  expect_equal(
    fresh$value,
    criterion_value(update(fit_resp_main, . ~ . - gender), "QICu"),
    tolerance = 1e-10
  )
  expect_true(exists(
    build_step_criterion_key(fit_resp_main, "- gender"),
    envir = cache,
    inherits = FALSE
  ))
})

test_that("step_criterion requires an improvement beyond numerical noise", {
  expect_equal(compute_step_tolerance(0.5), sqrt(.Machine$double.eps))
  expect_equal(
    compute_step_tolerance(-200),
    200 * sqrt(.Machine$double.eps)
  )
})

test_that("step_criterion keeps the terms of the lower scope", {
  out <- step_criterion(
    fit_resp_main,
    scope = list(lower = ~ baseline + treatment + gender),
    criterion = "QICu"
  )
  expect_true(all(
    c("baseline", "treatment", "gender") %in%
      attr(terms(formula(out)), "term.labels")
  ))
})

test_that("step_criterion uses cov_type only for QIC and QICC", {
  naive <- step_criterion(
    fit_resp_main,
    criterion = "QICu",
    cov_type = "naive",
    steps = 1
  )
  expect_equal(naive$candidates[[1L]], qicu_search$candidates[[1L]])
  qic_naive <- step_criterion(
    fit_resp_main,
    criterion = "QIC",
    cov_type = "naive",
    steps = 0
  )
  expect_equal(
    qic_naive$anova$QIC,
    criterion_value(fit_resp_main, "QIC", "naive"),
    tolerance = 1e-10
  )
})

test_that("step_criterion stops at the step limit or without a forward scope", {
  out <- step_criterion(fit_resp_main, steps = 0)
  expect_identical(nrow(out$anova), 1L)
  expect_equal(coef(out), coef(fit_resp_main))
  expect_error(
    step_criterion(fit_resp_small, direction = "forward"),
    "'scope' must be supplied for direction = \"forward\"",
    fixed = TRUE
  )
})

test_that("step_criterion stops when the initial criterion is unavailable", {
  d <- ref_sim(
    "cc",
    rep(list(1:4), 4),
    family = gaussian(),
    link = "identity"
  )
  fit <- geewa(
    formula = y ~ x1 + x2,
    data = d$df,
    id = id,
    repeated = visit,
    family = gaussian(link = "identity"),
    corstr = "exchangeable"
  )
  expect_error(
    step_criterion(fit, criterion = "QICC"),
    "step_criterion failed: QICC cannot be evaluated for the initial model"
  )
})

test_that("step_criterion skips candidates whose criterion is unavailable", {
  d <- ref_sim(
    "cc",
    rep(list(1:4), 6),
    family = gaussian(),
    link = "identity"
  )
  d$df$z <- rep(c(-1, 0, 1), length.out = nrow(d$df))
  fit <- geewa(
    formula = y ~ x1 + x2,
    data = d$df,
    id = id,
    repeated = visit,
    family = gaussian(link = "identity"),
    corstr = "exchangeable"
  )
  expect_true(is.finite(criterion_value(fit, "QICC")))
  out <- step_criterion(
    fit,
    scope = ~ x1 + x2 + z,
    direction = "forward",
    criterion = "QICC"
  )
  expect_identical(nrow(out$anova), 1L)
  expect_identical(attr(terms(formula(out)), "term.labels"), c("x1", "x2"))
})

test_that("step_criterion validates its arguments", {
  expect_error(
    step_criterion(fit_resp_main, criterion = "CIC"),
    "CIC is meant for selecting the working association structure"
  )
  expect_error(
    step_criterion(fit_resp_main, criterion = "qicu"),
    "'criterion' must be one of: QICu, QIC, QICC, AGPC, SGPC"
  )
  expect_error(
    step_criterion(fit_resp_main, direction = "sideways"),
    "'direction' must be one of"
  )
  expect_error(
    step_criterion(fit_resp_main, trace = NA),
    "'trace' must be TRUE or FALSE"
  )
})
