testthat::local_edition(3)

fit_resp_full_indep <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = test_data$respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)

fit_resp_full_exch <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = test_data$respiratory2,
  link = "probit",
  orstr = "exchangeable",
  method = "pgee-jeffreys"
)

fit_resp_lower <- geewa_binary(
  formula = status ~ baseline + treatment,
  id = id,
  repeated = visit,
  data = test_data$respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)

fit_resp_mid <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit,
  id = id,
  repeated = visit,
  data = test_data$respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)

fit_bin_full_indep <- geewa(
  formula = ecg ~ treatment + factor(period),
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = test_data$cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_bin_lower <- geewa(
  formula = ecg ~ treatment,
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = test_data$cerebrovascular,
  corstr = "independence",
  method = "gee"
)

test_that("step_p rejects invalid scalar inputs", {
  expect_error(
    step_p(fit_resp_full_indep, p_enter = 0, direction = "backward"),
    "p_enter"
  )
  expect_error(
    step_p(fit_resp_full_indep, p_remove = "0.1", direction = "backward"),
    "p_remove"
  )
  expect_error(
    step_p(fit_resp_full_indep, steps = 1.5, direction = "backward"),
    "steps"
  )
})

test_that("step_p rejects invalid scope specifications", {
  expect_error(
    step_p(
      fit_resp_full_indep,
      scope = list(upper = ~ baseline + treatment),
      direction = "both"
    ),
    "scope"
  )

  expect_error(
    step_p(
      fit_resp_full_indep,
      scope = list(
        lower = ~ baseline + treatment + age,
        upper = ~ baseline + treatment
      ),
      direction = "both"
    ),
    "scope"
  )
})

test_that("step_p rejects working-lrt for non-independence association structures", {
  expect_error(
    step_p(
      fit_resp_full_exch,
      direction = "backward",
      test = "working-lrt"
    ),
    "independence working model"
  )
})

test_that("step_p backward returns a valid path within the step limit", {
  start_terms <- term_labels(fit_resp_full_indep)

  out <- step_p(
    fit_resp_full_indep,
    direction = "backward",
    test = "wald",
    cov_type = "robust",
    p_remove = 0.20,
    steps = 3
  )

  expect_step_p_path(out)

  final_terms <- term_labels(out)
  expect_true(all(final_terms %in% start_terms))
  expect_true(length(final_terms) <= length(start_terms))
  expect_true(nrow(out$anova) <= 4L)
})

test_that("step_p backward respects the lower scope", {
  out <- step_p(
    object = fit_bin_full_indep,
    scope = list(lower = ~ treatment, upper = ~ treatment + factor(period)),
    direction = "backward",
    test = "wald",
    cov_type = "naive",
    p_remove = 0.99,
    steps = 3
  )

  expect_step_p_path(out)
  expect_true("treatment" %in% term_labels(out))
  expect_true(all(term_labels(out) %in% c("treatment", "factor(period)")))
})

test_that("step_p forward respects the upper scope", {
  out <- step_p(
    fit_resp_lower,
    scope = ~ baseline + treatment + gender + visit + age,
    direction = "forward",
    test = "score",
    cov_type = "robust",
    p_enter = 0.20,
    steps = 3
  )

  expect_step_p_path(out)

  final_terms <- term_labels(out)
  allowed_terms <- c("baseline", "treatment", "gender", "visit", "age")

  expect_true(all(c("baseline", "treatment") %in% final_terms))
  expect_true(all(final_terms %in% allowed_terms))
  expect_true(nrow(out$anova) <= 4L)
})

test_that("step_p forward can stop without taking a step", {
  out <- step_p(
    object = fit_bin_lower,
    scope = list(lower = ~ treatment, upper = ~ treatment + factor(period)),
    direction = "forward",
    test = "wald",
    cov_type = "naive",
    p_enter = 1e-10,
    steps = 3
  )

  expect_no_step_taken(out)
  expect_equal(term_labels(out), term_labels(fit_bin_lower))
})

test_that("step_p both-direction respects explicit lower and upper scope", {
  out <- step_p(
    fit_resp_mid,
    scope = list(
      lower = ~ baseline + treatment,
      upper = ~ baseline + treatment + gender + visit + age
    ),
    direction = "both",
    test = "working-wald",
    cov_type = "robust",
    pmethod = "rao-scott",
    p_enter = 0.20,
    p_remove = 0.20,
    steps = 4
  )

  expect_step_p_path(out)

  final_terms <- term_labels(out)
  lower_terms <- c("baseline", "treatment")
  upper_terms <- c("baseline", "treatment", "gender", "visit", "age")

  expect_true(all(lower_terms %in% final_terms))
  expect_true(all(final_terms %in% upper_terms))
  expect_true(nrow(out$anova) <= 5L)
})

test_that("step_p returns the original model when steps = 0", {
  out <- step_p(
    fit_resp_full_indep,
    scope = ~ baseline + treatment + gender + visit + age,
    direction = "backward",
    test = "wald",
    cov_type = "robust",
    p_remove = 0.20,
    steps = 0
  )

  expect_no_step_taken(out)
  expect_equal(term_labels(out), term_labels(fit_resp_full_indep))
})
