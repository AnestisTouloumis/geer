testthat::local_edition(3)

data("respiratory", package = "geer")
respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
cerebrovascular <- test_data$cerebrovascular

fit_resp_full_indep <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)

fit_resp_full_exch <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = respiratory2,
  link = "probit",
  orstr = "exchangeable",
  method = "pgee-jeffreys"
)

fit_resp_lower <- geewa_binary(
  formula = status ~ baseline + treatment,
  id = id,
  repeated = visit,
  data = respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)

fit_resp_mid <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit,
  id = id,
  repeated = visit,
  data = respiratory2,
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
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_bin_lower <- geewa(
  formula = ecg ~ treatment,
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
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
    step_p(fit_resp_full_indep, scope = 1, direction = "both"),
    "scope"
  )
  expect_error(
    step_p(
      fit_resp_full_indep,
      scope = list(foo = ~ baseline + treatment),
      direction = "both"
    ),
    "scope"
  )
  expect_error(
    step_p(
      fit_resp_full_indep,
      scope = list(lower = 1),
      direction = "both"
    ),
    "scope\\$lower"
  )
})


test_that("step_p accepts a list scope with upper only", {
  out <- step_p(
    fit_resp_lower,
    scope = list(upper = ~ baseline + treatment + gender + visit + age),
    direction = "forward",
    test = "score",
    cov_type = "robust",
    p_enter = 0.20,
    steps = 3
  )
  expect_step_p_result(out)
  expect_true(all(c("baseline", "treatment") %in% term_labels(out)))
  expect_true(all(term_labels(out) %in% c("baseline", "treatment", "gender", "visit", "age")))
})


test_that("step_p treats missing scope and explicit NULL scope equivalently", {
  out_missing <- step_p(
    fit_resp_full_indep,
    direction = "backward",
    test = "wald",
    cov_type = "robust",
    p_remove = 0.20,
    steps = 3
  )
  out_null <- step_p(
    fit_resp_full_indep,
    scope = NULL,
    direction = "backward",
    test = "wald",
    cov_type = "robust",
    p_remove = 0.20,
    steps = 3
  )
  expect_step_p_result(out_missing)
  expect_step_p_result(out_null)
  expect_equal(term_labels(out_missing), term_labels(out_null))
  expect_equal(out_missing$anova$Step, out_null$anova$Step)
  expect_equal(out_missing$anova$Df, out_null$anova$Df)
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
  expect_step_p_result(out)
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
  expect_step_p_result(out)
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
  expect_step_p_result(out)
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
  expect_step_p_result(out)
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


test_that("step_p stores a readable anova heading on the returned fit", {
  out <- step_p(
    fit_resp_full_indep,
    direction = "backward",
    test = "wald",
    cov_type = "robust",
    p_remove = 0.20,
    steps = 1
  )
  expect_step_p_result(out)
  heading <- attr(out$anova, "heading")
  expect_type(heading, "character")
  expect_true(length(heading) >= 4L)
  expect_true(any(grepl("Initial Model:", heading, fixed = TRUE)))
  expect_true(any(grepl("Final Model:", heading, fixed = TRUE)))
})


test_that("step_p stores numbered steps and descriptive row names", {
  out <- step_p(
    fit_resp_full_indep,
    direction = "backward",
    test = "wald",
    cov_type = "robust",
    p_remove = 0.20,
    steps = 3
  )
  expect_s3_class(out$anova, "anova")
  expect_identical(rownames(out$anova)[1L], "Initial model")
  expect_identical(out$anova$Step[1L], "")

  if (nrow(out$anova) > 1L) {
    expect_identical(
      out$anova$Step[-1L],
      as.character(seq_len(nrow(out$anova) - 1L))
    )
    expect_true(all(grepl("^Step [0-9]+: [+-] ", rownames(out$anova)[-1L])))
  }
})


test_that("step_p preserves the original data call for downstream update", {
  fitted_model <- geewa_binary(
    formula = status ~ (I(treatment == "active") + gender + visit + age + baseline)^2,
    id = id,
    repeated = visit,
    data = respiratory2,
    link = "probit",
    orstr = "independence",
    method = "pgee-jeffreys"
  )
  out <- step_p(
    fitted_model,
    scope = list(lower = status ~ 1, upper = status ~ .^2),
    test = "wald",
    cov_type = "bias-corrected",
    p_remove = 0.1,
    direction = "backward"
  )
  expect_s3_class(out, "geer")
  expect_no_error(update(out, orstr = "exchangeable"))
})
