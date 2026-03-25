testthat::local_edition(3)

data("respiratory", package = "geer")
respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]

fit_resp_full_indep <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)


test_that("scope validator accepts valid specifications", {
  scope_formula <- ~ baseline + treatment
  scope_char <- "~ baseline + treatment"
  scope_list <- list(lower = ~ baseline, upper = ~ baseline + treatment)
  expect_null(.step_p_validate_scope(NULL))
  expect_identical(.step_p_validate_scope(scope_formula), scope_formula)
  expect_identical(.step_p_validate_scope(scope_char), scope_char)
  expect_identical(.step_p_validate_scope(scope_list), scope_list)
})


test_that("scope validator rejects malformed scope objects", {
  expect_error(.step_p_validate_scope(1), "scope")
  expect_error(.step_p_validate_scope(list()), "lower.*upper|upper.*lower")
  expect_error(.step_p_validate_scope(list(foo = ~ baseline)), "scope")
  expect_error(.step_p_validate_scope(list(lower = 1)), "scope\\$lower")
  expect_error(.step_p_validate_scope(list(upper = 1)), "scope\\$upper")
})


test_that("scope term helper returns factor matrices and empty output for NULL", {
  base_formula <- stats::formula(fit_resp_full_indep)
  expect_identical(.step_p_scope_terms(base_formula, NULL), numeric())
  factors_formula <- .step_p_scope_terms(base_formula, ~ baseline + treatment)
  expect_true(is.matrix(factors_formula))
  expect_true(all(c("baseline", "treatment") %in% colnames(factors_formula)))
  factors_char <- .step_p_scope_terms(base_formula, "~ baseline + treatment")
  expect_equal(dim(factors_char), dim(factors_formula))
  expect_equal(colnames(factors_char), colnames(factors_formula))
})


test_that("scope factor helper handles NULL, formula, and list scopes", {
  model_terms <- stats::terms(fit_resp_full_indep)
  model_factors <- attr(model_terms, "factors")
  out_null <- .step_p_scope_factors(fit_resp_full_indep, NULL, model_terms)
  expect_equal(out_null$add, model_factors)
  expect_identical(out_null$drop, numeric())
  out_formula <- .step_p_scope_factors(
    fit_resp_full_indep,
    ~ baseline + treatment,
    model_terms
  )
  expect_true(is.matrix(out_formula$add))
  expect_true(all(c("baseline", "treatment") %in% colnames(out_formula$add)))
  expect_identical(out_formula$drop, numeric())
  out_list <- .step_p_scope_factors(
    fit_resp_full_indep,
    list(lower = ~ baseline, upper = ~ baseline + treatment),
    model_terms
  )
  expect_true(is.matrix(out_list$add))
  expect_true(is.matrix(out_list$drop))
  expect_true("baseline" %in% colnames(out_list$add))
  expect_true("baseline" %in% colnames(out_list$drop))
})


test_that("step-path initialiser stores the initial model state", {
  out <- .step_p_init_models(fit_resp_full_indep, cov_type = "robust", steps = 3)
  expect_identical(out$steps, 3L)
  expect_length(out$models, 4L)
  expect_identical(out$models_no, 1L)
  expect_identical(out$models[[1L]]$Step, "")
  expect_true(is.na(out$models[[1L]]$Df))
  expect_true(is.finite(out$models[[1L]]$CIC))
})


test_that("append-model helper increments the path correctly", {
  init <- .step_p_init_models(fit_resp_full_indep, cov_type = "robust", steps = 2)
  out <- .step_p_append_model(
    models = init$models,
    models_no = init$models_no,
    step = "- gender",
    df = 1,
    chi = 3.5,
    pval = 0.06,
    cic = 12.3
  )
  expect_identical(out$models_no, 2L)
  expect_identical(out$models[[2L]]$Step, "- gender")
  expect_equal(out$models[[2L]]$Df, 1)
  expect_equal(out$models[[2L]]$Chi, 3.5)
  expect_equal(out$models[[2L]]$`Pr(>Chi)`, 0.06)
  expect_equal(out$models[[2L]]$CIC, 12.3)
})


test_that("result helper stores an anova table on the fitted model", {
  init <- .step_p_init_models(fit_resp_full_indep, cov_type = "robust", steps = 2)
  out_models <- .step_p_append_model(
    models = init$models,
    models_no = init$models_no,
    step = "- gender",
    df = 1,
    chi = 3.5,
    pval = 0.06,
    cic = 12.3
  )
  out <- .step_p_results(
    models = out_models$models[seq_len(out_models$models_no)],
    fit = fit_resp_full_indep,
    object = fit_resp_full_indep
  )
  expect_s3_class(out, "geer")
  expect_s3_class(out$anova, "anova")
  expect_true(is.data.frame(out$anova))
  expect_equal(nrow(out$anova), 2L)
  expect_true(all(c("Step", "Df", "Chi", "Pr(>Chi)", "CIC") %in% names(out$anova)))
  heading <- attr(out$anova, "heading")
  expect_type(heading, "character")
  expect_true(any(grepl("Initial Model:", heading, fixed = TRUE)))
  expect_true(any(grepl("Final Model:", heading, fixed = TRUE)))
})


test_that("update-fit helper applies the requested model change", {
  fit_resp_full_indep <- geewa_binary(
    formula = status ~ baseline + treatment + gender + visit + age,
    id = id,
    repeated = visit,
    data = respiratory2,
    link = "probit",
    orstr = "independence",
    method = "pgee-jeffreys"
  )
  out <- .step_p_update_fit(fit_resp_full_indep, "- gender", obs_no = fit_resp_full_indep$obs_no)
  expect_s3_class(out, "geer")
  expect_false("gender" %in% attr(stats::terms(out), "term.labels"))
  expect_identical(out$obs_no, fit_resp_full_indep$obs_no)
})

test_that("update-fit helper applies the requested model change", {
  fitted_model <- fit_geewa_pois_exch

  updated <- .step_p_update_fit(
    fitted_model = fitted_model,
    change = "- lnage",
    obs_no = fitted_model$obs_no
  )

  expect_s3_class(updated, "geer")
  expect_false("lnage" %in% attr(stats::terms(updated), "term.labels"))
  expect_true(all(c("treatment", "lnbaseline") %in%
                    attr(stats::terms(updated), "term.labels")))
  expect_identical(updated$obs_no, fitted_model$obs_no)
})


test_that("backward selector prioritizes zero-df drops", {
  aod <- data.frame(
    Df = c(NA, 1, 0, 1),
    Chi = c(NA, 2, 0, 3),
    `Pr(>Chi)` = c(NA, 0.40, NA, 0.90),
    check.names = FALSE
  )
  rownames(aod) <- c("<none>", "- x1", "- x2", "- x3")
  expect_identical(.step_p_selected_row_backward(aod), 3L)
})


test_that("backward selector otherwise uses the first maximum p-value", {
  aod <- data.frame(
    Df = c(NA, 1, 1, 1),
    Chi = c(NA, 2, 3, 4),
    `Pr(>Chi)` = c(NA, 0.70, 0.70, 0.20),
    check.names = FALSE
  )
  rownames(aod) <- c("<none>", "- a", "- b", "- c")
  expect_identical(.step_p_selected_row_backward(aod), 2L)
})


test_that("backward stop helper honors zero-df candidates and threshold logic", {
  aod_zdf <- data.frame(
    Df = c(NA, 0, 1),
    Chi = c(NA, 0, 2),
    `Pr(>Chi)` = c(NA, NA, 0.01),
    check.names = FALSE
  )
  rownames(aod_zdf) <- c("<none>", "- x0", "- x1")
  aod_stop <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.01),
    check.names = FALSE
  )
  rownames(aod_stop) <- c("<none>", "- x1")
  aod_continue <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.40),
    check.names = FALSE
  )
  rownames(aod_continue) <- c("<none>", "- x1")
  expect_false(.step_p_backward_should_stop(aod_zdf, 2L, pvalue = 0.15))
  expect_true(.step_p_backward_should_stop(aod_stop, 2L, pvalue = 0.15))
  expect_false(.step_p_backward_should_stop(aod_continue, 2L, pvalue = 0.15))
  expect_true(.step_p_backward_should_stop(aod_continue, NA_integer_, pvalue = 0.15))
})


test_that("forward selector uses the first minimum p-value", {
  aod <- data.frame(
    Df = c(NA, 1, 1, 1),
    Chi = c(NA, 4, 5, 6),
    `Pr(>Chi)` = c(NA, 0.03, 0.03, 0.20),
    check.names = FALSE
  )
  rownames(aod) <- c("<none>", "+ x1", "+ x2", "+ x3")
  expect_identical(.step_p_selected_row_forward(aod), 2L)
})


test_that("forward stop helper stops only when the selected p-value is above threshold", {
  aod_stop <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.30),
    check.names = FALSE
  )
  rownames(aod_stop) <- c("<none>", "+ x1")
  aod_continue <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.01),
    check.names = FALSE
  )
  rownames(aod_continue) <- c("<none>", "+ x1")
  row_stop <- .step_p_selected_row_forward(aod_stop)
  row_continue <- .step_p_selected_row_forward(aod_continue)
  expect_true(.step_p_forward_should_stop(aod_stop, row_stop, pvalue = 0.15))
  expect_false(.step_p_forward_should_stop(aod_continue, row_continue, pvalue = 0.15))
  expect_true(.step_p_forward_should_stop(aod_continue, NA_integer_, pvalue = 0.15))
})


test_that("step label and term helpers return the expected values", {
  aod <- data.frame(Df = 1, row.names = "+ x1")
  expect_identical(.step_p_selected_step_label(aod, 1L), "+ x1")
  expect_identical(.step_p_selected_step_label(aod, 2L), NA_character_)
  expect_identical(.step_p_selected_term("+ x1"), "x1")
  expect_identical(.step_p_selected_term("- x2"), "x2")
  expect_identical(.step_p_selected_term(NA_character_), NA_character_)
  expect_identical(.step_p_selected_term(""), NA_character_)
})


test_that("cycle detection depends on the underlying term", {
  expect_true(.step_p_is_cycle("+ x1", "- x1"))
  expect_true(.step_p_is_cycle("- x1", "+ x1"))
  expect_false(.step_p_is_cycle("+ x1", "- x2"))
  expect_false(.step_p_is_cycle("+ x1", NA_character_))
})
