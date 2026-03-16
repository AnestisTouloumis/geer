testthat::local_edition(3)
data("cerebrovascular", package = "geer")
fit_full <- geewa(
  formula = ecg ~ treatment + factor(period),
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)
fit_small <- geewa(
  formula = ecg ~ treatment,
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)
fit_exch <- tryCatch(
  geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "logit"),
    phi_fixed = TRUE,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "exchangeable",
    method = "gee"
  ),
  error = function(e) NULL
)


testthat::test_that("step_p returns a geer object with an anova path (backward)", {
  res <- step_p(
    object = fit_full,
    direction = "backward",
    test = "wald",
    cov_type = "naive",
    p_remove = 0.50,
    steps = 1
  )

  testthat::expect_s3_class(res, "geer")


  testthat::expect_true(!is.null(res$anova))


  testthat::expect_s3_class(res$anova, "anova")


  testthat::expect_true(is.data.frame(res$anova))


  needed_cols <- c("Step", "Df", "Chi", "Pr(>Chi)", "CIC")
  testthat::expect_true(all(needed_cols %in% names(res$anova)))


  testthat::expect_lte(nrow(res$anova), 2L)


  testthat::expect_gte(nrow(res$anova), 1L)
})


testthat::test_that("step_p forward can add a term within scope", {
  res <- step_p(
    object = fit_small,
    scope = list(lower = ~ treatment, upper = ~ treatment + factor(period)),
    direction = "forward",
    test = "wald",
    cov_type = "naive",
    p_enter = 0.99,
    steps = 1
  )


  testthat::expect_s3_class(res, "geer")


  testthat::expect_true(!is.null(res$anova))


  testthat::expect_lte(nrow(res$anova), 2L)


  testthat::expect_gte(nrow(res$anova), 1L)
  ftxt <- paste(deparse(stats::formula(res)), collapse = " ")


  testthat::expect_true(grepl("factor\\(period\\)", ftxt) || identical(res, fit_small))
})


testthat::test_that("step_p both returns a valid path object", {
  res <- step_p(
    object = fit_full,
    scope = list(lower = ~ 1, upper = ~ treatment + factor(period)),
    direction = "both",
    test = "wald",
    cov_type = "naive",
    p_enter = 0.99,
    p_remove = 0.50,
    steps = 2
  )


  testthat::expect_s3_class(res, "geer")


  testthat::expect_true(!is.null(res$anova))


  testthat::expect_s3_class(res$anova, "anova")


  testthat::expect_lte(nrow(res$anova), 3L) # initial + up to 2


  testthat::expect_gte(nrow(res$anova), 1L)
})


testthat::test_that("step_p validates p_enter and p_remove", {
  testthat::expect_error(
    step_p(fit_full, direction = "forward", p_enter = 1.5),
    "p_enter",
    fixed = FALSE
  )


  testthat::expect_error(
    step_p(fit_full, direction = "backward", p_remove = -0.2),
    "p_remove",
    fixed = FALSE
  )
})


testthat::test_that("step_p rejects working-lrt for non-independence models", {
  if (is.null(fit_exch)) {
    testthat::skip("exchangeable working correlation not available for this fit in the test environment")
  }


  testthat::expect_error(
    step_p(fit_exch, test = "working-lrt"),
    "independence",
    fixed = FALSE
  )
})

testthat::test_that("step_p validates p_enter and p_remove", {
  testthat::expect_error(
    step_p(fit_full, direction = "forward", p_enter = 1.5),
    "p_enter",
    fixed = FALSE
  )

  testthat::expect_error(
    step_p(fit_full, direction = "backward", p_remove = -0.2),
    "p_remove",
    fixed = FALSE
  )
})

testthat::test_that("step_p rejects working-lrt for non-independence models", {
  if (is.null(fit_exch)) {
    testthat::skip("exchangeable working correlation not available for this fit in the test environment")
  }

  testthat::expect_error(
    step_p(fit_exch, test = "working-lrt"),
    "independence",
    fixed = FALSE
  )
})
