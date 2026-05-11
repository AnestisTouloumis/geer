testthat::local_edition(3)

cerebrovascular <- test_data$cerebrovascular
skip_if_no_emmeans <- function() {
  testthat::skip_if_not_installed("emmeans")
}


test_that("emmeans support methods are available for geer objects", {
  skip_if_no_emmeans()
  expect_true(
    exists("recover_data.geer", where = asNamespace("geer"), inherits = FALSE)
  )
  expect_true(
    exists("emm_basis.geer", where = asNamespace("geer"), inherits = FALSE)
  )
})


test_that("emmeans recovers the predictor data for geer fits", {
  skip_if_no_emmeans()
  rec <- emmeans::recover_data(fit_geewa_bin_exch)
  expect_s3_class(rec, "data.frame")
  expect_true(all(c("treatment", "period") %in% names(rec)))
  expect_false("ecg" %in% names(rec))
})


test_that("emmeans works for geer objects on link and response scales", {
  skip_if_no_emmeans()
  emm_link <- emmeans::emmeans(fit_geewa_bin_exch, ~ treatment)
  emm_resp <- emmeans::emmeans(
    fit_geewa_bin_exch,
    ~ treatment,
    type = "response"
  )
  link_df <- as.data.frame(emm_link)
  resp_df <- as.data.frame(emm_resp)
  expect_s4_class(emm_link, "emmGrid")
  expect_s4_class(emm_resp, "emmGrid")
  expect_true("emmean" %in% names(link_df))
  expect_true(any(c("prob", "response") %in% names(resp_df)))
  expect_true(all(is.finite(link_df$SE)))
  expect_true(all(is.finite(resp_df$SE)))
  expect_true(all(is.infinite(link_df$df)))
})


test_that("naive emmeans from an independence geer fit match glm results", {
  skip_if_no_emmeans()
  fit_geer <- geewa_binary(
    formula = ecg ~ treatment + factor(period),
    id = id,
    data = cerebrovascular,
    link = "probit",
    orstr = "independence",
    method = "gee"
  )
  fit_glm <- glm(
    formula = ecg ~ treatment + factor(period),
    data = cerebrovascular,
    family = binomial(link = "probit")
  )
  emm_geer <- emmeans::emmeans(
    fit_geer,
    ~ treatment,
    vcov.method = "naive"
  )
  emm_glm <- emmeans::emmeans(fit_glm, ~ treatment)
  geer_df <- as.data.frame(emm_geer)
  glm_df <- as.data.frame(emm_glm)
  expect_equal(geer_df$emmean, glm_df$emmean, tolerance = 1e-6)
  expect_equal(geer_df$SE, glm_df$SE, tolerance = 1e-6)
})


test_that("emmeans support honors covariance aliases and custom matrices", {
  skip_if_no_emmeans()
  custom_vcov <- vcov(fit_geewa_bin_exch, cov_type = "naive")
  emm_custom_method <- emmeans::emmeans(
    fit_geewa_bin_exch,
    ~ treatment,
    vcov.method = custom_vcov
  )
  emm_custom_dot <- emmeans::emmeans(
    fit_geewa_bin_exch,
    ~ treatment,
    vcov. = custom_vcov
  )
  custom_method_df <- as.data.frame(emm_custom_method)
  custom_dot_df <- as.data.frame(emm_custom_dot)
  expect_equal(custom_method_df$emmean, custom_dot_df$emmean, tolerance = 1e-10)
  expect_equal(custom_method_df$SE, custom_dot_df$SE, tolerance = 1e-10)
})
