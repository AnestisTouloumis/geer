local({
  data("respiratory", package = "geer")
  fmla <- status ~ baseline + treatment + gender + visit + age + center
  fitted_pgee_corr <- geewa(
    formula = fmla,
    id = id,
    repeated = visit,
    family = stats::binomial(link = "logit"),
    data = respiratory,
    phi_fixed = TRUE,
    phi_value = 1,
    corstr = "independence",
    method = "pgee-jeffreys"
  )
  fitted_pgee_or <- geewa_binary(
    formula = fmla,
    id = id,
    repeated = visit,
    link = "logit",
    data = respiratory,
    orstr = "independence",
    method = "pgee-jeffreys"
  )
  fitted_brnaive_or <- update(fitted_pgee_or, method = "brgee-naive")


  test_that("jeffreys = naive under independence (OR matches COR and brgee-naive)", {
    expect_equal(coef(fitted_pgee_or), coef(fitted_pgee_corr), tolerance = 1e-5)
    expect_equal(coef(fitted_pgee_or), coef(fitted_brnaive_or), tolerance = 1e-5)
  })
})
