testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch
plot_runs_method <- utils::getS3method("plot", "geer_runs_test")


test_that("plot.geer_runs_test only accepts runs-test results", {
  expect_error(
    plot_runs_method(list(signs = c(1, -1, 1))),
    "'x' must be a 'geer_runs_test' object"
  )
  expect_error(
    plot_runs_method(structure(list(), class = "htest")),
    "'x' must be a 'geer_runs_test' object"
  )
})


test_that("plot.geer_runs_test detects a malformed cluster vector", {
  out <- runs_test(count_fit)
  out$cluster <- out$cluster[-1L]

  expect_error(
    plot(out),
    "'x' is malformed: 'cluster' must have one value per tested residual"
  )
})


test_that("plot.geer_runs_test draws results with omitted zero residuals", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  fit <- count_fit
  fit$residuals <- rep(c(1, 0, -1, 0, -1), length.out = fit$obs_no)
  out <- runs_test(fit)
  expect_gt(out$zero, 0L)

  expect_silent(plot(out))
  expect_silent(plot(out, sub = "custom subtitle"))
})


test_that("plot.geer_runs_test requires every extra argument to be named", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- runs_test(count_fit)

  ## The ten formals after 'x' consume the first ten positional values; an
  ## eleventh would be an unnamed graphical parameter.
  expect_error(
    plot_runs_method(
      out, NULL, NULL, c("red", "blue", "green", "orange"), 1L,
      NULL, NULL, NULL, NULL, NULL, "unnamed"
    ),
    "arguments in '...' must all be named"
  )
})


test_that("plot.geer_runs_test passes named graphical parameters through", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- runs_test(count_fit)

  expect_silent(plot(out, main = "Runs", xlab = "Index", ylab = "Sign", las = 1))
})
