plot_runs_method <- getS3method("plot", "geer_runs_test")

test_that("plot only accepts runs-test results", {
  expect_error(
    plot_runs_method(list(signs = c(1, -1, 1))),
    "'x' must be a 'geer_runs_test' object"
  )
  expect_error(
    plot_runs_method(structure(list(), class = "htest")),
    "'x' must be a 'geer_runs_test' object"
  )
})

test_that("plot detects a malformed cluster vector", {
  out <- runs_test(fit_geewa_pois_exch)
  out$cluster <- out$cluster[-1L]

  expect_error(
    plot(out),
    "'x' is malformed: 'cluster' must have one value per tested residual"
  )
})

test_that("plot draws results with omitted zero residuals", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  fit <- fit_geewa_pois_exch
  fit$residuals <- rep(c(1, 0, -1, 0, -1), length.out = fit$obs_no)
  out <- runs_test(fit)
  expect_gt(out$zero, 0L)

  expect_silent(plot(out))
  expect_silent(plot(out, sub = "custom subtitle"))
})

test_that("plot requires every extra argument to be named", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  out <- runs_test(fit_geewa_pois_exch)

  expect_error(
    plot_runs_method(
      out, NULL, NULL, c("red", "blue", "green", "orange"), 1L,
      NULL, NULL, NULL, NULL, NULL, "unnamed"
    ),
    "arguments in '...' must all be named"
  )
})

test_that("plot passes named graphical parameters through", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  out <- runs_test(fit_geewa_pois_exch)

  expect_silent(plot(out, main = "Runs", xlab = "Index", ylab = "Sign", las = 1))
})

test_that("plot accepts point arguments", {
  out <- runs_test(fit_geewa_pois_exch)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot(out, col = "black", pch = 20L, cex = 0.5))
  expect_silent(plot(out, run_colors = c("red", "blue", "green")))
  expect_error(
    plot(out, run_colors = character(0)),
    "must contain at least three colors"
  )
  expect_error(
    plot(out, run_colors = c("red", "blue")),
    "must contain at least three colors"
  )
  expect_silent(plot(out, col = "black", run_colors = "red"))
})

test_that("plot guards cluster breaks and its input", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  by_fitted <- runs_test(fit_geewa_pois_exch, order_by = "fitted")
  expect_false(by_fitted$natural_order)
  expect_silent(plot(by_fitted))
  expect_silent(plot(by_fitted, ylim = c(-2, 2)))
  expect_silent(plot(by_fitted, yaxt = "s"))
  expect_error(plot(by_fitted, type = "p"), "'type' must not be supplied")
  expect_error(
    plot(by_fitted, cluster_breaks = TRUE),
    paste0(
      "'cluster_breaks' must be FALSE or NULL unless 'x' uses the natural ",
      "cluster/repeated ordering"
    ),
    fixed = TRUE
  )

  natural <- runs_test(fit_geewa_pois_exch)
  expect_error(
    plot(natural, cluster_breaks = "yes"),
    "must be NULL, TRUE or FALSE"
  )

  stripped <- natural
  stripped$signs <- NULL
  expect_error(plot(stripped), "tested sign sequence")
})

test_that("plot returns the plotted sequence invisibly", {
  out <- runs_test(fit_geewa_pois_exch)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  plotted <- plot(out)
  expect_s3_class(plotted, "data.frame")
  expect_identical(nrow(plotted), out$nonzero)
  expect_named(plotted, c("position", "sign", "run", "cluster"))
  expect_identical(plotted$sign, out$signs)
  expect_identical(plotted$position, seq_len(out$nonzero))
  expect_identical(max(plotted$run), out$runs)
  expect_identical(plotted$run[[1L]], 1L)
  expect_true(all(diff(plotted$run) %in% c(0L, 1L)))
  expect_identical(
    which(diff(plotted$run) == 1L),
    which(plotted$sign[-1L] != plotted$sign[-nrow(plotted)])
  )
})
