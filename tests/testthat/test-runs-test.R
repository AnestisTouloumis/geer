testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch


count_runs <- function(signs) {
  signs <- signs[signs != 0]
  1L + sum(signs[-1L] != signs[-length(signs)])
}


test_that("runs-test calculation matches the Chang sign-sequence example", {
  signs <- c(1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1)
  out <- geer:::compute_runs_statistics(signs, "two.sided")

  expected_runs <- 71 / 11
  variance_runs <- 294 / 121
  expected_z <- (6 - expected_runs) / sqrt(variance_runs)

  expect_identical(out$runs, 6L)
  expect_identical(out$positive, 5L)
  expect_identical(out$negative, 6L)
  expect_identical(out$zero, 0L)
  expect_equal(out$expected_runs, expected_runs, tolerance = 1e-12)
  expect_equal(out$variance_runs, variance_runs, tolerance = 1e-12)
  expect_equal(out$statistic, expected_z, tolerance = 1e-12)
  expect_equal(
    out$p_value,
    2 * stats::pnorm(-abs(expected_z)),
    tolerance = 1e-12
  )
})


test_that("the statistic does not depend on the residual type", {
  ## Working, Pearson and deviance residuals differ only by strictly positive
  ## scale factors (or a signed square root), so they share a sign sequence.
  ord <- order(
    count_fit$id,
    count_fit$repeated,
    seq_len(count_fit$obs_no)
  )
  signs <- lapply(
    c("working", "pearson", "deviance"),
    function(type) sign(residuals(count_fit, type = type)[ord])
  )
  expect_identical(signs[[2L]], signs[[1L]])
  expect_identical(signs[[3L]], signs[[1L]])

  stats <- lapply(
    signs,
    function(x) geer:::compute_runs_statistics(x, "two.sided")
  )
  expect_identical(stats[[2L]], stats[[1L]])
  expect_identical(stats[[3L]], stats[[1L]])
})


test_that("runs_test uses natural cluster/repeated ordering by default", {
  out <- runs_test(count_fit)
  residual_values <- residuals(count_fit, type = "working")
  ord <- order(
    count_fit$id,
    count_fit$repeated,
    seq_along(residual_values)
  )
  signs <- sign(residual_values[ord])

  expect_s3_class(out, "htest")
  expect_identical(out$alternative, "two.sided")
  expect_identical(out$order_by, "natural cluster/repeated order")
  expect_identical(out$parameter, c(n_p = out$positive, n_n = out$negative))
  expect_identical(out$estimate, c("number of runs" = as.numeric(out$runs)))
  expect_identical(out$null.value, c("number of runs" = out$expected_runs))
  expect_identical(out$runs, as.integer(count_runs(signs)))
  expect_identical(out$positive, as.integer(sum(signs > 0)))
  expect_identical(out$negative, as.integer(sum(signs < 0)))
  expect_identical(out$zero, as.integer(sum(signs == 0)))
  expect_equal(
    unname(out$statistic),
    (out$runs - out$expected_runs) / sqrt(out$variance_runs),
    tolerance = 1e-12
  )
  expect_equal(
    out$p.value,
    2 * stats::pnorm(-abs(unname(out$statistic))),
    tolerance = 1e-12
  )
})


test_that("one-sided alternatives use the corresponding normal tail", {
  two_sided <- runs_test(count_fit)
  less <- runs_test(count_fit, alternative = "less")
  greater <- runs_test(count_fit, alternative = "greater")
  z <- unname(two_sided$statistic)

  expect_identical(less$alternative, "less")
  expect_identical(greater$alternative, "greater")
  expect_identical(unname(less$statistic), z)
  expect_identical(unname(greater$statistic), z)
  expect_equal(less$p.value, stats::pnorm(z), tolerance = 1e-12)
  expect_equal(
    greater$p.value,
    stats::pnorm(z, lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(less$p.value + greater$p.value, 1, tolerance = 1e-12)
  expect_equal(
    two_sided$p.value,
    2 * min(less$p.value, greater$p.value),
    tolerance = 1e-12
  )
})


test_that("the Hardin and Hilbe worked example is reproduced one-sided", {
  ## Hardin and Hilbe (2013), equations (4.35)-(4.41): n_p = 42, n_n = 38,
  ## T = 44, E(T) = 40.9, V(T) = 19.65, Z = 0.6993 and a one-sided p of .2422.
  signs <- c(rep(c(1, -1), times = 38), rep(1, 4))
  out <- geer:::compute_runs_statistics(signs, "greater")

  expect_identical(out$positive, 42L)
  expect_identical(out$negative, 38L)
  expect_equal(out$expected_runs, 40.9, tolerance = 1e-12)
  expect_equal(out$variance_runs, 19.65, tolerance = 1e-3)

  z <- (44 - out$expected_runs) / sqrt(out$variance_runs)
  expect_equal(z, 0.6993, tolerance = 1e-3)
  expect_equal(stats::pnorm(z, lower.tail = FALSE), 0.2422, tolerance = 1e-3)
})


test_that("compute_runs_statistics rejects an unknown alternative", {
  expect_error(
    geer:::compute_runs_statistics(c(1, -1, 1), "left.sided"),
    "'alternative' must be one of",
    fixed = TRUE
  )
  expect_error(
    geer:::compute_runs_statistics(c(1, -1, 1), c("less", "greater")),
    "'alternative' must be a single character value",
    fixed = TRUE
  )
  expect_error(
    geer:::compute_runs_statistics(c(1, -1, 1), NA_character_),
    "'alternative' must be a single character value",
    fixed = TRUE
  )
})


test_that("compute_runs_statistics returns the retained signs it counted", {
  out <- geer:::compute_runs_statistics(c(1, 0, -1, -1, 2), "two.sided")
  expect_identical(out$retained, c(TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_identical(out$signs, c(1L, -1L, -1L, 1L))
  expect_identical(out$zero, 1L)
  expect_identical(out$nonzero, 4L)
  expect_identical(out$runs, 3L)
})


test_that("runs_test can order residuals by fitted values", {
  out <- runs_test(count_fit, order_by = "fitted")
  residual_values <- residuals(count_fit, type = "working")
  natural <- order(count_fit$id, count_fit$repeated, seq_along(residual_values))
  ord <- natural[order(count_fit$fitted.values[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "fitted values")
  expect_identical(out$runs, as.integer(count_runs(signs)))
  expect_identical(out$positive, as.integer(sum(signs > 0)))
  expect_identical(out$negative, as.integer(sum(signs < 0)))
})


test_that("runs_test can order residuals by a model-matrix covariate", {
  out <- runs_test(count_fit, order_by = "lnage")
  residual_values <- residuals(count_fit, type = "working")
  natural <- order(count_fit$id, count_fit$repeated, seq_along(residual_values))
  key <- count_fit$x[, "lnage"]
  ord <- natural[order(key[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "model-matrix column 'lnage'")
  expect_identical(out$runs, as.integer(count_runs(signs)))
})


test_that("order_by can name a covariate outside the model matrix", {
  ## 'visit' supplies the repeated index and is not a model-matrix column, so
  ## it has to be resolved from the data used to fit the model. With
  ## repeated = visit the resolved ordering must agree with the stored
  ## repeated index.
  fit <- geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy,
    id = id,
    repeated = visit,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_false("visit" %in% colnames(fit$x))

  out <- runs_test(fit, order_by = "visit")
  expect_identical(out$order_by, "variable 'visit'")

  reference <- runs_test(fit, order_by = fit$repeated)
  expect_identical(out$runs, reference$runs)
  expect_identical(unname(out$statistic), unname(reference$statistic))
})


test_that("runs_test accepts a supplied ordering vector", {
  key <- rep(c(2, 1, 3), length.out = count_fit$obs_no)
  out <- runs_test(count_fit, order_by = key)
  residual_values <- residuals(count_fit, type = "working")
  natural <- order(count_fit$id, count_fit$repeated, seq_along(residual_values))
  ord <- natural[order(key[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "supplied ordering vector")
  expect_identical(out$runs, as.integer(count_runs(signs)))
})


test_that("zero residuals are omitted from the sign sequence", {
  fit <- count_fit
  fit$residuals <- rep(c(1, 0, -1, 0), length.out = fit$obs_no)
  out <- runs_test(fit)

  expect_identical(out$zero, as.integer(sum(fit$residuals == 0)))
  expect_identical(out$nonzero, fit$obs_no - out$zero)
  expect_identical(out$positive + out$negative, out$nonzero)
})


test_that("the result carries the tested sign sequence", {
  out <- runs_test(count_fit)

  expect_s3_class(out, "geer_runs_test")
  expect_s3_class(out, "htest")
  expect_true(out$natural_order)
  expect_length(out$signs, out$nonzero)
  expect_length(out$cluster, out$nonzero)
  expect_true(all(out$signs %in% c(-1L, 1L)))
  expect_identical(sum(out$signs > 0L), out$positive)
  expect_identical(sum(out$signs < 0L), out$negative)
  expect_identical(
    1L + sum(out$signs[-1L] != out$signs[-out$nonzero]),
    out$runs
  )
  expect_identical(
    as.numeric(sort(unique(out$cluster))),
    as.numeric(sort(unique(count_fit$id)))
  )
})


test_that("plot.geer_runs_test accepts point arguments", {
  out <- runs_test(count_fit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot(out, col = "black", pch = 20L, cex = 0.5))
  expect_silent(plot(out, run_colors = c("red", "blue", "green")))
  expect_error(plot(out, run_colors = character(0)), "at least one colour")
})


test_that("natural_order reports the ordering obtained, not the one requested", {
  ## An ordering with a single distinct value is entirely tied, so the
  ## tie-break returns the natural order and the flag must say so.
  tied <- runs_test(count_fit, order_by = rep(1, count_fit$obs_no))
  natural <- runs_test(count_fit)

  expect_true(tied$natural_order)
  expect_identical(tied$signs, natural$signs)
  expect_identical(tied$runs, natural$runs)
  expect_identical(tied$order_by, "supplied ordering vector")

  ## An ordering that genuinely reorders the residuals must not.
  expect_false(runs_test(count_fit, order_by = "fitted")$natural_order)
})


test_that("plot.geer_runs_test guards cluster breaks and its input", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  by_fitted <- runs_test(count_fit, order_by = "fitted")
  expect_false(by_fitted$natural_order)
  expect_silent(plot(by_fitted))
  expect_warning(
    plot(by_fitted, cluster_breaks = TRUE),
    "only interpretable under the natural"
  )

  natural <- runs_test(count_fit)
  expect_error(
    plot(natural, cluster_breaks = "yes"),
    "single non-missing logical"
  )

  stripped <- natural
  stripped$signs <- NULL
  expect_error(plot(stripped), "tested sign sequence")
})


test_that("runs_test validates ordering and residual sign requirements", {
  expect_error(
    runs_test(count_fit, order_by = "not-a-variable"),
    "unknown 'order_by'"
  )
  expect_error(
    runs_test(count_fit, order_by = seq_len(count_fit$obs_no - 1L)),
    "one value per fitted observation"
  )

  bad_order <- seq_len(count_fit$obs_no)
  bad_order[1] <- NA_real_
  expect_error(
    runs_test(count_fit, order_by = bad_order),
    "finite and non-missing"
  )

  fit <- count_fit
  fit$residuals <- rep(1, fit$obs_no)
  expect_error(
    runs_test(fit),
    "at least one positive and one negative residual"
  )
})


test_that("large sign counts keep the normal approximation", {
  out <- runs_test(count_fit)
  expect_equal(
    out$p.value,
    2 * stats::pnorm(-abs(unname(out$statistic))),
    tolerance = 1e-12
  )
})


test_that("runs_test warns when the normal approximation is based on small sign counts", {
  fit <- count_fit
  fit$residuals <- c(
    rep(1, 10),
    rep(-1, 10),
    rep(0, fit$obs_no - 20)
  )

  expect_warning(
    out <- runs_test(fit),
    "normal approximation may be unreliable"
  )
  expect_identical(out$positive, 10L)
  expect_identical(out$negative, 10L)
})


test_that("plot.geer_runs_test returns the plotted sequence invisibly", {
  out <- runs_test(count_fit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted <- plot(out)
  expect_s3_class(plotted, "data.frame")
  expect_identical(nrow(plotted), out$nonzero)
  expect_identical(names(plotted), c("position", "sign", "run", "cluster"))
  expect_identical(plotted$sign, out$signs)
  expect_identical(plotted$position, seq_len(out$nonzero))
  ## Runs are numbered consecutively and there are exactly T of them.
  expect_identical(max(plotted$run), out$runs)
  expect_identical(plotted$run[[1L]], 1L)
  expect_true(all(diff(plotted$run) %in% c(0L, 1L)))
  ## A colour change happens exactly where the sign changes.
  expect_identical(
    which(diff(plotted$run) == 1L),
    which(plotted$sign[-1L] != plotted$sign[-nrow(plotted)])
  )
})
