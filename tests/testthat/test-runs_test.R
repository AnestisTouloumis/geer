list2env(get_cached_test_fixtures("runs_test"), environment())

test_that("runs-test calculation matches the Chang sign-sequence example", {
  signs <- c(1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1)
  out <- compute_runs_statistics(signs, "two.sided")

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
    2 * pnorm(-abs(expected_z)),
    tolerance = 1e-12
  )
})

test_that("the statistic does not depend on the residual type", {
  ord <- order(
    fit_geewa_pois_exch$id,
    fit_geewa_pois_exch$repeated,
    seq_len(fit_geewa_pois_exch$obs_no)
  )
  signs <- lapply(
    c("working", "pearson", "deviance"),
    function(type) sign(residuals(fit_geewa_pois_exch, type = type)[ord])
  )
  expect_identical(signs[[2L]], signs[[1L]])
  expect_identical(signs[[3L]], signs[[1L]])

  stats <- lapply(
    signs,
    function(x) compute_runs_statistics(x, "two.sided")
  )
  expect_identical(stats[[2L]], stats[[1L]])
  expect_identical(stats[[3L]], stats[[1L]])
})

test_that("runs_test uses natural cluster/repeated ordering by default", {
  out <- runs_test(fit_geewa_pois_exch)
  residual_values <- residuals(fit_geewa_pois_exch, type = "working")
  ord <- order(
    fit_geewa_pois_exch$id,
    fit_geewa_pois_exch$repeated,
    seq_along(residual_values)
  )
  signs <- sign(residual_values[ord])

  expect_s3_class(out, "htest")
  expect_identical(out$alternative, "two.sided")
  expect_identical(out$order_by, "natural cluster/repeated order")
  expect_identical(out$parameter, c(n_p = out$positive, n_n = out$negative))
  expect_identical(out$estimate, c("number of runs" = as.numeric(out$runs)))
  expect_identical(out$null.value, c("number of runs" = out$expected.runs))
  expect_identical(out$runs, as.integer(count_runs(signs)))
  expect_identical(out$positive, as.integer(sum(signs > 0)))
  expect_identical(out$negative, as.integer(sum(signs < 0)))
  expect_identical(out$zero, as.integer(sum(signs == 0)))
  expect_equal(
    unname(out$statistic),
    (out$runs - out$expected.runs) / sqrt(out$variance.runs),
    tolerance = 1e-12
  )
  expect_equal(
    out$p.value,
    2 * pnorm(-abs(unname(out$statistic))),
    tolerance = 1e-12
  )
})

test_that("one-sided alternatives use the corresponding normal tail", {
  two_sided <- runs_test(fit_geewa_pois_exch)
  less <- runs_test(fit_geewa_pois_exch, alternative = "less")
  greater <- runs_test(fit_geewa_pois_exch, alternative = "greater")
  z <- unname(two_sided$statistic)

  expect_identical(less$alternative, "less")
  expect_identical(greater$alternative, "greater")
  expect_identical(unname(less$statistic), z)
  expect_identical(unname(greater$statistic), z)
  expect_equal(less$p.value, pnorm(z), tolerance = 1e-12)
  expect_equal(
    greater$p.value,
    pnorm(z, lower.tail = FALSE),
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
  signs <- c(rep(c(1, -1), times = 38), rep(1, 4))
  out <- compute_runs_statistics(signs, "greater")

  expect_identical(out$positive, 42L)
  expect_identical(out$negative, 38L)
  expect_equal(out$expected_runs, 40.9, tolerance = 1e-12)
  n_p <- 42
  n_n <- 38
  n <- n_p + n_n
  variance_expected <- 2 * n_p * n_n * (2 * n_p * n_n - n) / (n^2 * (n - 1))
  expect_equal(out$variance_runs, variance_expected, tolerance = 1e-12)
  expect_published(out$variance_runs, 19.65, digits = 2L)

  z <- (44 - out$expected_runs) / sqrt(variance_expected)
  expect_equal(z, 0.6993, tolerance = 2e-4)
  expect_published(pnorm(z, lower.tail = FALSE), 0.2422, digits = 4L)
})

test_that("compute_runs_statistics returns the retained signs it counted", {
  out <- compute_runs_statistics(c(1, 0, -1, -1, 2), "two.sided")
  expect_identical(out$retained, c(TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_identical(out$signs, c(1L, -1L, -1L, 1L))
  expect_identical(out$zero, 1L)
  expect_identical(out$nonzero, 4L)
  expect_identical(out$runs, 3L)
})

test_that("runs_test can order residuals by fitted values", {
  out <- runs_test(fit_geewa_pois_exch, order_by = "fitted")
  residual_values <- residuals(fit_geewa_pois_exch, type = "working")
  natural <- order(fit_geewa_pois_exch$id, fit_geewa_pois_exch$repeated, seq_along(residual_values))
  ord <- natural[order(fit_geewa_pois_exch$fitted.values[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "fitted values")
  expect_identical(out$runs, as.integer(count_runs(signs)))
  expect_identical(out$positive, as.integer(sum(signs > 0)))
  expect_identical(out$negative, as.integer(sum(signs < 0)))
})

test_that("runs_test can order residuals by a model-matrix covariate", {
  out <- runs_test(fit_geewa_pois_exch, order_by = "lnage")
  residual_values <- residuals(fit_geewa_pois_exch, type = "working")
  natural <- order(fit_geewa_pois_exch$id, fit_geewa_pois_exch$repeated, seq_along(residual_values))
  key <- fit_geewa_pois_exch$x[, "lnage"]
  ord <- natural[order(key[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "model-matrix column 'lnage'")
  expect_identical(out$runs, as.integer(count_runs(signs)))
})

test_that("order_by can name a covariate outside the model matrix", {
  fit <- geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
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
  key <- rep(c(2, 1, 3), length.out = fit_geewa_pois_exch$obs_no)
  out <- runs_test(fit_geewa_pois_exch, order_by = key)
  residual_values <- residuals(fit_geewa_pois_exch, type = "working")
  natural <- order(fit_geewa_pois_exch$id, fit_geewa_pois_exch$repeated, seq_along(residual_values))
  ord <- natural[order(key[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "supplied ordering vector")
  expect_identical(out$runs, as.integer(count_runs(signs)))
})

test_that("zero residuals are omitted from the sign sequence", {
  fit <- fit_geewa_pois_exch
  fit$residuals <- rep(c(1, 0, -1, 0), length.out = fit$obs_no)
  out <- runs_test(fit)

  expect_identical(out$zero, as.integer(sum(fit$residuals == 0)))
  expect_identical(out$nonzero, fit$obs_no - out$zero)
  expect_identical(out$positive + out$negative, out$nonzero)
})

test_that("the result carries the tested sign sequence", {
  out <- runs_test(fit_geewa_pois_exch)

  expect_s3_class(out, "geer_runs_test")
  expect_s3_class(out, "htest")
  expect_s3_class(out, "geer_htest")
  expect_true(out$natural.order)
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
    as.numeric(sort(unique(fit_geewa_pois_exch$id)))
  )
})

test_that("natural.order reports the ordering obtained, not requested", {
  tied <- runs_test(fit_geewa_pois_exch, order_by = rep(1, fit_geewa_pois_exch$obs_no))
  natural <- runs_test(fit_geewa_pois_exch)

  expect_true(tied$natural.order)
  expect_identical(tied$signs, natural$signs)
  expect_identical(tied$runs, natural$runs)
  expect_identical(tied$order_by, "supplied ordering vector")

  by_fitted <- runs_test(fit_geewa_pois_exch, order_by = "fitted")
  expect_false(by_fitted$natural.order)
})

test_that("runs_test validates the ordering", {
  expect_error(
    runs_test(fit_geewa_pois_exch, order_by = "not-a-variable"),
    "unknown 'order_by'"
  )
  expect_error(
    runs_test(fit_geewa_pois_exch, order_by = seq_len(fit_geewa_pois_exch$obs_no - 1L)),
    "one value per fitted observation"
  )

  bad_order <- seq_len(fit_geewa_pois_exch$obs_no)
  bad_order[1] <- NA_real_
  expect_error(
    runs_test(fit_geewa_pois_exch, order_by = bad_order),
    "finite and non-missing"
  )

})

test_that("large sign counts keep the normal approximation", {
  out <- runs_test(fit_geewa_pois_exch)
  expect_equal(
    out$p.value,
    2 * pnorm(-abs(unname(out$statistic))),
    tolerance = 1e-12
  )
})

test_that("runs_test warns when the normal approximation is based on small sign counts", {
  fit <- fit_geewa_pois_exch
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
