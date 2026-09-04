testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch


count_runs <- function(signs) {
  signs <- signs[signs != 0]
  1L + sum(signs[-1L] != signs[-length(signs)])
}


test_that("runs-test calculation matches the Chang sign-sequence example", {
  signs <- c(1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1)
  out <- geer:::compute_runs_statistics(signs)

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
    function(x) geer:::compute_runs_statistics(x)
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
    "should be one of"
  )
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


test_that("the exact null distribution matches full enumeration", {
  ## Every arrangement of n_p positive and n_n negative signs is enumerated and
  ## the resulting tail probabilities compared with the closed-form expression.
  enumerate_runs <- function(positive_no, negative_no) {
    n <- positive_no + negative_no
    positions <- utils::combn(n, positive_no)
    apply(positions, 2L, function(pos) {
      s <- rep(-1, n)
      s[pos] <- 1
      1L + sum(s[-1L] != s[-n])
    })
  }

  for (counts in list(c(4L, 3L), c(5L, 5L), c(2L, 6L), c(7L, 3L))) {
    np <- counts[[1L]]
    nn <- counts[[2L]]
    all_runs <- enumerate_runs(np, nn)
    for (t in sort(unique(all_runs))) {
      expect_equal(
        geer:::compute_runs_exact_p(np, nn, t, "less"),
        mean(all_runs <= t),
        tolerance = 1e-12
      )
      expect_equal(
        geer:::compute_runs_exact_p(np, nn, t, "greater"),
        mean(all_runs >= t),
        tolerance = 1e-12
      )
      expect_equal(
        geer:::compute_runs_exact_p(np, nn, t, "two.sided"),
        min(1, 2 * min(mean(all_runs <= t), mean(all_runs >= t))),
        tolerance = 1e-12
      )
    }
  }
})


test_that("the exact distribution matches the continuity-corrected normal", {
  ## runs_test applies no continuity correction, following Chang (2000), so the
  ## exact tail probabilities agree with the *corrected* normal approximation
  ## and differ from the uncorrected one by roughly that correction even at
  ## large sign counts.
  positive_no <- 500L
  negative_no <- 500L
  n <- positive_no + negative_no
  expected_runs <- 2 * positive_no * negative_no / n + 1
  variance_runs <- 2 * positive_no * negative_no *
    (2 * positive_no * negative_no - n) / (n^2 * (n - 1))

  for (t in c(510L, 520L, 530L)) {
    corrected <- stats::pnorm(
      (t - 0.5 - expected_runs) / sqrt(variance_runs),
      lower.tail = FALSE
    )
    expect_equal(
      geer:::compute_runs_exact_p(positive_no, negative_no, t, "greater"),
      corrected,
      tolerance = 1e-3
    )
  }

  ## The uncorrected approximation, which is what exact = FALSE reports, is
  ## visibly different.
  uncorrected <- stats::pnorm(
    (520 - expected_runs) / sqrt(variance_runs),
    lower.tail = FALSE
  )
  expect_gt(
    abs(
      geer:::compute_runs_exact_p(positive_no, negative_no, 520L, "greater") -
        uncorrected
    ),
    0.005
  )
})


test_that("small sign counts use the exact distribution instead of warning", {
  fit <- count_fit
  fit$residuals <- c(
    rep(1, 10),
    rep(-1, 10),
    rep(0, fit$obs_no - 20)
  )

  expect_silent(out <- runs_test(fit))
  expect_true(out$exact)
  expect_true(is.na(out$nperm))
  expect_match(out$method, "exact null distribution")
  expect_identical(out$positive, 10L)
  expect_identical(out$negative, 10L)
  expect_equal(
    out$p.value,
    geer:::compute_runs_exact_p(10L, 10L, out$runs, "two.sided"),
    tolerance = 1e-12
  )

  ## The normal approximation, and its warning, remain available on request.
  expect_warning(
    approximate <- runs_test(fit, exact = FALSE),
    "normal approximation may be unreliable"
  )
  expect_false(approximate$exact)
  expect_equal(
    approximate$p.value,
    2 * stats::pnorm(-abs(unname(approximate$statistic))),
    tolerance = 1e-12
  )
})


test_that("exact and nperm cannot both be requested", {
  expect_error(
    runs_test(count_fit, nperm = 99, exact = TRUE),
    "must not both be requested"
  )
  expect_error(runs_test(count_fit, exact = "yes"), "single non-missing logical")
  expect_error(runs_test(count_fit, exact = NA), "single non-missing logical")
})


test_that("large sign counts keep the normal approximation by default", {
  out <- runs_test(count_fit)
  expect_false(out$exact)
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
    out <- runs_test(fit, exact = FALSE),
    "normal approximation may be unreliable"
  )
  expect_identical(out$positive, 10L)
  expect_identical(out$negative, 10L)
})


test_that("permutation p-values are Monte Carlo tests on the same statistic", {
  set.seed(101)
  out <- runs_test(count_fit, nperm = 199)
  analytic <- runs_test(count_fit)

  expect_identical(out$nperm, 199L)
  expect_true(is.na(analytic$nperm))
  expect_identical(out$runs, analytic$runs)
  ## The statistic follows the null actually used, so it is standardized by
  ## the permutation moments rather than by E(T) and V(T).
  expect_equal(
    unname(out$statistic),
    (out$runs - out$permuted_mean) / out$permuted_sd,
    tolerance = 1e-12
  )
  expect_identical(unname(out$null.value), out$permuted_mean)
  expect_identical(out$expected_runs, analytic$expected_runs)
  expect_identical(out$variance_runs, analytic$variance_runs)
  expect_true(is.na(analytic$permuted_mean))
  expect_true(is.na(analytic$permuted_sd))
  expect_match(out$method, "199 within-cluster permutations")
  expect_gt(out$p.value, 0)
  expect_lte(out$p.value, 1)
  ## Attainable p-values are multiples of 1 / (nperm + 1), doubled and capped
  ## for the two-sided case.
  expect_equal(
    min(1, round(out$p.value * 200 / 2) * 2 / 200),
    out$p.value,
    tolerance = 1e-12
  )
})


test_that("permutation p-values are reproducible under set.seed", {
  set.seed(202)
  first <- runs_test(count_fit, nperm = 99)
  set.seed(202)
  second <- runs_test(count_fit, nperm = 99)
  set.seed(303)
  third <- runs_test(count_fit, nperm = 99)

  expect_identical(first$p.value, second$p.value)
  expect_false(is.null(third$p.value))
})


test_that("within-cluster permutation does not reject cluster lopsidedness", {
  ## Every cluster is internally constant in sign, so the analytic test sees
  ## far too few runs while no within-cluster permutation can change T.
  cluster <- rep(seq_len(20L), each = 5L)
  signs <- rep(rep(c(1, -1), length.out = 20L), each = 5L)

  analytic <- geer:::compute_runs_statistics(signs, "less")
  expect_identical(analytic$runs, 20L)
  expect_lt(analytic$statistic, -4)
  expect_lt(analytic$p_value, 1e-4)

  set.seed(404)
  expect_warning(
    permuted <- geer:::compute_runs_permutation(
      residual_values = signs,
      cluster = cluster,
      runs = analytic$runs,
      alternative = "less",
      nperm = 99
    ),
    "same number of runs"
  )
  expect_true(permuted$degenerate)
  expect_identical(permuted$p_value, 1)
  expect_identical(permuted$sd, 0)
  expect_identical(permuted$mean, as.numeric(analytic$runs))
})


test_that("a degenerate permutation null gives a missing statistic", {
  ## Standardizing by a zero permutation standard deviation is undefined, so
  ## no statistic is reported rather than the misleading analytic one.
  fit <- count_fit
  cluster_sizes <- as.integer(table(fit$id))
  fit$residuals <- unlist(
    Map(
      function(size, s) rep(s, size),
      cluster_sizes,
      rep(c(1, -1), length.out = length(cluster_sizes))
    ),
    use.names = FALSE
  )

  set.seed(707)
  expect_warning(out <- runs_test(fit, nperm = 99), "same number of runs")
  expect_true(is.na(unname(out$statistic)))
  expect_identical(out$p.value, 1)
  expect_identical(out$permuted_sd, 0)
})


test_that("within-cluster permutation still detects a within-cluster pattern", {
  ## A common sign pattern inside every cluster is destroyed by permutation.
  cluster <- rep(seq_len(20L), each = 5L)
  signs <- rep(c(-1, -1, -1, 1, 1), times = 20L)
  runs_observed <- geer:::compute_runs_statistics(signs, "less")$runs

  set.seed(505)
  permuted <- geer:::compute_runs_permutation(
    residual_values = signs,
    cluster = cluster,
    runs = runs_observed,
    alternative = "less",
    nperm = 199
  )
  expect_false(permuted$degenerate)
  expect_lt(permuted$p_value, 0.05)
})


test_that("permutation ignores zero residuals and validates nperm", {
  fit <- count_fit
  fit$residuals <- rep(c(1, 0, -1, 0), length.out = fit$obs_no)
  set.seed(606)
  out <- runs_test(fit, nperm = 99)
  expect_identical(out$nonzero, fit$obs_no - out$zero)

  expect_error(runs_test(count_fit, nperm = 0), "positive whole number")
  expect_error(runs_test(count_fit, nperm = 2.5), "positive whole number")
  expect_error(runs_test(count_fit, nperm = c(10, 20)), "positive whole number")
  expect_error(runs_test(count_fit, nperm = NA), "positive whole number")
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


test_that("plot.geer_runs_test accepts point arguments", {
  out <- runs_test(count_fit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot(out, col = "black", pch = 20L, cex = 0.5))
  expect_silent(plot(out, run_colors = c("red", "blue", "green")))
  expect_error(plot(out, run_colors = character(0)), "at least one colour")
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
