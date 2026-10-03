expect_test_result <- function(x) {
  expect_type(x, "list")
  expect_true(all(c("test_stat", "test_df", "test_p") %in% names(x)))
  expect_true(is.numeric(x$test_stat))
  expect_true(is.numeric(x$test_df))
  expect_true(is.numeric(x$test_p))
  expect_length(x$test_stat, 1L)
  expect_length(x$test_df, 1L)
  expect_length(x$test_p, 1L)
  expect_true(is.finite(x$test_stat))
  expect_true(is.finite(x$test_df))
  expect_true(is.finite(x$test_p))
  expect_gte(x$test_stat, 0)
  expect_gt(x$test_df, 0)
  expect_gte(x$test_p, 0)
  expect_lte(x$test_p, 1)
}

jackknife_delete_estimates_by_refit <- function(fit, data, refit) {
  # Brute-force leave-one-cluster-out estimates; refit() receives the data with
  # one cluster removed and must return a fitted model.
  ids <- unique(data$id)
  out <- t(vapply(
    ids,
    function(id_out) coef(refit(data[data$id != id_out, ])),
    numeric(length(coef(fit)))
  ))
  dimnames(out) <- list(as.character(ids), names(coef(fit)))
  out
}

expect_jackknife_vcov <- function(fit, delete_estimates) {
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  expected <- ((nrow(delete_estimates) - 1L) / nrow(delete_estimates)) *
    crossprod(centered)
  dimnames(expected) <- list(names(coef(fit)), names(coef(fit)))
  observed <- vcov(fit, cov_type = "jackknife")
  expect_equal(compute_jackknife_delete_estimates(fit), delete_estimates,
               tolerance = 1e-7)
  expect_equal(observed, expected, tolerance = 1e-7)
  expect_equal(observed, t(observed), tolerance = 1e-12)
}
