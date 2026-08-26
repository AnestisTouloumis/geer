testthat::local_edition(3)


fit_binary_indep <- geewa_binary(
  formula = ecg ~ period * treatment,
  id = id,
  data = test_data$cerebrovascular,
  link = "logit",
  orstr = "independence",
  method = "gee"
)

fit_binary_exch <- update(
  fit_binary_indep,
  orstr = "exchangeable"
)

cerebrovascular_small <- test_data$cerebrovascular[
  seq_len(nrow(test_data$cerebrovascular) - 1L),
  ,
  drop = FALSE
]

fit_binary_small <- geewa_binary(
  formula = ecg ~ period * treatment,
  id = id,
  data = cerebrovascular_small,
  link = "logit",
  orstr = "independence",
  method = "gee"
)


test_that("geecriteria returns the expected structure for single and multiple models", {
  out_single <- geecriteria(fit_geewa_pois_exch)
  expect_geecriteria_table(out_single, n_rows = 1L)
  expect_equal(out_single$Parameters, length(coef(fit_geewa_pois_exch)))
  out_multi <- geecriteria(
    fit_binary_indep,
    fit_binary_exch,
    cov_type = "robust",
    digits = 3
  )
  expect_geecriteria_table(
    out_multi,
    n_rows = 2L,
    row_names = c("fit_binary_indep", "fit_binary_exch")
  )
})


test_that("geecriteria warns when models differ in number of observations", {
  expect_warning(
    geecriteria(fit_binary_indep, fit_binary_small),
    regexp = "same number of observations",
    ignore.case = TRUE
  )
})


test_that("geecriteria rejects invalid inputs", {
  lm_fit <- lm(seizures ~ treatment, data = test_data$epilepsy)
  expect_error(
    geecriteria(1),
    regexp = "geer",
    ignore.case = TRUE
  )
  expect_error(
    geecriteria(fit_geewa_pois_exch, lm_fit),
    regexp = "geer|Only 'geer' objects are supported",
    ignore.case = TRUE
  )
  expect_error(
    geecriteria(fit_binary_indep, cov_type = "not-a-type"),
    regexp = "cov_type|arg",
    ignore.case = TRUE
  )
  expect_error(
    geecriteria(fit_binary_indep, digits = 1.2),
    regexp = "digits",
    ignore.case = TRUE
  )
})


test_that("geecriteria works for representative model and cov_type variants", {
  expect_geecriteria_table(geecriteria(fit_geewa_bin_exch), n_rows = 1L)
  for (cov_type in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    out <- geecriteria(fit_geewa_pois_exch, cov_type = cov_type)
    expect_geecriteria_table(out, n_rows = 1L)
  }
})


test_that("geecriteria defaults to the classical robust covariance", {
  out_default <- geecriteria(fit_geewa_pois_exch, digits = 15)
  out_robust <- geecriteria(
    fit_geewa_pois_exch,
    cov_type = "robust",
    digits = 15
  )
  expect_equal(out_default, out_robust, tolerance = 1e-12)
})


test_that("QICHH agrees with QIC under working independence for ordinary GEE", {
  out <- geecriteria(fit_geewa_pois_indep, digits = 15)
  expect_equal(out$QICHH, out$QIC, tolerance = 1e-6)
})


test_that("QICHH uses independence estimates in its quasi-likelihood and penalty", {
  object <- fit_geewa_pois_exch
  quantities <- geer:::compute_independence_gee_quantities(object)
  covariance <- vcov(object, cov_type = "robust")
  penalty <- sum(quantities$naive_inverse * covariance)
  expected <- 2 * (penalty - quantities$quasi_loglikelihood)
  observed <- geecriteria(object, digits = 15)$QICHH
  expect_equal(observed, expected, tolerance = 1e-10)
})


test_that("QICC follows the Hardin-Hilbe finite-cluster correction", {
  object <- fit_geewa_pois_exch
  out <- geecriteria(object, digits = 15)
  p <- length(coef(object))
  m <- geer:::compute_n_estimated_association_parameters(object)
  n_clusters <- object$clusters_no
  correction <- 2 * (p + m) * (p + m + 1) /
    (n_clusters - p - m - 1)
  expected <- out$QIC - correction
  expect_equal(out$QICC, expected, tolerance = 1e-10)
})


test_that("QICC counts only estimated working-association parameters", {
  expect_equal(
    geer:::compute_n_estimated_association_parameters(list(
      association_structure = "independence",
      alpha = 0
    )),
    0L
  )
  expect_equal(
    geer:::compute_n_estimated_association_parameters(list(
      association_structure = "fixed",
      alpha = c(0.2, 0.1, 0.3)
    )),
    0L
  )
  expect_equal(
    geer:::compute_n_estimated_association_parameters(list(
      association_structure = "exchangeable",
      alpha = 0.2
    )),
    1L
  )
})


test_that("QICC accounts for working-association dimensionality", {
  qic <- 100
  p <- 4
  n_clusters <- 50
  qicc_ind <- geer:::compute_qicc(
    qic,
    p = p,
    association_params_no = 0,
    clusters_no = n_clusters
  )
  qicc_exch <- geer:::compute_qicc(
    qic,
    p = p,
    association_params_no = 1,
    clusters_no = n_clusters
  )
  expect_equal(qicc_ind, 100 - 2 * 4 * 5 / (50 - 4 - 1))
  expect_equal(qicc_exch, 100 - 2 * 5 * 6 / (50 - 4 - 1 - 1))
  expect_false(isTRUE(all.equal(qicc_ind, qicc_exch)))
})


test_that("QICC is undefined when the finite-cluster denominator is nonpositive", {
  expect_true(is.na(geer:::compute_qicc(
    qic = 100,
    p = 4,
    association_params_no = 1,
    clusters_no = 6
  )))
  expect_true(is.na(geer:::compute_qicc(
    qic = 100,
    p = 4,
    association_params_no = 1,
    clusters_no = 5
  )))
})


test_that("EQIC uses the adjusted extended quasi-likelihood", {
  object <- fit_geewa_pois_exch
  k <- 1 / 6
  mu <- object$fitted.values
  y <- object$y
  weights <- object$prior.weights
  z <- y + k
  mu_adjusted <- mu + k
  deviance_contributions <- 2 * (
    z * log(z / mu_adjusted) - (z - mu_adjusted)
  )
  deviance <- sum(weights * deviance_contributions)
  phi <- deviance / object$obs_no
  log_variance_term <- sum(
    weights * log(2 * pi * phi * mu_adjusted)
  )
  independence_inverse <- geer:::compute_independence_naive_inverse(
    object,
    phi = phi
  )
  covariance <- vcov(object, cov_type = "robust")
  penalty <- sum(independence_inverse * covariance)
  expected <- deviance / phi + log_variance_term + 2 * penalty
  observed <- geecriteria(object, digits = 15)$EQIC
  expect_equal(observed, expected, tolerance = 1e-10)
})


test_that("EQIC adjustment is finite at discrete-response boundaries", {
  k <- 1 / 6
  poisson_dev <- geer:::compute_eqic_adjusted_deviance(
    y = c(0, 1, 3),
    mu = c(0.2, 1.1, 2.5),
    family_name = "poisson",
    k = k
  )
  binomial_dev <- geer:::compute_eqic_adjusted_deviance(
    y = c(0, 1),
    mu = c(0.1, 0.9),
    family_name = "binomial",
    k = k
  )
  binomial_var <- geer:::compute_eqic_adjusted_variance(
    mu = c(0, 1),
    family_name = "binomial",
    k = k
  )
  expect_true(all(is.finite(poisson_dev)))
  expect_true(all(poisson_dev >= 0))
  expect_true(all(is.finite(binomial_dev)))
  expect_true(all(binomial_dev >= 0))
  expect_true(all(is.finite(binomial_var)))
  expect_true(all(binomial_var > 0))

  for (family_name in c("gaussian", "poisson", "Gamma", "inverse.gaussian")) {
    mu <- c(0.5, 1.5, 2.5)
    deviance_at_fit <- geer:::compute_eqic_adjusted_deviance(
      y = mu,
      mu = mu,
      family_name = family_name,
      k = k
    )
    expect_equal(deviance_at_fit, rep(0, length(mu)), tolerance = 1e-12)
  }
  binomial_mu <- c(0.1, 0.5, 0.9)
  expect_equal(
    geer:::compute_eqic_adjusted_deviance(
      y = binomial_mu,
      mu = binomial_mu,
      family_name = "binomial",
      k = k
    ),
    rep(0, length(binomial_mu)),
    tolerance = 1e-12
  )
})



test_that("AGPC and SGPC use the conventional penalized Gaussian pseudo-likelihood", {
  object <- fit_geewa_pois_exch
  out <- geecriteria(object, digits = 15)
  p <- length(coef(object))
  q <- geer:::compute_n_estimated_association_parameters(object)
  gaussian_deviance <- object$obs_no * log(2 * pi) - 2 * out$GPC
  expect_equal(
    out$AGPC,
    gaussian_deviance + 2 * (p + q),
    tolerance = 1e-10
  )
  expect_equal(
    out$SGPC,
    gaussian_deviance + log(object$clusters_no) * (p + q),
    tolerance = 1e-10
  )
})


test_that("penalized Gaussian pseudo-likelihood does not count fixed association parameters", {
  object <- fit_geewa_pois_exch
  gpc <- -25
  p <- length(coef(object))
  base <- geer:::compute_gaussian_pseudolikelihood_criteria(
    object = object,
    gpc = gpc,
    p = p,
    association_params_no = 0
  )
  with_one_association_parameter <- geer:::compute_gaussian_pseudolikelihood_criteria(
    object = object,
    gpc = gpc,
    p = p,
    association_params_no = 1
  )
  expect_equal(with_one_association_parameter$AGPC - base$AGPC, 2)
  expect_equal(
    with_one_association_parameter$SGPC - base$SGPC,
    log(object$clusters_no)
  )
})


test_that("GHYC and PAC agree with direct covariance calculations under independence", {
  object <- fit_binary_indep
  repeated_max <- max(as.integer(object$repeated))
  empirical_sum <- matrix(0, repeated_max, repeated_max)
  working_sum <- matrix(0, repeated_max, repeated_max)
  pair_counts <- matrix(0, repeated_max, repeated_max)
  residuals <- object$y - object$fitted.values

  for (indices in split(seq_along(object$id), object$id)) {
    repeated <- as.integer(object$repeated[indices])
    mu <- object$fitted.values[indices]
    weights <- object$prior.weights[indices]
    empirical_sum[repeated, repeated] <-
      empirical_sum[repeated, repeated, drop = FALSE] +
      tcrossprod(residuals[indices])
    working_sum[repeated, repeated] <-
      working_sum[repeated, repeated, drop = FALSE] +
      diag(mu * (1 - mu) / weights, nrow = length(indices))
    pair_counts[repeated, repeated] <-
      pair_counts[repeated, repeated, drop = FALSE] + 1
  }

  empirical_mean <- empirical_sum / pair_counts
  working_mean <- working_sum / pair_counts
  discrepancy <- empirical_mean %*% solve(working_mean) - diag(repeated_max)
  expected_ghyc <- sum(diag(discrepancy %*% discrepancy))
  expected_pac <- abs(det(empirical_mean) / det(working_mean) - 1)
  out <- geecriteria(object, digits = 15)

  expect_equal(out$GHYC, expected_ghyc, tolerance = 1e-10)
  expect_equal(out$PAC, expected_pac, tolerance = 1e-10)
})


test_that("GHYC and PAC use pairwise-available repeated positions", {
  out <- geecriteria(fit_binary_small, digits = 15)
  expect_true(is.finite(out$GHYC))
  expect_true(is.finite(out$PAC))
})


test_that("working covariance helper supports odds-ratio GEE", {
  object <- fit_binary_exch
  indices <- split(seq_along(object$id), object$id)[[1L]]
  covariance <- geer:::compute_working_covariance_for_criteria(object, indices)
  expect_equal(nrow(covariance), length(indices))
  expect_equal(ncol(covariance), length(indices))
  expect_equal(covariance, t(covariance), tolerance = 1e-12)
  expect_true(all(is.finite(covariance)))
  expect_true(all(diag(covariance) > 0))
})


test_that("working covariance helper supports correlation GEE", {
  object <- fit_geewa_pois_exch
  indices <- split(seq_along(object$id), object$id)[[1L]]
  repeated <- as.integer(object$repeated[indices])
  repeated_max <- max(as.integer(object$repeated))
  correlation <- geer:::get_correlation_matrix(
    object$association_structure,
    object$alpha,
    repeated_max
  )
  mu <- object$fitted.values[indices]
  weights <- object$prior.weights[indices]
  marginal_sd <- sqrt(object$phi * object$family$variance(mu) / weights)
  expected <- correlation[repeated, repeated, drop = FALSE] *
    tcrossprod(marginal_sd)
  observed <- geer:::compute_working_covariance_for_criteria(object, indices)
  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("odds-ratio pair indexing matches upper-triangular ordering", {
  expect_equal(
    vapply(
      list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4)),
      function(pair) geer:::compute_upper_triangular_pair_index(
        pair[[1L]], pair[[2L]], 4L
      ),
      integer(1)
    ),
    seq_len(6L)
  )
})
