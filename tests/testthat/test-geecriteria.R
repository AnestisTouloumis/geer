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
  seq_len(nrow(test_data$cerebrovascular) - 1L), ,
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

test_that("geecriteria stops when models differ in number of observations", {
  expect_error(
    geecriteria(fit_binary_indep, fit_binary_small),
    "models must be fitted to the same observations",
    fixed = TRUE
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
  finite_columns <- setdiff(geecriteria_cols(), undefined_geecriteria_cols)
  out <- geecriteria(fit_geewa_bin_exch)
  expect_true(all(is.finite(as.matrix(out[finite_columns]))))
  for (cov_type in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    out <- geecriteria(fit_geewa_pois_exch, cov_type = cov_type)
    expect_true(
      all(is.finite(as.matrix(out[finite_columns]))),
      label = cov_type
    )
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

test_that("QICHH agrees with QIC for binomial geewa fits with estimated phi", {
  fit <- geewa(
    ecg ~ period * treatment,
    family = binomial("logit"),
    data = test_data$cerebrovascular,
    id = id,
    corstr = "independence"
  )
  expect_false(isTRUE(all.equal(fit$phi, 1)))
  out <- geecriteria(fit, criteria = c("QIC", "QICHH"), digits = 15)
  expect_equal(out$QICHH, out$QIC, tolerance = 1e-6)
  quantities <- compute_independence_fit(fit)
  expect_equal(quantities$phi, fit$phi, tolerance = 1e-6)

  fit_fixed <- update(fit, phi_fixed = TRUE, phi_value = 1)
  out_fixed <- geecriteria(fit_fixed, criteria = c("QIC", "QICHH"), digits = 15)
  expect_equal(out_fixed$QICHH, out_fixed$QIC, tolerance = 1e-6)
  expect_identical(compute_independence_fit(fit_fixed)$phi, 1)
})

test_that("QICHH uses independence estimates in its quasi-likelihood and penalty", {
  object <- fit_geewa_pois_exch
  quantities <- compute_independence_fit(object)
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
  m <- get_association_count(object)
  n_clusters <- object$clusters_no
  correction <- 2 * (p + m) * (p + m + 1) /
    (n_clusters - p - m - 1)
  expected <- out$QIC - correction
  expect_equal(out$QICC, expected, tolerance = 1e-10)
})

test_that("QICC under independence uses the data's cluster count and design width", {
  fit_ind <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence"
  )
  out <- geecriteria(fit_ind, digits = 15)
  n_clusters <- length(unique(test_data$epilepsy$id))
  p <- 4L
  correction <- 2 * p * (p + 1) / (n_clusters - p - 1)
  expect_equal(out$QICC, out$QIC - correction, tolerance = 1e-10)
})

test_that("QICC counts only estimated working-association parameters", {
  expect_equal(
    get_association_count(list(
      association_structure = "independence",
      alpha = 0
    )),
    0L
  )
  expect_equal(
    get_association_count(list(
      association_structure = "fixed",
      alpha = c(0.2, 0.1, 0.3)
    )),
    0L
  )
  expect_equal(
    get_association_count(list(
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
  qicc_ind <- compute_qicc(
    qic,
    p = p,
    association_params_no = 0,
    clusters_no = n_clusters
  )
  qicc_exch <- compute_qicc(
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
  expect_true(is.na(compute_qicc(
    qic = 100,
    p = 4,
    association_params_no = 1,
    clusters_no = 6
  )))
  expect_true(is.na(compute_qicc(
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
  independence_inverse <- compute_independence_info(
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
  poisson_dev <- compute_eqic_adjusted_deviance(
    y = c(0, 1, 3),
    mu = c(0.2, 1.1, 2.5),
    family_name = "poisson",
    k = k
  )
  binomial_dev <- compute_eqic_adjusted_deviance(
    y = c(0, 1),
    mu = c(0.1, 0.9),
    family_name = "binomial",
    k = k
  )
  binomial_var <- compute_eqic_adjusted_variance(
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
    deviance_at_fit <- compute_eqic_adjusted_deviance(
      y = mu,
      mu = mu,
      family_name = family_name,
      k = k
    )
    expect_equal(deviance_at_fit, rep(0, length(mu)), tolerance = 1e-12)
  }
  binomial_mu <- c(0.1, 0.5, 0.9)
  expect_equal(
    compute_eqic_adjusted_deviance(
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
  q <- get_association_count(object)
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
  base <- compute_gpc_penalized(
    object = object,
    gpc = gpc,
    p = p,
    association_params_no = 0
  )
  with_one_association_parameter <- compute_gpc_penalized(
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

test_that("all complexity penalties count only estimated association parameters", {
  exch <- update(fit_geewa_pois_exch, corstr = "exchangeable")
  fixed_fit <- update(
    fit_geewa_pois_exch,
    corstr = "fixed",
    alpha_vector = rep(exch$alpha, choose(4L, 2L))
  )
  expect_identical(
    get_association_count(fixed_fit),
    0L
  )
  expect_identical(
    get_association_count(exch),
    1L
  )

  p <- length(coef(fixed_fit))
  out <- geecriteria(fixed_fit, digits = 15)
  gaussian_deviance <- fixed_fit$obs_no * log(2 * pi) - 2 * out$GPC
  expect_equal(out$AGPC, gaussian_deviance + 2 * p, tolerance = 1e-10)
  expect_equal(
    out$SGPC,
    gaussian_deviance + log(fixed_fit$clusters_no) * p,
    tolerance = 1e-10
  )
  expect_equal(
    out$QICC,
    out$QIC - 2 * p * (p + 1) / (fixed_fit$clusters_no - p - 1),
    tolerance = 1e-10
  )

  out_exch <- geecriteria(exch, digits = 15)
  expect_equal(
    out$GESSC / out_exch$GESSC,
    (exch$obs_no - p - 1) / (fixed_fit$obs_no - p),
    tolerance = 1e-8
  )
})

test_that("GHYC and PAC agree with direct covariance calculations under independence", {
  object <- fit_binary_indep
  repeated_max <- max(as.integer(object$repeated))
  cluster_indices <- split(seq_along(object$id), object$id)
  empirical_sum <- matrix(0, repeated_max, repeated_max)
  working_sum <- matrix(0, repeated_max, repeated_max)
  residuals <- object$y - object$fitted.values

  for (indices in cluster_indices) {
    repeated <- as.integer(object$repeated[indices])
    mu <- object$fitted.values[indices]
    weights <- object$prior.weights[indices]
    empirical_sum[repeated, repeated] <-
      empirical_sum[repeated, repeated, drop = FALSE] +
      tcrossprod(residuals[indices])
    working_sum[repeated, repeated] <-
      working_sum[repeated, repeated, drop = FALSE] +
      diag(mu * (1 - mu) / weights, nrow = length(indices))
  }

  empirical_mean <- empirical_sum / length(cluster_indices)
  working_mean <- working_sum / length(cluster_indices)
  discrepancy <- empirical_mean %*% solve(working_mean) - diag(repeated_max)
  expected_ghyc <- sum(diag(discrepancy %*% discrepancy))
  expected_pac <- abs(det(empirical_mean) / det(working_mean) - 1)
  out <- geecriteria(object, digits = 15)

  expect_equal(out$GHYC, expected_ghyc, tolerance = 1e-10)
  expect_equal(out$PAC, expected_pac, tolerance = 1e-10)
})

test_that("an unavailable criterion is reported as NA, not an error", {
  fit <- fit_geewa_pois_exch
  broken <- fit
  broken$naive_covariance <- NULL

  expect_identical(
    compute_criterion_or_na(
      signal_criterion_unavailable("nope")
    ),
    NA_real_
  )
  expect_error(
    compute_criterion_or_na(stop("a real bug")),
    "a real bug"
  )
})

test_that("the quasi-likelihood clamp keeps boundary fits finite", {
  eps <- sqrt(.Machine$double.eps)
  y <- c(0, 1, 2, 3)
  weights <- rep(1, 4L)

  clamped <- compute_quasi_loglikelihood(
    y = y,
    mu = c(0, 1, 2, 3),
    weights = weights,
    family_name = "poisson",
    phi = 1
  )
  expect_true(is.finite(clamped))
  expect_equal(
    clamped,
    sum(weights * (y * log(pmax(c(0, 1, 2, 3), eps)) - pmax(c(0, 1, 2, 3), eps))),
    tolerance = 1e-10
  )

  binomial_value <- compute_quasi_loglikelihood(
    y = c(0, 1),
    mu = c(0, 1),
    weights = c(1, 1),
    family_name = "binomial",
    phi = 1
  )
  expect_true(is.finite(binomial_value))
})

test_that("criteria selects and orders the reported columns", {
  full <- geecriteria(fit_geewa_pois_exch, digits = 15)

  selected <- geecriteria(
    fit_geewa_pois_exch,
    criteria = c("CIC", "QIC"),
    digits = 15
  )
  expect_named(selected, c("CIC", "QIC", "Parameters"))
  expect_equal(selected$CIC, full$CIC, tolerance = 1e-12)
  expect_equal(selected$QIC, full$QIC, tolerance = 1e-12)
  expect_identical(selected$Parameters, full$Parameters)

  one <- geecriteria(fit_geewa_pois_exch, criteria = "PT", digits = 15)
  expect_s3_class(one, "data.frame")
  expect_named(one, c("PT", "Parameters"))

  expect_named(
    geecriteria(fit_geewa_pois_exch, criteria = c("cic", "CIC", "qic")),
    c("CIC", "QIC", "Parameters")
  )

  expect_named(
    geecriteria(fit_geewa_pois_exch, criteria = "all"),
    names(full)
  )
})

test_that("criteria selection works with several models", {
  fit_ind <- update(fit_geewa_pois_exch, corstr = "independence")
  out <- geecriteria(
    fit_geewa_pois_exch,
    fit_ind,
    criteria = c("QIC", "CIC")
  )
  expect_named(out, c("QIC", "CIC", "Parameters"))
  expect_equal(nrow(out), 2L)
  expect_identical(
    rownames(out),
    c("fit_geewa_pois_exch", "fit_ind")
  )
})

test_that("criteria rejects invalid selections", {
  expect_error(
    geecriteria(fit_geewa_pois_exch, criteria = "QICX"),
    "unknown entries in 'criteria'"
  )
  expect_error(
    geecriteria(fit_geewa_pois_exch, criteria = c("all", "QIC")),
    "not both"
  )
  expect_error(
    geecriteria(fit_geewa_pois_exch, criteria = character(0)),
    "non-missing character vector"
  )
  expect_error(
    geecriteria(fit_geewa_pois_exch, criteria = 1),
    "non-missing character vector"
  )
  expect_error(
    geecriteria(fit_geewa_pois_exch, criteria = NA_character_),
    "non-missing character vector"
  )
})

test_that("restricting the criteria does not change the values reported", {
  full <- geecriteria(fit_geewa_pois_exch, digits = 15)
  for (criterion in geer_criteria_columns) {
    one <- geecriteria(
      fit_geewa_pois_exch,
      criteria = criterion,
      digits = 15
    )
    expect_named(one, c(criterion, "Parameters"))
    expect_equal(one[[criterion]], full[[criterion]], tolerance = 1e-10)
  }
})

test_that("unrequested criteria are not computed", {
  out <- compute_criteria(
    fit_geewa_pois_exch,
    cov_type = "robust",
    criteria = c("GESSC", "GPC")
  )
  expect_named(out, c(geer_criteria_columns, "Parameters"))
  expect_true(all(is.finite(c(out$GESSC, out$GPC))))
  skipped <- setdiff(geer_criteria_columns, c("GESSC", "GPC"))
  expect_true(all(is.na(unlist(out[skipped]))))
})

test_that("repeated models do not collide in the row names", {
  fit <- fit_geewa_pois_exch
  out <- geecriteria(fit, fit)
  expect_equal(nrow(out), 2L)
  expect_false(anyDuplicated(rownames(out)) > 0L)
  expect_equal(out[1L, "CIC"], out[2L, "CIC"])
})

test_that("GESSC and GPC are addressed by name, not position", {
  object <- fit_geewa_pois_exch
  stats <- get_working_covariance_criteria_cc(
    object$y,
    object$id,
    object$repeated,
    object$family$family,
    object$fitted.values,
    object$association_structure,
    object$alpha,
    object$phi,
    object$prior.weights
  )
  expect_named(stats, c("sc", "gp"))

  p <- length(coef(object))
  m <- get_association_count(object)
  out <- geecriteria(object, digits = 15)
  expect_equal(
    out$GESSC,
    stats$sc / (object$obs_no - p - m),
    tolerance = 1e-10
  )
  expect_equal(out$GPC, stats$gp, tolerance = 1e-10)
})

test_that("PT, WR and RMR match a direct generalized eigenvalue calculation", {
  object <- fit_geewa_pois_exch
  out <- geecriteria(object, digits = 15)

  omega_i <- compute_independence_info(object)
  covariance <- vcov(object, cov_type = "robust")
  lambda <- sort(Re(eigen(covariance %*% omega_i, only.values = TRUE)$values))
  expect_true(all(lambda > 0))
  ratios <- lambda / (1 + lambda)

  expect_equal(out$PT, sum(ratios), tolerance = 1e-8)
  expect_equal(out$WR, prod(ratios), tolerance = 1e-8)
  expect_equal(out$RMR, max(ratios), tolerance = 1e-8)

  reference <- covariance %*% solve(covariance + solve(omega_i))
  expect_equal(out$PT, sum(diag(reference)), tolerance = 1e-8)
  expect_equal(out$WR, det(reference), tolerance = 1e-8)
})

test_that("the eigenvalue criteria satisfy their structural bounds", {
  out <- geecriteria(fit_geewa_pois_exch, digits = 15)
  p <- length(coef(fit_geewa_pois_exch))

  expect_gt(out$RMR, 0)
  expect_lt(out$RMR, 1)
  expect_gt(out$WR, 0)
  expect_lt(out$WR, 1)
  expect_lte(out$RMR, out$PT)
  expect_lt(out$PT, p)
  expect_lte(out$WR, out$RMR)
})

test_that("the eigenvalue criteria respond to cov_type", {
  robust <- geecriteria(fit_geewa_pois_exch, cov_type = "robust", digits = 15)
  naive <- geecriteria(fit_geewa_pois_exch, cov_type = "naive", digits = 15)

  expect_false(isTRUE(all.equal(robust$PT, naive$PT)))
  expect_false(isTRUE(all.equal(robust$WR, naive$WR)))
  expect_false(isTRUE(all.equal(robust$RMR, naive$RMR)))
})

test_that("the eigenvalue criteria return NA for a degenerate reference", {
  expect_identical(
    compute_jang_criteria(
      beta_covariance = diag(2),
      independence_inverse = matrix(c(1, 2, 2, 1), 2L, 2L)
    ),
    list(PT = NA_real_, WR = NA_real_, RMR = NA_real_)
  )
})

test_that("GHYC and PAC are reported only for balanced designs", {
  balanced <- geecriteria(fit_binary_indep, digits = 15)
  expect_true(is.finite(balanced$GHYC))
  expect_true(is.finite(balanced$PAC))

  cluster_sizes <- as.integer(table(fit_binary_small$id))
  expect_gt(length(unique(cluster_sizes)), 1L)
  unbalanced <- geecriteria(fit_binary_small, digits = 15)
  expect_true(is.na(unbalanced$GHYC))
  expect_true(is.na(unbalanced$PAC))

  expect_true(is.finite(unbalanced$QIC))
  expect_true(is.finite(unbalanced$CIC))
  expect_true(is.finite(unbalanced$RJC))
})

test_that("working covariance helper supports correlation GEE", {
  object <- fit_geewa_pois_exch
  indices <- split(seq_along(object$id), object$id)[[1L]]
  expect_equal(
    get_working_covariance(object, indices),
    ref_fit_covariance(object, indices),
    tolerance = 1e-12
  )
})

test_that("odds-ratio pair indexing matches upper-triangular ordering", {
  expect_equal(
    vapply(
      list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4)),
      function(pair) {
        get_pair_index(
          pair[[1L]], pair[[2L]], 4L
        )
      },
      integer(1)
    ),
    seq_len(6L)
  )
})

test_that("the quasi-log-likelihood differs from its value at mu = y by minus half the unit deviance", {
  cases <- list(
    gaussian = list(family = gaussian(), y = 1.3, mu = 0.7),
    binomial = list(family = binomial(), y = 0.3, mu = 0.6),
    poisson = list(family = poisson(), y = 2.5, mu = 4),
    Gamma = list(family = Gamma(), y = 1.3, mu = 2.2),
    inverse.gaussian = list(family = inverse.gaussian(), y = 1.3, mu = 2.2)
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    q_mu <- compute_quasi_loglikelihood(
      y = case$y, mu = case$mu, weights = 1, family_name = name, phi = 1
    )
    q_y <- compute_quasi_loglikelihood(
      y = case$y, mu = case$y, weights = 1, family_name = name, phi = 1
    )
    deviance <- case$family$dev.resids(case$y, case$mu, 1)
    expect_equal(q_mu - q_y, -deviance / 2, tolerance = 1e-10, info = name)
    expect_lt(q_mu, q_y)
  }
})

test_that("geecriteria stops when models differ in family or response", {
  fit_pois <- fit_geewa_pois_exch
  fit_gauss <- geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = gaussian(link = "identity"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_error(
    geecriteria(fit_pois, fit_gauss),
    "models differ in the family (poisson, gaussian)",
    fixed = TRUE
  )
  fit_other_y <- fit_pois
  fit_other_y$y <- fit_other_y$y + 1
  expect_error(
    geecriteria(fit_pois, fit_other_y),
    "models differ in the response",
    fixed = TRUE
  )
  expect_no_warning(geecriteria(fit_pois, fit_pois))
})

test_that("odds-ratio working covariance matches high-precision values", {
  cases <- data.frame(
    p = c(1e-6, 1e-2, 1e-6, 0.3, 0.9),
    q = c(1e-6, 1e-6, 1e-2, 0.3, 0.8),
    psi = c(1 + 1e-6, 1.01, 2, 0.5, 0.01),
    reference = c(
      9.99997999916733467e-19,
      9.89900010298930725e-11,
      9.80196098035314672e-09,
      -2.84226894136091687e-02,
      -1.97156188675706827e-02
    )
  )
  for (i in seq_len(nrow(cases))) {
    expect_equal(
      compute_plackett_covariance(cases$p[i], cases$q[i], cases$psi[i]),
      cases$reference[i],
      tolerance = 1e-12
    )
  }
  expect_identical(compute_plackett_covariance(0.3, 0.4, 1), 0)
  expect_identical(compute_plackett_covariance(1e-9, 0.4, 5), 0)
})

test_that("odds-ratio working covariance in R matches the C++ engine", {
  fit <- fit_geewa_bin_exch
  clusters <- split(seq_along(fit$id), fit$id)
  sc <- 0
  log_det <- 0
  for (indices in clusters) {
    v_matrix <- compute_working_covariance_or(fit, indices)
    residual <- fit$y[indices] - fit$fitted.values[indices]
    sc <- sc + sum(residual * solve(v_matrix, residual))
    log_det <- log_det +
      as.numeric(determinant(v_matrix, logarithm = TRUE)$modulus)
  }
  expected <- get_working_covariance_criteria_or(
    fit$y,
    fit$id,
    fit$repeated,
    fit$fitted.values,
    get_alpha_or(fit),
    fit$prior.weights
  )
  expect_equal(sc, expected$sc, tolerance = 1e-10)
  expect_equal(-0.5 * (sc + log_det), expected$gp, tolerance = 1e-10)
})

test_that("geecriteria labels models supplied through do.call by position", {
  fit_ar1 <- update(fit_geewa_pois_exch, corstr = "ar1")
  out <- do.call(
    geecriteria,
    list(fit_geewa_pois_exch, fit_ar1, criteria = "QIC")
  )
  expect_identical(rownames(out), c("Model 1", "Model 2"))
  out_named <- geecriteria(fit_geewa_pois_exch, fit_ar1, criteria = "QIC")
  expect_identical(rownames(out_named), c("fit_geewa_pois_exch", "fit_ar1"))
})
