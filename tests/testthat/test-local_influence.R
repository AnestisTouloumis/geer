test_that("the observed derivative of the estimating equations is exact", {
  for (fit in list(fit_gamma_weighted, fit_geewa_bin_exch)) {
    quantities <- compute_influence_quantities(fit)
    jacobian <- ref_jac(
      function(beta) ref_estimating_equations(fit, beta),
      coef(fit)
    )
    expect_equal(quantities$sensitivity, -jacobian, tolerance = 1e-6)
  }
})

test_that("local influence agrees with a direct implementation", {
  for (fit in list(fit_gamma_weighted, fit_geewa_bin_exch)) {
    reference <- ref_influence(fit)
    sensitivity <- compute_influence_quantities(fit)$sensitivity
    inverse <- solve(sensitivity)
    deltas <- list(
      "cluster-weights" = reference$cluster_delta,
      "observation-weights" = reference$observation_delta,
      response = reference$response_delta
    )
    for (perturbation in names(deltas)) {
      delta <- deltas[[perturbation]]
      expect_equal(
        unname(local_influence(fit, perturbation = perturbation)),
        unname(2 * abs(rowSums((delta %*% inverse) * delta))),
        tolerance = 1e-8
      )
    }
    delta <- deltas[["cluster-weights"]]
    curvature <- delta %*% inverse %*% t(delta)
    curvature <- 0.5 * (curvature + t(curvature))
    decomposition <- eigen(curvature, symmetric = TRUE)
    leading <- decomposition$vectors[, which.max(abs(decomposition$values))]
    expect_equal(
      unname(local_influence(fit, type = "local")),
      abs(leading),
      tolerance = 1e-8
    )
  }
})

test_that("local influence for a subset of coefficients removes the others", {
  reference <- ref_influence(fit_geewa_pois_exch)
  sensitivity <- compute_influence_quantities(
    fit_geewa_pois_exch
  )$sensitivity
  nuisance <- c(1L, 3L, 4L)
  inverse <- solve(sensitivity)
  inverse[nuisance, nuisance] <- inverse[nuisance, nuisance] -
    solve(sensitivity[nuisance, nuisance])
  delta <- reference$cluster_delta
  out <- local_influence(fit_geewa_pois_exch, parm = "treatmentprogabide")
  expect_equal(
    unname(out),
    unname(2 * abs(rowSums((delta %*% inverse) * delta))),
    tolerance = 1e-8
  )
})

test_that("local influence validates its arguments", {
  expect_error(local_influence(fit_geewa_pois_exch, perturbation = "cw"),
               "'perturbation' must be one of")
  expect_error(local_influence(fit_geewa_pois_exch, type = "max"),
               "'type' must be one of: total, local")
})

test_that("rank-deficient perturbations retain all coefficient columns", {
  fit <- geewa(formula = y ~ x + z, data = test_data$influence_deficient,
               id = id, phi_fixed = TRUE)
  design <- model.matrix(fit)
  delta <- rowsum(design * (fit$y - fitted(fit)), group = fit$id)
  expect_lt(qr(delta)$rank, ncol(delta))
  curvature <- delta %*% solve(crossprod(design)) %*% t(delta)
  decomposition <- eigen(curvature, symmetric = TRUE)
  expected <- abs(decomposition$vectors[
    , which.max(abs(decomposition$values))
  ])
  expect_equal(unname(local_influence(fit, type = "local")), expected,
               tolerance = 1e-6)
})

test_that("a zero perturbation matrix has no local direction", {
  expect_error(
    compute_local_direction(matrix(0, 3L, 2L), diag(2L)),
    "local_influence failed: the perturbation matrix has zero rank",
    fixed = TRUE
  )
})
