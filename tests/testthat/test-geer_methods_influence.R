gaussian_sim <- ref_sim(
  "cc",
  rep(list(1:4), 20),
  family = gaussian(),
  link = "identity"
)

fit_gaussian_ar1 <- geewa(
  formula = y ~ x1 + x2,
  data = gaussian_sim$df,
  id = id,
  repeated = visit,
  family = gaussian(link = "identity"),
  corstr = "ar1",
  control = geer_control(tolerance = 1e-10)
)

expect_influence_reference <- function(fit) {
  reference <- ref_influence(fit)
  cluster_hat <- vapply(split(reference$hat, fit$id), sum, numeric(1))
  expect_equal(unname(hatvalues(fit)), unname(cluster_hat), tolerance = 1e-10)
  expect_equal(
    unname(hatvalues(fit, unit = "observations")),
    reference$hat,
    tolerance = 1e-10
  )
  expect_equal(
    unname(dfbeta(fit)),
    unname(reference$cluster_dfbeta),
    tolerance = 1e-10
  )
  expect_equal(
    unname(dfbeta(fit, unit = "observations")),
    unname(reference$observation_dfbeta),
    tolerance = 1e-10
  )
  bread_inverse <- solve(reference$bread)
  expect_equal(
    unname(cooks.distance(fit)),
    unname(rowSums((reference$cluster_dfbeta %*% bread_inverse) *
                     reference$cluster_dfbeta)) / ncol(bread_inverse),
    tolerance = 1e-10
  )
}

test_that("leverages and one-step deletions match a direct implementation", {
  expect_influence_reference(fit_geewa_pois_exch)
  expect_influence_reference(fit_gamma_weighted)
  expect_influence_reference(fit_geewa_bin_exch)
})

test_that("one-step deletion estimates are exact for a gaussian identity fit", {
  expect_equal(
    dfbeta(fit_gaussian_ar1),
    dfbeta(fit_gaussian_ar1, method = "exact"),
    tolerance = 1e-6
  )
  expect_true(exists("jackknife_estimates", envir = fit_gaussian_ar1$cache))
  observation_dfbeta <- dfbeta(fit_gaussian_ar1, unit = "observations")
  for (j in c(1L, 6L, 80L)) {
    keep <- rep(TRUE, fit_gaussian_ar1$obs_no)
    keep[j] <- FALSE
    refit <- refit_jackknife_cc(fit_gaussian_ar1, keep, "observation")
    expect_equal(
      unname(observation_dfbeta[j, ]),
      unname(coef(fit_gaussian_ar1) - refit),
      tolerance = 1e-6
    )
  }
})

test_that("dfbetas and cooks.distance use the requested covariance matrix", {
  dfbeta_matrix <- dfbeta(fit_geewa_pois_exch, unit = "observations")
  robust <- vcov(fit_geewa_pois_exch, cov_type = "robust")
  expect_equal(
    dfbetas(fit_geewa_pois_exch, unit = "observations", cov_type = "robust"),
    sweep(dfbeta_matrix, 2L, sqrt(diag(robust)), `/`),
    tolerance = 1e-10
  )
  parm <- c("treatmentprogabide", "lnage")
  selected <- dfbeta_matrix[, parm]
  expect_equal(
    cooks.distance(
      fit_geewa_pois_exch,
      unit = "observations",
      cov_type = "robust",
      parm = parm
    ),
    rowSums((selected %*% solve(robust[parm, parm])) * selected) / 2,
    tolerance = 1e-10
  )
  expect_equal(
    cooks.distance(fit_geewa_pois_exch, parm = c(2L, 4L)),
    cooks.distance(fit_geewa_pois_exch, parm = parm)
  )
})

test_that("influence measures are named by the original identifiers", {
  data <- test_data$epilepsy
  data$patient <- paste0("p", data$id)
  fit <- geewa(
    formula = seizures ~ treatment + lnbaseline,
    data = data,
    id = patient,
    repeated = visit,
    family = poisson(link = "log"),
    corstr = "exchangeable"
  )
  expect_identical(fit$cluster_labels, sort(unique(data$patient)))
  expect_identical(fit$repeated_labels, as.character(sort(unique(data$visit))))
  codes <- as.numeric(names(split(fit$id, fit$id)))
  cluster_labels <- fit$cluster_labels[codes]
  expect_identical(names(hatvalues(fit)), cluster_labels)
  expect_identical(
    dimnames(dfbeta(fit)),
    list(cluster_labels, names(coef(fit)))
  )
  expect_identical(names(residuals(fit, type = "mahalanobis")),
                   cluster_labels)
  observation_labels <- paste0(
    fit$cluster_labels[fit$id], "(", fit$repeated_labels[fit$repeated], ")"
  )
  expect_identical(names(hatvalues(fit, unit = "observations")),
                   observation_labels)
  expect_identical(
    names(local_influence(fit, perturbation = "response")),
    observation_labels
  )
  legacy <- fit
  legacy$cluster_labels <- NULL
  legacy$repeated_labels <- NULL
  expect_identical(names(hatvalues(legacy)), as.character(codes))
  expect_identical(
    names(hatvalues(legacy, unit = "observations")),
    paste0(fit$id, "(", fit$repeated, ")")
  )
  no_repeated <- update(fit, repeated = NULL)
  expect_null(no_repeated$repeated_labels)
  expect_identical(
    names(cooks.distance(no_repeated, unit = "observations"))[1:2],
    c("p1(1)", "p1(2)")
  )
})

test_that("deletions that leave coefficients inestimable give NA", {
  data <- test_data$epilepsy
  data$rare <- factor(ifelse(data$id == 1, "a", "b"))
  data$spike <- as.numeric(seq_len(nrow(data)) == 2L)
  fit <- geewa(
    formula = seizures ~ rare + spike + lnbaseline,
    data = data,
    id = id,
    repeated = visit,
    family = poisson(link = "log"),
    corstr = "exchangeable"
  )
  expect_warning(
    out <- dfbeta(fit),
    paste0(
      "deleting cluster 1 leaves coefficients that cannot be estimated, so ",
      "the corresponding values are NA"
    ),
    fixed = TRUE
  )
  expect_true(all(is.na(out["1", ])))
  expect_false(anyNA(out[-1L, ]))
  expect_warning(
    out <- cooks.distance(fit, unit = "observations"),
    "deleting observation 1(2) leaves coefficients that cannot be estimated",
    fixed = TRUE
  )
  expect_identical(unname(which(is.na(out))), 2L)
  brgee <- update(fit, method = "brgee-robust")
  expect_warning(
    exact <- dfbeta(brgee, method = "exact"),
    "deleting cluster 1 leaves coefficients that cannot be estimated",
    fixed = TRUE
  )
  expect_true(all(is.na(exact["1", ])))
  expect_false(anyNA(exact[-1L, ]))
})

test_that("exact deletion errors name the calling function", {
  fit <- suppressWarnings(update(
    fit_geewa_pois_exch,
    formula = seizures ~ 1,
    data = test_data$epilepsy[test_data$epilepsy$id == 1, ]
  ))
  expect_error(
    dfbeta(fit, method = "exact"),
    "dfbeta.geer requires at least two clusters",
    fixed = TRUE
  )
})

test_that("influence measures stop for estimators without a one-step form", {
  fit <- update(fit_geewa_pois_exch, method = "pgee-jeffreys")
  expect_error(
    dfbeta(fit),
    paste0(
      "dfbeta.geer failed: method = \"one-step\" requires a fit with ",
      "method = \"gee\"; use method = \"exact\" with unit = \"clusters\""
    ),
    fixed = TRUE
  )
  expect_error(
    cooks.distance(fit, unit = "observations", method = "exact"),
    paste0(
      "cooks.distance.geer failed: method = \"exact\" is available only for ",
      "unit = \"clusters\""
    ),
    fixed = TRUE
  )
  expect_error(
    local_influence(fit),
    paste0(
      "local_influence failed: local influence requires a fit with ",
      "method = \"gee\""
    ),
    fixed = TRUE
  )
  expect_length(hatvalues(fit), fit$clusters_no)
})

test_that("influence measures validate their arguments", {
  expect_error(hatvalues(fit_geewa_pois_exch, unit = "cluster"),
               "'unit' must be one of: clusters, observations")
  expect_error(dfbeta(fit_geewa_pois_exch, method = "full"),
               "'method' must be one of: one-step, exact")
  expect_error(cooks.distance(fit_geewa_pois_exch, coefs = "age"),
               "does not use the argument")
})
