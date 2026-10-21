fit_agresti <- geewa(
  formula = abnormal ~ time + diagnosis + treatment,
  data = test_data$agresti,
  id = id,
  repeated = time,
  family = binomial(link = "logit"),
  corstr = "independence",
  phi_fixed = TRUE,
  phi_value = 1
)

fit_resp_or <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  data = test_data$respiratory2,
  id = id,
  repeated = visit,
  link = "logit",
  orstr = "exchangeable"
)

ref_working_covariance <- function(fit, i) {
  mu <- fit$fitted.values[i]
  if (!identical(fit$association_structure, "ar1")) {
    psi <- fit$alpha[1L]
    pairs <- outer(mu, mu, "+")
    s <- 1 + pairs * (psi - 1)
    joint <- (s - sqrt(s^2 - 4 * psi * (psi - 1) * outer(mu, mu))) /
      (2 * (psi - 1))
    out <- joint - outer(mu, mu)
    diag(out) <- mu * (1 - mu)
    return(out)
  }
  positions <- fit$repeated[i]
  fit$phi * outer(sqrt(mu * (1 - mu)), sqrt(mu * (1 - mu))) *
    fit$alpha^abs(outer(positions, positions, "-"))
}

ref_gof_residual <- function(fit, statistic, covariance) {
  mu <- fit$fitted.values
  r <- fit$y - mu
  v <- mu * (1 - mu)
  d_matrix <- fit$x * fit$family$mu.eta(fit$linear.predictors)
  clusters <- split(seq_along(fit$id), fit$id)
  n <- length(mu)
  positions_no <- max(fit$repeated)
  standardized <- matrix(NA_real_, length(clusters), positions_no)
  for (k in seq_along(clusters)) {
    i <- clusters[[k]]
    standardized[k, fit$repeated[i]] <- r[i] / sqrt(v[i])
  }
  pooled <- matrix(0, positions_no, positions_no)
  for (j in seq_len(positions_no)) {
    for (l in seq_len(positions_no)) {
      pooled[j, l] <- mean(standardized[, j] * standardized[, l], na.rm = TRUE)
    }
  }
  v_inverse <- matrix(0, n, n)
  cov_y <- matrix(0, n, n)
  for (i in clusters) {
    v_inverse[i, i] <- solve(ref_working_covariance(fit, i))
    positions <- fit$repeated[i]
    cov_y[i, i] <- if (identical(covariance, "empirical")) {
      tcrossprod(r[i])
    } else {
      diag(sqrt(v[i])) %*% pooled[positions, positions] %*% diag(sqrt(v[i]))
    }
  }
  sensitivity <- crossprod(d_matrix, v_inverse %*% d_matrix)
  h_matrix <- d_matrix %*% solve(sensitivity, t(d_matrix) %*% v_inverse)
  if (identical(statistic, "pearson")) {
    value <- sum(r^2 / v)
    mean_value <- n
    weights <- (1 - 2 * mu) / v
  } else {
    value <- sum(r^2)
    mean_value <- sum(v)
    weights <- 1 - 2 * mu
  }
  adjusted <- crossprod(diag(n) - h_matrix, weights)
  variance <- drop(crossprod(adjusted, cov_y %*% adjusted))
  z <- (value - mean_value) / sqrt(variance)
  c(z = z, p = 2 * pnorm(-abs(z)))
}

test_that("the covariate-region test reproduces Barnhart and Williamson", {
  main <- gof_score_test(
    fit_agresti,
    partition = "covariates",
    regions = ~ diagnosis + treatment,
    cov_type = "naive"
  )
  expect_lt(abs(unname(main$statistic) - 33.83), 0.005)
  expect_identical(unname(main$parameter), 8L)
  interaction <- gof_score_test(
    update(fit_agresti, . ~ . + treatment:time),
    partition = "covariates",
    regions = ~ diagnosis + treatment,
    cov_type = "naive"
  )
  expect_lt(abs(unname(interaction$statistic) - 4.45), 0.005)
  expect_identical(unname(interaction$parameter), 7L)
  expect_equal(sum(main$table$observations), fit_agresti$obs_no)
  expect_equal(sum(main$table$observed), sum(fit_agresti$y))
  expect_equal(sum(main$table$expected), sum(fit_agresti$fitted.values))
})

test_that("the partition tests equal score tests of the alternative model", {
  risk <- gof_score_test(fit_agresti, n_groups = 6)
  indicators <- model.matrix(~ 0 + risk$partition)[, -1L]
  colnames(indicators) <- paste0("g", seq_len(ncol(indicators)))
  augmented_data <- cbind(test_data$agresti, indicators)
  augmented <- update(
    fit_agresti,
    formula = as.formula(paste(
      "abnormal ~ time + diagnosis + treatment +",
      paste(colnames(indicators), collapse = " + ")
    )),
    data = augmented_data
  )
  restricted <- update(fit_agresti, data = augmented_data)
  expect_equal(
    unname(gof_score_test(fit_agresti, n_groups = 6,
                          cov_type = "robust")$statistic),
    anova(restricted, augmented, test = "score", cov_type = "robust")$Chi[2L],
    tolerance = 1e-8
  )
  mu <- fit_agresti$fitted.values
  expect_identical(
    levels(risk$partition),
    levels(cut(mu, unique(quantile(mu, 0:6 / 6)), include.lowest = TRUE))
  )
  by_formula <- gof_score_test(fit_agresti, partition = "covariates",
                               regions = ~ diagnosis + treatment)
  regions <- with(
    test_data$agresti,
    interaction(diagnosis, treatment, lex.order = TRUE)
  )
  expect_equal(
    gof_score_test(fit_agresti, partition = "covariates", regions = regions),
    by_formula,
    ignore_attr = TRUE
  )
  expect_equal(
    gof_score_test(fit_agresti, partition = "covariates",
                   regions = ~ diagnosis:treatment),
    by_formula
  )
  local_regions <- function(fit) {
    cut_points <- c(-0.5, 0.5, 1.5)
    gof_score_test(fit, partition = "covariates",
                   regions = ~ diagnosis + cut(treatment, cut_points))
  }
  expect_equal(local_regions(fit_agresti)$statistic, by_formula$statistic)
})

test_that("the residual tests agree with a direct implementation", {
  probit <- geewa(
    formula = status ~ baseline + treatment + age,
    data = test_data$respiratory2,
    id = id,
    repeated = visit,
    family = binomial(link = "probit"),
    corstr = "ar1"
  )
  unbalanced <- update(fit_resp_or, data = test_data$respiratory2[-c(1, 6), ])
  cases <- list(
    list(fit_resp_or, "pearson", "unstructured"),
    list(fit_resp_or, "pearson", "empirical"),
    list(fit_resp_or, "uss", "unstructured"),
    list(fit_resp_or, "uss", "empirical"),
    list(probit, "pearson", "unstructured"),
    list(unbalanced, "pearson", "unstructured"),
    list(unbalanced, "pearson", "empirical")
  )
  for (case in cases) {
    out <- gof_residual_test(case[[1L]], statistic = case[[2L]],
                             covariance = case[[3L]])
    expected <- ref_gof_residual(case[[1L]], case[[2L]], case[[3L]])
    expect_equal(unname(out$statistic), unname(expected["z"]),
                 tolerance = 1e-8)
  }
  expect_equal(out$p.value, unname(expected["p"]), tolerance = 1e-8)
})

test_that("gof_tests collects the individual tests", {
  out <- gof_tests(fit_resp_or, regions = ~treatment)
  expect_s3_class(out, "geer_anova")
  expect_named(out, c("Statistic", "Df", "Mean", "Variance", "Pr(Test)"))
  expect_identical(nrow(out), 5L)
  expect_equal(
    out[1L, "Statistic"],
    unname(gof_score_test(fit_resp_or)$statistic)
  )
  expect_equal(
    out[2L, "Pr(Test)"],
    gof_score_test(fit_resp_or, partition = "covariates",
                   regions = ~treatment, cov_type = "naive")$p.value
  )
  pearson <- gof_residual_test(fit_resp_or)
  expect_equal(out[4L, c("Statistic", "Mean", "Variance")],
               data.frame(Statistic = unname(pearson$estimate),
                          Mean = pearson$mean,
                          Variance = pearson$variance),
               ignore_attr = TRUE)
  expect_identical(nrow(gof_tests(fit_resp_or, cov_type = "naive")), 3L)
  expect_identical(
    attr(out, "heading")[1L],
    "Goodness-of-Fit Tests for Binary Responses\n"
  )
})

test_that("goodness-of-fit tests check the fit and their arguments", {
  fit_brgee <- update(fit_resp_or, method = "brgee-robust")
  expect_warning(
    out_brgee <- gof_tests(fit_brgee),
    "were developed for method = \"gee\"; for method = \"brgee-robust\""
  )
  expect_true(all(is.finite(out_brgee[["Pr(Test)"]])))
  expect_error(
    gof_score_test(fit_brgee, regions = ~treatment),
    "'regions' and 'by_repeated' are used only with partition = \"covariates\"",
    fixed = TRUE
  )
  expect_warning(
    gof_score_test(fit_geewa_bin_exch),
    "only 4 groups of fitted probabilities could be formed",
    fixed = TRUE
  )
  components <- compute_gof_residual_parts(fit_resp_or, "empirical")
  components$covariances <- lapply(components$covariances, `*`, 0)
  expect_error(
    compute_gof_residual_test(components, "pearson", "gof_residual_test"),
    "the approximate variance of the statistic is not positive",
    fixed = TRUE
  )
  sparse <- list(
    id = rep(1:60, each = 2),
    repeated = rep(1:2, 60),
    y = rep(c(1, 0), 60)
  )
  expect_warning(
    check_gof_regions(sparse, factor(rep(1:2, each = 60)), TRUE),
    "some region by position cells have no events or no non-events",
    fixed = TRUE
  )
  expect_silent(
    check_gof_regions(sparse, factor(rep(1:2, each = 60)), FALSE)
  )
  expect_error(
    gof_residual_test(fit_geewa_pois_exch),
    "require a binary (0/1) response with unit prior weights",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, cov_type = "jackknife"),
    "cov_type = \"jackknife\" is not available",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates"),
    "'regions' must be supplied",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates", n_groups = 5,
                   regions = ~treatment),
    "'n_groups' is used only with partition = \"risk\"",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates", regions = 1:3),
    "'regions' must be a one-sided formula or a vector",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates",
                   regions = ~treatment, by_repeated = FALSE),
    "the partition adds no column that is linearly independent",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates", regions = ~1),
    "'regions' must contain at least one variable",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates",
                   regions = ~unknown_variable),
    "'regions' cannot be evaluated in the data used to fit the model",
    fixed = TRUE
  )
  incomplete <- test_data$respiratory2
  incomplete$region <- ifelse(incomplete$id == incomplete$id[1L], NA, 1)
  fit_incomplete <- update(fit_resp_or, data = incomplete)
  expect_error(
    gof_score_test(fit_incomplete, partition = "covariates",
                   regions = ~region),
    "the variables in 'regions' have missing values for observations used",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(fit_resp_or, partition = "covariates",
                   regions = rep(1, fit_resp_or$obs_no)),
    "'regions' must define at least two regions",
    fixed = TRUE
  )
  expect_error(
    gof_tests(fit_resp_or, by_repeated = FALSE),
    "'by_repeated' is used only when 'regions' is supplied",
    fixed = TRUE
  )
  expect_error(
    gof_tests(fit_resp_or, by_repeated = NA),
    "'by_repeated' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    gof_score_test(update(fit_resp_or, . ~ 1)),
    "the fitted probabilities take too few distinct values",
    fixed = TRUE
  )
  expect_warning(
    gof_score_test(fit_resp_or, partition = "covariates",
                   regions = ~ treatment + gender),
    "some regions have fewer than 10 clusters"
  )
})
