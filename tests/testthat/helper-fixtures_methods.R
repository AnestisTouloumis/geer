test_methods_fixtures <- list(geer_methods_add1_drop1_anova = function() {
  fit_null <- geewa(
    formula = seizures ~ 1, data = test_data$epilepsy,
    id = id, family = poisson(link = "log"), corstr = "independence",
    method = "gee"
  )
  fit_trt <- geewa(
    formula = seizures ~ treatment,
    data = test_data$epilepsy,
    id = id, family = poisson(link = "log"), corstr = "independence",
    method = "gee"
  )
  as.list(environment(), all.names = TRUE)
}, geer_methods_emmeans = function() {
  skip_if_no_emmeans <- function() {
    skip_if_not_installed("emmeans", minimum_version = "2.0.0")
  }
  as.list(environment(), all.names = TRUE)
}, geer_methods_influence = function() {
  gaussian_sim <- ref_sim("cc", rep(list(1:4), 20),
    family = gaussian(),
    link = "identity"
  )
  fit_gaussian_ar1 <- geewa(
    formula = y ~ x1 + x2, data = gaussian_sim$df,
    id = id, repeated = visit, family = gaussian(link = "identity"),
    corstr = "ar1", control = geer_control(tolerance = 1e-10)
  )
  expect_influence_reference <- function(fit) {
    reference <- ref_influence(fit)
    cluster_hat <- vapply(split(reference$hat, fit$id), sum, numeric(1))
    expect_equal(unname(hatvalues(fit)), unname(cluster_hat),
      tolerance = 1e-10
    )
    expect_equal(unname(hatvalues(fit, unit = "observations")),
      reference$hat,
      tolerance = 1e-10
    )
    expect_equal(unname(dfbeta(fit)), unname(reference$cluster_dfbeta),
      tolerance = 1e-10
    )
    expect_equal(unname(dfbeta(fit, unit = "observations")),
      unname(reference$observation_dfbeta),
      tolerance = 1e-10
    )
    bread_inverse <- solve(reference$bread)
    expect_equal(unname(cooks.distance(fit)),
      unname(rowSums(
        (reference$cluster_dfbeta %*% bread_inverse) *
          reference$cluster_dfbeta
      )) / ncol(bread_inverse),
      tolerance = 1e-10
    )
  }
  as.list(environment(), all.names = TRUE)
}, geer_methods_marginaleffects = function() {
  skip_if_no_marginaleffects <- function() {
    skip_if_not_installed("marginaleffects", minimum_version = "0.32.0")
  }
  as.list(environment(), all.names = TRUE)
}, geer_methods_plot_runs = function() {
  plot_runs_method <- getS3method("plot", "geer_runs_test")
  as.list(environment(), all.names = TRUE)
}, geer_output_snapshots = function() {
  fit_full <- geewa(seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id, repeated = visit
  )
  fit_reduced <- geewa(seizures ~ treatment,
    data = test_data$epilepsy,
    id = id, repeated = visit
  )
  as.list(environment(), all.names = TRUE)
})
