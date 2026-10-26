test_estimation_fixtures <- list(estimation_bias_reduction = function() {
  adjusted_residual <- function(d, method, ...) {
    fit <- ref_fit(d, method, ...)
    b <- unname(coef(fit))
    lam <- ref_lambda(d, method, b)
    list(converged = fit$converged, res = ref_U(d, b) + lam, lam = lam)
  }
  expect_adjusted_equation <- function(d, method, label, ...) {
    r <- adjusted_residual(d, method, ...)
    expect_true(r$converged, label = paste(label, "converged"))
    expect_lt(max(abs(r$res)), 1e-05 * max(1, max(abs(r$lam))))
  }
  cc_methods <- c(
    "brgee-robust", "brgee-naive", "brgee-empirical",
    "pgee-jeffreys"
  )
  joint_binary <- function(mu, lam) {
    m <- length(mu)
    Y <- as.matrix(expand.grid(rep(list(0:1), m)))
    pr <- combn(m, 2)
    pair <- rowSums(Y[, pr[1, ], drop = FALSE] * Y[, pr[2, ], drop = FALSE])
    th <- qlogis(mu)
    for (it in 1:500) {
      lp <- drop(Y %*% th) + lam * pair
      P <- exp(lp - max(lp))
      P <- P / sum(P)
      marg <- colSums(Y * P)
      th <- th + 0.8 * (qlogis(mu) - qlogis(pmin(
        pmax(marg, 1e-12),
        1 - 1e-12
      )))
    }
    list(Y = Y, P = drop(P))
  }
  as.list(environment(), all.names = TRUE)
}, estimation_convergence = function() {
  conv_df <- withr::with_seed(11, local({
    n <- 20
    df <- data.frame(
      id = rep(1:n, each = 4), visit = rep(1:4, n),
      x = rnorm(4 * n)
    )
    u <- rep(rnorm(n, sd = 0.7), each = 4)
    df$y <- rbinom(nrow(df), 1, plogis(-0.3 + 0.8 * df$x + u))
    df
  }))
  as.list(environment(), all.names = TRUE)
}, estimation_invariances = function() {
  all_methods <- c(
    "gee", "brgee-robust", "brgee-naive", "brgee-empirical",
    "bcgee-robust", "bcgee-naive", "bcgee-empirical", "pgee-jeffreys"
  )
  times6 <- list(1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3)
  as.list(environment(), all.names = TRUE)
}, geewa_binary_theory = function() {
  independence_formula <- status ~ baseline + treatment + gender + visit +
    age + center
  fit_independence_pair <- function(
    method, link = "probit",
    control = ctl_tight
  ) {
    list(
      or = geewa_binary(
        formula = independence_formula, id = id,
        repeated = visit,
        data = test_data$respiratory, link = link, orstr = "independence",
        method = method, control = control
      ),
      cor = geewa(
        formula = independence_formula,
        id = id, repeated = visit, family = binomial(link = link),
        data = test_data$respiratory,
        phi_fixed = TRUE, phi_value = 1, corstr = "independence",
        method = method,
        control = control
      )
    )
  }
  as.list(environment(), all.names = TRUE)
}, geewa_fit = function() {
  rinse <- test_data$rinse
  fit_geewa_gauss_indep <- geewa(
    formula = score ~ treatment + baseline +
      time, data = rinse, id = id, family = gaussian(link = "identity"),
    corstr = "independence", method = "gee"
  )
  as.list(environment(), all.names = TRUE)
}, geewa_theory = function() {
  fit_gee_probit <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "probit"), phi_fixed = TRUE, phi_value = 1,
    data = test_data$cerebrovascular, corstr = "independence",
    method = "gee",
    control = ctl_tight
  )
  fit_glm_probit <- glm(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "probit"),
    data = test_data$cerebrovascular, control = list(
      epsilon = 1e-14,
      maxit = 100
    )
  )
  fit_gee_logit <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "logit"), phi_fixed = TRUE, phi_value = 1,
    data = test_data$cerebrovascular, corstr = "independence",
    method = "gee",
    control = ctl_tight
  )
  fit_glm_logit <- glm(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"),
    data = test_data$cerebrovascular, control = list(
      epsilon = 1e-14,
      maxit = 100
    )
  )
  fit_brglm <- function(object, type) {
    update(object,
      method = brglm2::brglmFit,
      control = list(
        epsilon = 1e-14,
        maxit = 100, type = type
      )
    )
  }
  fit_brgee_naive_logit <- update(fit_gee_logit, method = "brgee-naive")
  fit_pgee_jeffreys_logit <- update(fit_brgee_naive_logit,
    method = "pgee-jeffreys"
  )
  fit_brgee_naive_probit <- update(fit_gee_logit,
    family = binomial(link = "probit"),
    method = "brgee-naive"
  )
  fit_pgee_jeffreys_probit <- update(fit_brgee_naive_probit,
    method = "pgee-jeffreys"
  )
  as.list(environment(), all.names = TRUE)
})
