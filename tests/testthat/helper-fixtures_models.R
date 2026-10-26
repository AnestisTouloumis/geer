test_models_fixtures <- list(geecriteria = function() {
  fit_binary_indep <- geewa_binary(
    formula = ecg ~ period * treatment,
    id = id, data = test_data$cerebrovascular, link = "logit",
    orstr = "independence",
    method = "gee"
  )
  fit_binary_exch <- update(fit_binary_indep, orstr = "exchangeable")
  cerebrovascular_small <- test_data$cerebrovascular[
    seq_len(nrow(test_data$cerebrovascular) - 1L), ,
    drop = FALSE
  ]
  fit_binary_small <- geewa_binary(
    formula = ecg ~ period * treatment,
    id = id, data = cerebrovascular_small, link = "logit",
    orstr = "independence",
    method = "gee"
  )
  as.list(environment(), all.names = TRUE)
}, geer_constructor = function() {
  make_valid_geer_input <- function() {
    dat <- data.frame(y = c(1, 2, 3), x1 = c(0, 1, 0), x2 = c(
      1, 0.5,
      -0.5
    ), id = c(1, 1, 2), repeated = c(1, 2, 1))
    form <- y ~ x1 + x2
    xmat <- model.matrix(form, dat)
    beta <- c(`(Intercept)` = 0.2, x1 = 0.5, x2 = -0.1)
    fitted <- as.numeric(xmat %*% beta)
    resid <- dat$y - fitted
    V <- diag(c(0.1, 0.2, 0.3))
    dimnames(V) <- list(names(beta), names(beta))
    list(
      coefficients = beta, residuals = resid, fitted.values = fitted,
      qr = qr(xmat), rank = ncol(xmat), family = gaussian(),
      linear.predictors = fitted,
      iter = 3L, prior.weights = rep(1, nrow(dat)),
      df.residual = nrow(dat) -
        ncol(xmat), y = dat$y, x = xmat, id = dat$id,
      repeated = dat$repeated,
      converged = TRUE, call = quote(geewa(formula = y ~ x1 + x2)),
      formula = form, terms = terms(form), data = dat, offset = NULL,
      control = list(), method = "gee", naive_covariance = V,
      robust_covariance = V,
      bias_corrected_covariance = V, association_structure = "exchangeable",
      alpha = 0.2, phi = 1, obs_no = nrow(dat),
      clusters_no = length(unique(dat$id)),
      min_cluster_size = 1L, max_cluster_size = 2L
    )
  }
  make_valid_geer_object <- function() {
    new_geer(make_valid_geer_input())
  }
  as.list(environment(), all.names = TRUE)
}, geer_fit_inputs = function() {
  fit_inputs_data <- data.frame(
    id = rep(1:3, each = 2), t = rep(1:2,
      times = 3
    ),
    y = c(1.2, 0.4, 2.2, 1.9, 0.7, 3.1), x = c(
      0.1, 0.5, 0.2, 0.9, 0.4,
      0.8
    )
  )
  fit_or_offset_formula <- geewa_binary(
    formula = ecg ~ treatment + offset(period),
    link = "logit", id = id, data = test_data$cerebrovascular,
    orstr = "independence",
    method = "gee"
  )
  fit_or_offset_argument <- geewa_binary(
    formula = ecg ~ treatment,
    link = "logit",
    id = id, offset = period, data = test_data$cerebrovascular,
    orstr = "independence",
    method = "gee"
  )
  fit_cc_offset_formula <- geewa(
    formula = ecg ~ treatment + offset(period),
    family = binomial(link = "logit"), id = id,
    data = test_data$cerebrovascular,
    corstr = "independence", method = "gee", phi_fixed = TRUE,
    phi_value = 1
  )
  as.list(environment(), all.names = TRUE)
}, geer_fit_output = function() {
  finalize_fit <- function(
    criterion = 1e-08, fitted = c(0.2, 0.5, 0.8, 0.3, 0.6, 0.7),
    method = "gee", family = binomial(), structure = "exchangeable",
    repeated = c(1, 2, 3, 1, 2, 3), alpha = 0.3, fit_function = "geewa",
    tolerance = 1e-06, iter = 1L, failure = NULL
  ) {
    build_geer_final_fit(
      fit = list(iter = iter, fitted.values = fitted),
      geesolver_fit = list(
        criterion = criterion, alpha = alpha,
        failure = failure
      ),
      tolerance = tolerance, method = method, family = family,
      association_structure = structure,
      repeated = repeated, fit_function = fit_function
    )
  }
  as.list(environment(), all.names = TRUE)
}, geer_hypothesis_tests = function() {
  fit_bin_trt <- geewa(
    formula = ecg ~ treatment, id = id,
    family = binomial(link = "logit"),
    phi_fixed = TRUE, phi_value = 1, data = test_data$cerebrovascular,
    corstr = "independence", method = "gee"
  )
  fit_bin_period <- geewa(
    formula = ecg ~ factor(period), id = id,
    family = binomial(link = "logit"),
    phi_fixed = TRUE, phi_value = 1, data = test_data$cerebrovascular,
    corstr = "independence", method = "gee"
  )
  fit_bin_full <- geewa(
    formula = ecg ~ treatment + factor(period), id = id,
    family = binomial(link = "logit"), phi_fixed = TRUE, phi_value = 1,
    data = test_data$cerebrovascular, corstr = "independence",
    method = "gee"
  )
  fit_bin_full_exch <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id, family = binomial(link = "logit"), phi_fixed = TRUE,
    phi_value = 1,
    data = test_data$cerebrovascular, corstr = "exchangeable",
    method = "gee"
  )
  fixed_null_fit <- function(d, formula) {
    geewa(formula,
      data = d$df, id = id, repeated = visit,
      family = d$family,
      corstr = "fixed", alpha_vector = ref_corr(
        d$struct, d$alpha,
        d$tmax
      )[lower.tri(diag(d$tmax))], phi_fixed = TRUE,
      phi_value = d$phi,
      weights = w, control = ctl_tight
    )
  }
  as.list(environment(), all.names = TRUE)
}, geer_jackknife = function() {
  jackknife_gaussian_data <- data.frame(
    id = rep(seq_len(8), each = 3),
    time = rep(seq_len(3), times = 8), group = rep(rep(c(0, 1), each = 4),
      each = 3
    ), y = c(
      1.2, 1.7, 2.4, 2, 2.2, 2.9, 1.5, 2.1, 2.5, 2.3,
      2.7, 3.2, 2.1, 2.8, 3.4, 2.7, 3, 3.8, 2.4, 2.6, 3.5, 3, 3.5,
      3.9
    )
  )
  jackknife_binary_data <- data.frame(
    id = rep(seq_len(10), each = 3),
    time = rep(seq_len(3), times = 10), group = rep(rep(c(0, 1), each = 5),
      each = 3
    ), y = c(
      0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0,
      0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0
    )
  )
  as.list(environment(), all.names = TRUE)
}, geer_small_sample_covariances = function() {
  expect_small_sample_reference <- function(fit) {
    reference <- ref_small_sample_covariances(fit)
    for (cov_type in geer_ss_cov_types) {
      expect_equal(unname(vcov(fit, cov_type = cov_type)),
        unname(reference[[cov_type]]),
        tolerance = 1e-10, label = cov_type
      )
    }
  }
  as.list(environment(), all.names = TRUE)
}, geer_start_values = function() {
  start_model_matrix <- function() {
    model.matrix(seizures ~ treatment + lnbaseline + lnage,
      data = test_data$epilepsy
    )
  }
  start_values <- function(
    beta_start = NULL, family = poisson(),
    link = "log",
    method = "gee", y = test_data$epilepsy$seizures, offset = rep(
      0,
      nrow(test_data$epilepsy)
    )
  ) {
    compute_geer_start_values(
      model_matrix = start_model_matrix(), y = y,
      family = family, weights = rep(1, nrow(test_data$epilepsy)),
      offset = offset, method = method, link = link,
      beta_start = beta_start,
      control = geer_control(), control_glm = list()
    )
  }
  binary_start_values <- function(
    beta_start = NULL, method = "gee",
    offset = rep(
      0,
      nrow(test_data$cerebrovascular)
    )
  ) {
    cer <- test_data$cerebrovascular
    compute_binary_start_values(
      model_matrix = model.matrix(ecg ~ period + treatment, data = cer),
      y = cer$ecg, family = binomial(),
      weights = rep(
        1,
        nrow(cer)
      ), offset = offset, method = method, beta_start = beta_start,
      control_glm = list(), tolerance = 1e-08, jeffreys_power = 0.5
    )
  }
  as.list(environment(), all.names = TRUE)
}, step_criterion = function() {
  fit_resp_main <- geewa_binary(
    formula = status ~ baseline + treatment +
      gender + visit + age, id = id, repeated = visit,
    data = test_data$respiratory2,
    link = "probit", orstr = "exchangeable", method = "pgee-jeffreys"
  )
  fit_resp_small <- update(fit_resp_main, formula = status ~ baseline)
  qicu_trace <- capture.output(
    qicu_search <- step_criterion(
      fit_resp_main, criterion = "QICu", trace = TRUE
    )
  )
  as.list(environment(), all.names = TRUE)
}, step_p = function() {
  fit_resp_full_exch <- geewa_binary(
    formula = status ~ baseline + treatment +
      gender + visit + age, id = id, repeated = visit,
    data = test_data$respiratory2,
    link = "probit", orstr = "exchangeable", method = "pgee-jeffreys"
  )
  fit_resp_lower <- geewa_binary(
    formula = status ~ baseline + treatment,
    id = id, repeated = visit, data = test_data$respiratory2,
    link = "probit",
    orstr = "independence", method = "pgee-jeffreys"
  )
  fit_resp_mid <- geewa_binary(
    formula = status ~ baseline + treatment +
      gender + visit, id = id, repeated = visit,
    data = test_data$respiratory2,
    link = "probit", orstr = "independence", method = "pgee-jeffreys"
  )
  fit_bin_full_indep <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id, family = binomial(link = "logit"), phi_fixed = TRUE,
    phi_value = 1,
    data = test_data$cerebrovascular, corstr = "independence",
    method = "gee"
  )
  fit_bin_lower <- geewa(
    formula = ecg ~ treatment, id = id,
    family = binomial(link = "logit"),
    phi_fixed = TRUE, phi_value = 1, data = test_data$cerebrovascular,
    corstr = "independence", method = "gee"
  )
  as.list(environment(), all.names = TRUE)
})
