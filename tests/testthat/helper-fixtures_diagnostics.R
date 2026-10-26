test_diagnostics_fixtures <- list(frechet_bounds_cor = function() {
  fit_bin_exch_T2 <- geewa(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"), data = test_data$cerebrovascular,
    id = id, corstr = "exchangeable", method = "gee"
  )
  fit_bin_exch_T4 <- geewa(
    formula = status ~ baseline + treatment,
    family = binomial(link = "logit"),
    data = test_data$respiratory2, id = id, repeated = visit,
    corstr = "exchangeable",
    method = "gee"
  )
  fit_bin_unstr_T4 <- geewa(
    formula = status ~ baseline + treatment,
    family = binomial(link = "logit"),
    data = test_data$respiratory2, id = id, repeated = visit,
    corstr = "unstructured",
    method = "gee"
  )
  fit_bin_ar1_T4 <- geewa(
    formula = status ~ baseline + treatment,
    family = binomial(link = "logit"),
    data = test_data$respiratory2, id = id, repeated = visit,
    corstr = "ar1",
    method = "gee"
  )
  as.list(environment(), all.names = TRUE)
}, gof_tests = function() {
  fit_agresti <- geewa(
    formula = abnormal ~ time + diagnosis + treatment,
    data = test_data$agresti, id = id, repeated = time,
    family = binomial(link = "logit"),
    corstr = "independence", phi_fixed = TRUE, phi_value = 1
  )
  fit_resp_or <- geewa_binary(
    formula = status ~ baseline + treatment +
      gender + visit + age, data = test_data$respiratory2, id = id,
    repeated = visit,
    link = "logit", orstr = "exchangeable"
  )
  as.list(environment(), all.names = TRUE)
}, mcar_homoscedasticity_test = function() {
  make_jj_screening_data <- function() {
    withr::local_seed(9104)
    x <- matrix(rnorm(180), nrow = 60, ncol = 3)
    colnames(x) <- c("visit1", "visit2", "visit3")
    x[31:45, 3] <- NA_real_
    x[46:60, 2:3] <- NA_real_
    x
  }
  as.list(environment(), all.names = TRUE)
}, mcar_internals = function() {
  mcar_internals_matrix <- function() {
    cbind(y1 = c(1, 2.1, 2.9, 4.2, 5.1, 5.8, 7.2, 8, 9.1, 9.9), y2 = c(
      2.2,
      1.7, 3.9, 3.1, 6.3, 5.2, 7.9, 6.8, 9.6, 8.7
    ), y3 = c(
      0.5, 1.9,
      1.1, 3.4, 2.8, 4.1, 3.6, 5.9, 5.2, 6.4
    ))
  }
  mcar_fit <- function(data = test_data$epilepsy, repeated = TRUE, ...) {
    if (repeated) {
      geewa(
        formula = seizures ~ treatment + lnbaseline + lnage,
        data = data,
        id = id, repeated = visit, family = poisson(link = "log"),
        corstr = "independence", method = "gee", ...
      )
    } else {
      geewa(
        formula = seizures ~ treatment + lnbaseline + lnage,
        data = data,
        id = id, family = poisson(link = "log"), corstr = "independence",
        method = "gee", ...
      )
    }
  }
  as.list(environment(), all.names = TRUE)
}, mcar_logistic_test = function() {
  make_mcar_dropout_data <- function() {
    withr::local_seed(7401)
    clusters <- 80L
    occasions <- 5L
    id <- rep(seq_len(clusters), each = occasions)
    visit <- rep(seq_len(occasions), times = clusters)
    treatment_subject <- rep(c(0, 1), length.out = clusters)
    treatment <- rep(treatment_subject, each = occasions)
    age_subject <- seq(35, 70, length.out = clusters)
    age <- rep(age_subject, each = occasions)
    y_complete <- 2.5 + 0.5 * treatment + 0.025 * age + 0.18 * visit +
      rnorm(length(id), sd = 0.8)
    y <- y_complete
    for (i in seq_len(clusters)) {
      idx <- which(id == i)
      for (j in 2:occasions) {
        eta <- -3.2 + 0.55 * treatment[idx[j]] + 0.45 *
          y_complete[idx[j - 1L]] + 0.12 * (j - 2L)
        if (rbinom(1L, 1L, plogis(eta)) == 1L) {
          y[idx[j:occasions]] <- NA_real_
          break
        }
      }
    }
    data.frame(id, visit, treatment, age, y, y_complete)
  }
  make_direct_dropout_risk_set <- function(dat) {
    rows <- vector("list", 0L)
    counter <- 0L
    for (i in unique(dat$id)) {
      idx <- which(dat$id == i)
      idx <- idx[order(dat$visit[idx])]
      if (length(idx) < 2L) next
      for (j in 2:length(idx)) {
        if (!is.na(dat$y[idx[j - 1L]])) {
          counter <- counter + 1L
          rows[[counter]] <- data.frame(
            missing = as.integer(is.na(dat$y[idx[j]])),
            id = dat$id[idx[j]], visit = dat$visit[idx[j]],
            treatment = dat$treatment[idx[j]],
            age = dat$age[idx[j]], lag_response = dat$y[idx[j - 1L]]
          )
        }
      }
    }
    do.call(rbind, rows)
  }
  as.list(environment(), all.names = TRUE)
}, runs_test_internals = function() {
  runs_order_data <- function() {
    data <- test_data$epilepsy
    data$date <- as.Date("2020-01-01") + 7 * data$visit
    data$label <- c("first", "second", "third", "fourth")[data$visit]
    data$late <- data$visit > 2
    data$stage <- factor(data$label, levels = c(
      "first", "second", "third",
      "fourth"
    ))
    data$complex_key <- complex(real = 1, imaginary = data$visit)
    data
  }
  fit_runs_order <- function(data) {
    geewa(
      formula = seizures ~ treatment + lnbaseline + lnage, data = data,
      id = id, family = poisson(link = "log"), corstr = "independence",
      method = "gee"
    )
  }
  as.list(environment(), all.names = TRUE)
}, runs_test = function() {
  count_runs <- function(signs) {
    signs <- signs[signs != 0]
    1L + sum(signs[-1L] != signs[-length(signs)])
  }
  as.list(environment(), all.names = TRUE)
})
