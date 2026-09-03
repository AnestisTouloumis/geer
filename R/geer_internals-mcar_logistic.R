make_mcar_covariate_formula <- function(object, formula) {
  if (is.null(formula)) {
    object_formula <- stats::formula(object)
    term_labels <- attr(object$terms, "term.labels")
    if (length(term_labels)) {
      formula <- stats::reformulate(term_labels)
    } else {
      formula <- stats::as.formula("~ 1")
    }
    environment(formula) <- environment(object_formula)
  } else {
    if (!inherits(formula, "formula")) {
      stop("'formula' must be a one-sided formula", call. = FALSE)
    }
    if (length(formula) != 2L) {
      stop(
        "'formula' must be one-sided, for example '~ treatment + age'",
        call. = FALSE
      )
    }
  }

  trms <- stats::terms(formula)
  if (!identical(attr(trms, "intercept"), 1L)) {
    stop(
      "'formula' must include an intercept; use '~ 1' when no covariates are required",
      call. = FALSE
    )
  }
  if (!is.null(attr(trms, "offset")) && length(attr(trms, "offset"))) {
    stop("offset terms are not supported in 'formula'", call. = FALSE)
  }

  formula
}


normalize_mcar_response <- function(response, object) {
  if (is.matrix(response) || is.data.frame(response)) {
    stop(
      "the regression-based MCAR diagnostic currently requires a univariate response at each occasion",
      call. = FALSE
    )
  }

  is_binomial_like <-
    identical(object$family$family, "binomial") ||
    identical(object$family$family, "quasibinomial") ||
    (identical(object$family$family, "quasi") &&
      !is.null(object$family$varfun) &&
      identical(object$family$varfun, "mu(1-mu)"))

  if (is_binomial_like && is.factor(response)) {
    if (nlevels(response) != 2L) {
      stop("a binomial factor response must have exactly two levels", call. = FALSE)
    }
    out <- as.numeric(response == levels(response)[2L])
    out[is.na(response)] <- NA_real_
    return(out)
  }
  if (is_binomial_like && is.character(response)) {
    observed <- response[!is.na(response)]
    lev <- levels(factor(observed))
    if (length(lev) != 2L) {
      stop("a binomial character response must have exactly two levels", call. = FALSE)
    }
    out <- as.numeric(response == lev[2L])
    out[is.na(response)] <- NA_real_
    return(out)
  }

  if (!is.numeric(response) && !is.logical(response)) {
    stop(
      "the response must be numeric or a two-level binary response for the regression-based MCAR diagnostic",
      call. = FALSE
    )
  }
  out <- as.numeric(response)
  if (any(!is.na(out) & !is.finite(out))) {
    stop("observed response values must be finite", call. = FALSE)
  }
  out
}


reconstruct_mcar_transition_frame <- function(object, formula, data = NULL) {
  if (is.null(data)) data <- object$data
  if (is.null(data)) {
    stop(
      "the original data are not available in the fitted object; supply them through 'data'",
      call. = FALSE
    )
  }

  object_formula <- stats::formula(object)
  formula_env <- environment(object_formula)
  if (is.null(formula_env)) formula_env <- parent.frame()
  eval_env <- new.env(parent = formula_env)
  assign(".geer_mcar_data", data, envir = eval_env)
  assign(".geer_original_formula", object_formula, envir = eval_env)

  response_call <- object$call
  response_call <- response_call[c(1L, match("data", names(response_call), 0L))]
  response_call[[1L]] <- quote(stats::model.frame)
  response_call$formula <- quote(.geer_original_formula)
  response_call$data <- quote(.geer_mcar_data)
  response_call$na.action <- quote(stats::na.pass)
  response_call$drop.unused.levels <- FALSE

  response_frame <- tryCatch(
    eval(response_call, envir = eval_env, enclos = formula_env),
    error = function(e) {
      stop(
        sprintf(
          "could not reconstruct the original response for the regression-based MCAR diagnostic: %s",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  response <- normalize_mcar_response(
    stats::model.response(response_frame),
    object
  )
  if (!length(response)) {
    stop("the original response could not be recovered", call. = FALSE)
  }
  assign(".geer_mcar_response", response, envir = eval_env)

  rhs <- formula[[2L]]
  analysis_formula <- as.call(
    list(as.name("~"), as.name(".geer_mcar_response"), rhs)
  )
  class(analysis_formula) <- "formula"
  environment(analysis_formula) <- eval_env
  assign(".geer_mcar_formula", analysis_formula, envir = eval_env)

  mcall <- object$call
  keep <- c("data", "subset", "id", "repeated")
  mcall <- mcall[c(1L, match(keep, names(mcall), 0L))]
  mcall[[1L]] <- quote(stats::model.frame)
  mcall$formula <- quote(.geer_mcar_formula)
  mcall$data <- quote(.geer_mcar_data)
  mcall$na.action <- quote(stats::na.pass)
  mcall$drop.unused.levels <- TRUE

  model_frame <- tryCatch(
    eval(mcall, envir = eval_env, enclos = formula_env),
    error = function(e) {
      stop(
        sprintf(
          "could not construct the risk-set data for the regression-based MCAR diagnostic: %s",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  y <- stats::model.response(model_frame)
  if (length(y) != nrow(model_frame)) {
    stop("the response could not be reconstructed correctly", call. = FALSE)
  }

  id_raw <- stats::model.extract(model_frame, "id")
  if (is.null(id_raw)) {
    stop("'id' could not be recovered from the fitted model", call. = FALSE)
  }
  if (anyNA(id_raw)) {
    stop("'id' cannot contain missing values in the missingness model", call. = FALSE)
  }
  id <- as.numeric(factor(id_raw))

  repeated_raw <- stats::model.extract(model_frame, "repeated")
  if (is.null(repeated_raw)) {
    repeated <- stats::ave(id, id, FUN = seq_along)
  } else {
    if (anyNA(repeated_raw)) {
      stop("'repeated' cannot contain missing values in the missingness model", call. = FALSE)
    }
    repeated <- as.numeric(factor(repeated_raw))
  }

  if (length(y) != length(id) || length(y) != length(repeated)) {
    stop("response, 'id', and 'repeated' must have the same length", call. = FALSE)
  }
  if (any(duplicated(data.frame(id = id, repeated = repeated)))) {
    stop("'repeated' must identify unique measurements within each 'id'", call. = FALSE)
  }

  trms <- stats::terms(analysis_formula, data = model_frame)
  design <- stats::model.matrix(trms, model_frame)
  intercept <- which(colnames(design) == "(Intercept)")
  if (length(intercept) != 1L) {
    stop("the covariate formula must contain exactly one intercept", call. = FALSE)
  }
  covariate_index <- setdiff(seq_len(ncol(design)), intercept)
  covariates <- design[, covariate_index, drop = FALSE]
  if (anyNA(covariates) || any(!is.finite(covariates))) {
    stop(
      "covariates in the missingness model must be fully observed and finite",
      call. = FALSE
    )
  }

  ord <- order(id, repeated)
  y <- as.numeric(y[ord])
  id <- id[ord]
  repeated <- repeated[ord]
  covariates <- covariates[ord, , drop = FALSE]

  current_index <- integer(0)
  previous_index <- integer(0)
  intermittent <- FALSE
  id_split <- split(seq_along(id), id)
  for (idx in id_split) {
    if (length(idx) < 2L) next
    missing_i <- is.na(y[idx])
    if (any(missing_i)) {
      first_missing <- which(missing_i)[1L]
      if (first_missing < length(idx) && any(!missing_i[(first_missing + 1L):length(idx)])) {
        intermittent <- TRUE
      }
    }
    for (k in 2:length(idx)) {
      if (!is.na(y[idx[k - 1L]])) {
        previous_index <- c(previous_index, idx[k - 1L])
        current_index <- c(current_index, idx[k])
      }
    }
  }

  if (!length(current_index)) {
    stop(
      "no transition risk set could be constructed; at least one cluster must have two occasions with the previous response observed",
      call. = FALSE
    )
  }

  list(
    missing = as.integer(is.na(y[current_index])),
    lag_response = y[previous_index],
    id = id[current_index],
    repeated = repeated[current_index],
    occasion = repeated[current_index],
    covariates = covariates[current_index, , drop = FALSE],
    covariate_names = colnames(covariates),
    formula = formula,
    intermittent = intermittent
  )
}


select_mcar_covariates <- function(covariates, occasion, lag_response) {
  occasion_factor <- factor(occasion)
  nuisance <- if (nlevels(occasion_factor) > 1L) {
    stats::model.matrix(~ occasion_factor)
  } else {
    matrix(1, nrow = length(occasion), ncol = 1L)
  }
  base <- nuisance
  rank_base <- qr(base)$rank

  keep <- logical(ncol(covariates))
  if (ncol(covariates)) {
    for (j in seq_len(ncol(covariates))) {
      candidate <- cbind(base, covariates[, j, drop = FALSE])
      rank_candidate <- qr(candidate)$rank
      if (rank_candidate > rank_base) {
        keep[j] <- TRUE
        base <- candidate
        rank_base <- rank_candidate
      }
    }
  }

  lag_candidate <- cbind(base, lag_response)
  if (qr(lag_candidate)$rank <= rank_base) {
    stop(
      "the previous response has no estimable effect after adjustment for occasion and the supplied covariates",
      call. = FALSE
    )
  }

  list(
    covariates = covariates[, keep, drop = FALSE],
    kept_names = colnames(covariates)[keep],
    dropped_names = colnames(covariates)[!keep],
    occasion_factor = occasion_factor
  )
}


fit_mcar_binary_model <- function(formula, analysis_data, orstr, control) {
  mcar_id <- analysis_data[["mcar_id"]]
  mcar_repeated <- analysis_data[["mcar_repeated"]]

  fit <- geewa_binary(
    formula = formula,
    link = "logit",
    data = analysis_data,
    id = mcar_id,
    repeated = mcar_repeated,
    control = control,
    orstr = orstr,
    method = "gee"
  )
  if (!isTRUE(fit$converged)) {
    stop("the binary GEE missingness model did not converge", call. = FALSE)
  }
  fit
}


mcar_nested_geer_test <- function(object0,
                                  object1,
                                  test,
                                  cov_type,
                                  pmethod) {
  value <- switch(
    test,
    wald = wald_test(object0, object1, cov_type),
    score = score_test(object0, object1, cov_type),
    `working-wald` = working_wald_test(
      object0, object1, cov_type, pmethod
    ),
    `working-score` = working_score_test(
      object0, object1, cov_type, pmethod
    ),
    `working-lrt` = working_lrt_test(
      object0, object1, cov_type, pmethod
    )
  )
  list(
    statistic = value$test_stat,
    df = as.numeric(value$test_df),
    p_value = value$test_p
  )
}


mcar_zero_test <- function() {
  list(statistic = 0, df = 0, p_value = 1)
}


mcar_coefficient_table <- function(fit, covariance, covariate_map) {
  beta <- stats::coef(fit)
  standard_error <- sqrt(diag(covariance))
  z_value <- as.numeric(beta) / standard_error
  display_names <- names(beta)
  if (length(covariate_map)) {
    for (safe in names(covariate_map)) {
      display_names[display_names == safe] <- covariate_map[[safe]]
    }
  }
  display_names[display_names == "mcar_lag_response"] <- "lag_response"
  display_names <- sub("^mcar_occasion", "occasion", display_names)

  data.frame(
    term = display_names,
    estimate = as.numeric(beta),
    std.error = standard_error,
    statistic = z_value,
    p.value = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE),
    row.names = NULL,
    check.names = FALSE
  )
}
