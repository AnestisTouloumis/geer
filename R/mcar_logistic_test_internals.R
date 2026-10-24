build_mcar_covariate_formula <- function(object, formula) {
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
      paste0(
        "'formula' must include an intercept; use '~ 1' when no covariates ",
        "are required"
      ),
      call. = FALSE
    )
  }
  if (!is.null(attr(trms, "offset")) && length(attr(trms, "offset"))) {
    stop("offset terms are not supported in 'formula'", call. = FALSE)
  }

  formula
}


normalize_mcar_response <- function(response, object) {
  check_mcar_univariate(response)
  is_binomial_like <- is_binomial_like_family(object$family)

  if (is_binomial_like && is.factor(response)) {
    if (nlevels(response) != 2L) {
      stop(
        paste0(
          "for binomial-type models, a factor response must have exactly two ",
          "levels"
        ),
        call. = FALSE
      )
    }
    out <- as.numeric(response == levels(response)[2L])
    out[is.na(response)] <- NA_real_
    return(out)
  }
  if (is_binomial_like && is.character(response)) {
    observed <- response[!is.na(response)]
    lev <- levels(factor(observed))
    if (length(lev) != 2L) {
      stop(
        paste0(
          "for binomial-type models, a character response must have exactly ",
          "two distinct values"
        ),
        call. = FALSE
      )
    }
    out <- as.numeric(response == lev[2L])
    out[is.na(response)] <- NA_real_
    return(out)
  }

  if (!is.numeric(response) && !is.logical(response)) {
    stop(
      "the response must be numeric or binary with two levels",
      call. = FALSE
    )
  }
  out <- as.numeric(response)
  if (any(!is.na(out) & !is.finite(out))) {
    stop("observed response values must be finite", call. = FALSE)
  }
  out
}


build_mcar_transition_frame <- function(object, formula, data = NULL) {
  data <- get_mcar_data(object, data)

  object_formula <- stats::formula(object)
  formula_env <- environment(object_formula)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }
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

  response_frame <- build_mcar_model_frame(
    response_call,
    eval_env,
    formula_env,
    "reconstruct the response used by the fit",
    context = "mcar_logistic_test"
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

  model_frame <- build_mcar_model_frame(
    mcall,
    eval_env,
    formula_env,
    "construct the risk-set data",
    context = "mcar_logistic_test"
  )

  y <- stats::model.response(model_frame)
  if (length(y) != nrow(model_frame)) {
    stop("the response could not be reconstructed correctly", call. = FALSE)
  }

  id_repeated <- extract_geer_id_repeated(model_frame, length(y))
  id <- id_repeated$id
  repeated <- id_repeated$repeated
  repeated_values <- stats::model.extract(model_frame, "repeated")
  if (is.null(repeated_values)) {
    repeated_values <- repeated
  }

  trms <- stats::terms(analysis_formula, data = model_frame)
  design <- stats::model.matrix(trms, model_frame)
  intercept <- which(colnames(design) == "(Intercept)")
  if (length(intercept) != 1L) {
    stop(
      "the covariate formula must contain exactly one intercept",
      call. = FALSE
    )
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
  repeated_values <- repeated_values[ord]
  covariates <- covariates[ord, , drop = FALSE]

  current_index <- integer(0)
  previous_index <- integer(0)
  intermittent <- FALSE
  id_split <- split(seq_along(id), id)
  incomplete_clusters <- any(lengths(id_split) < max(repeated))
  for (idx in id_split) {
    if (length(idx) < 2L) {
      next
    }
    missing_i <- is.na(y[idx])
    if (any(missing_i)) {
      first_missing <- which(missing_i)[1L]
      if (first_missing < length(idx) &&
          any(!missing_i[(first_missing + 1L):length(idx)])) {
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
      paste0(
        "no transition risk set could be constructed; at least one cluster ",
        "must have two occasions with the previous response observed"
      ),
      call. = FALSE
    )
  }

  list(
    missing = as.integer(is.na(y[current_index])),
    lag_response = y[previous_index],
    id = id[current_index],
    repeated = repeated[current_index],
    occasion = repeated[current_index],
    occasion_label = as.character(repeated_values[current_index]),
    covariates = covariates[current_index, , drop = FALSE],
    intermittent = intermittent,
    incomplete_clusters = incomplete_clusters
  )
}


normalize_mcar_risk_set <- function(prepared) {
  events <- tapply(prepared$missing, prepared$occasion, sum)
  totals <- tapply(prepared$missing, prepared$occasion, length)
  degenerate <- as.numeric(names(events))[events == 0L | events == totals]
  if (!length(degenerate)) {
    return(prepared)
  }
  keep <- !(prepared$occasion %in% degenerate)
  if (!any(keep)) {
    stop(
      paste0(
        "mcar_logistic_test failed: every occasion has either no dropouts ",
        "or no observed responses in the risk set; the missingness model ",
        "cannot be fitted"
      ),
      call. = FALSE
    )
  }
  omitted <- !keep
  labels <- unique(
    prepared$occasion_label[omitted][order(prepared$occasion[omitted])]
  )
  warning(
    sprintf(
      paste0(
        "occasion(s) %s were omitted from the risk set because they have no ",
        "dropouts or no observed responses, so their occasion effects cannot ",
        "be estimated"
      ),
      paste(labels, collapse = ", ")
    ),
    call. = FALSE
  )
  for (component in c("missing", "lag_response", "id", "repeated",
                      "occasion", "occasion_label")) {
    prepared[[component]] <- prepared[[component]][keep]
  }
  prepared$covariates <- prepared$covariates[keep, , drop = FALSE]
  prepared
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
      paste0(
        "mcar_logistic_test failed: the previous response has no estimable ",
        "effect after adjustment for occasion and the supplied covariates"
      ),
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
    stop(
      paste0(
        "mcar_logistic_test failed: the binary GEE fit of the missingness ",
        "model did not converge"
      ),
      call. = FALSE
    )
  }
  fit
}


compute_mcar_nested_test <- function(object0,
                                     object1,
                                     test,
                                     cov_type,
                                     pmethod) {
  value <- compute_nested_test(
    object0,
    object1,
    test,
    cov_type,
    pmethod,
    "mcar_logistic_test"
  )
  list(
    statistic = value$test_stat,
    df = as.numeric(value$test_df),
    p_value = value$test_p
  )
}


build_mcar_null_formula <- function(term_labels) {
  if (length(term_labels) == 0L) {
    return(stats::as.formula("mcar_missing ~ 1"))
  }
  stats::reformulate(term_labels, response = "mcar_missing")
}


build_mcar_zero_test <- function() {
  list(statistic = 0, df = 0, p_value = 1)
}


build_mcar_coefficient_table <- function(fit, covariance, covariate_map) {
  beta <- stats::coef(fit)
  standard_error <- compute_se_or_stop(
    diag(covariance),
    context = "mcar_logistic_test"
  )
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
