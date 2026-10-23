build_percent_labels <- function(probs, digits = 2, check = FALSE) {
  if (!is.numeric(probs)) {
    stop("'probs' must be numeric", call. = FALSE)
  }
  if (anyNA(probs) || any(!is.finite(probs))) {
    stop("'probs' must contain only finite numeric values", call. = FALSE)
  }
  if (!is_positive_integer_scalar(digits)) {
    stop("'digits' must be a positive integer", call. = FALSE)
  }
  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("'check' must be TRUE or FALSE", call. = FALSE)
  }
  if (check && any(probs < 0 | probs > 1)) {
    warning("some probabilities are outside [0, 1]", call. = FALSE)
  }
  formatted <- formatC(
    100 * probs,
    digits = digits,
    format = "fg",
    flag = "",
    drop0trailing = FALSE
  )
  paste0(formatted, "%")
}


build_test_label <- function(test) {
  out <- switch(
    test,
    wald = "Wald",
    score = "Score",
    `working-wald` = "Modified Working Wald",
    `working-score` = "Modified Working Score",
    `working-lrt` = "Modified Working LRT"
  )
  if (is.null(out)) {
    stop(
      paste0(
        "'test' must be one of: wald, score, working-wald, working-score, ",
        "working-lrt"
      ),
      call. = FALSE
    )
  }
  out
}


extract_formula_offset_terms <- function(object) {
  offset_index <- attr(object$terms, "offset")
  if (is.null(offset_index)) {
    return(character(0))
  }
  variables <- attr(object$terms, "variables")
  vapply(
    offset_index,
    function(i) {
      paste(
        deparse(variables[[i + 1L]], width.cutoff = 500L),
        collapse = " "
      )
    },
    character(1)
  )
}


refit_geer <- function(object, formula, context) {
  refit_call <- stats::update(object, formula. = formula, evaluate = FALSE)
  refit_call$beta_start <- NULL
  env <- environment(object$formula)
  if (is.null(env)) {
    env <- parent.frame()
  }
  evaluated_call <- refit_call
  frame_rows <- object$frame_rows
  if (is.data.frame(object$data) && !is.null(refit_call$data) &&
      is.numeric(frame_rows) && length(frame_rows) == NROW(object$y)) {
    keep <- logical(nrow(object$data))
    keep[frame_rows] <- TRUE
    evaluated_call$data <- object$data
    evaluated_call$subset <- keep
  }
  refit <- eval(evaluated_call, envir = env)
  refit$call <- refit_call
  if (isFALSE(refit$converged)) {
    stop(
      sprintf(
        "%s failed: the fit of '%s' did not converge",
        context,
        paste(deparse(stats::formula(refit), width.cutoff = 500L),
              collapse = " ")
      ),
      call. = FALSE
    )
  }
  refit
}


compute_se_or_stop <- function(variances, context, what = NULL) {
  labels <- names(variances)
  if (is.null(labels)) {
    labels <- as.character(seq_along(variances))
  }
  if (!is.null(what)) {
    labels <- paste(what, labels)
  }
  invalid <- !is.finite(variances)
  if (any(invalid)) {
    stop(
      sprintf(
        "%s failed: non-finite variance for %s",
        context,
        paste(labels[invalid], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  negative <- variances < 0
  if (any(negative)) {
    stop(
      sprintf(
        paste0(
          "%s failed: negative variance for %s; the covariance estimator is ",
          "not positive semi-definite, so no standard error is available"
        ),
        context,
        paste(labels[negative], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sqrt(variances)
}


get_cluster_labels <- function(object) {
  codes <- as.numeric(names(split(seq_along(object$id), object$id)))
  if (is.null(object$cluster_labels)) {
    return(as.character(codes))
  }
  object$cluster_labels[codes]
}


get_observation_labels <- function(object) {
  cluster_labels <- if (is.null(object$cluster_labels)) {
    as.character(object$id)
  } else {
    object$cluster_labels[object$id]
  }
  repeated_labels <- if (is.null(object$repeated_labels)) {
    as.character(object$repeated)
  } else {
    object$repeated_labels[object$repeated]
  }
  paste0(cluster_labels, "(", repeated_labels, ")")
}


is_balanced_design <- function(object) {
  positions <- seq_len(max(as.integer(object$repeated)))
  cluster_indices <- split(seq_along(object$id), object$id)
  all(vapply(
    cluster_indices,
    function(indices) {
      identical(sort(as.integer(object$repeated[indices])), positions)
    },
    logical(1)
  ))
}


build_fit_data_frames <- function(object, env) {
  if (is.null(object$call) || is.null(object$formula)) {
    return(NULL)
  }
  fit_env <- environment(object$formula)
  if (is.null(fit_env)) {
    fit_env <- env
  }
  model_call <- object$call
  model_call$formula <- object$formula
  fit_frame <- tryCatch(
    build_geer_model_frame(model_call, env = fit_env),
    error = function(e) NULL
  )
  if (is.null(fit_frame)) {
    return(NULL)
  }
  data <- NULL
  if (!is.null(model_call$data)) {
    data <- tryCatch(
      eval(model_call$data, envir = fit_env),
      error = function(e) NULL
    )
    if (is.null(data)) {
      return(NULL)
    }
  }
  list(fit_frame = fit_frame, data = data, fit_env = fit_env)
}


evaluate_fit_data_variables <- function(variables, data, enclos, rows) {
  variable_formula <- stats::as.formula(
    call("~", Reduce(function(x, y) call("+", x, y), variables)),
    env = enclos
  )
  variable_frame <- tryCatch(
    stats::model.frame(
      variable_formula,
      data = data,
      na.action = stats::na.pass
    ),
    error = function(e) NULL
  )
  if (is.null(variable_frame) || ncol(variable_frame) != length(variables)) {
    return(NULL)
  }
  index <- match(rows, rownames(variable_frame))
  if (anyNA(index)) {
    return(NULL)
  }
  lapply(variable_frame, function(x) x[index])
}


get_fit_frame_order <- function(object, fit_frame) {
  if (nrow(fit_frame) != object$obs_no || is.null(object$id) ||
      is.null(object$repeated)) {
    return(NULL)
  }
  id_repeated <- extract_geer_id_repeated(fit_frame, nrow(fit_frame))
  ord <- order(id_repeated$id, id_repeated$repeated)
  aligned <- isTRUE(all.equal(
    as.numeric(id_repeated$id[ord]),
    as.numeric(object$id),
    check.attributes = FALSE
  )) &&
    isTRUE(all.equal(
      as.numeric(id_repeated$repeated[ord]),
      as.numeric(object$repeated),
      check.attributes = FALSE
    ))
  if (aligned) ord else NULL
}


extract_fit_data_variables <- function(object, variables, enclos, env) {
  frames <- build_fit_data_frames(object, env)
  if (is.null(frames)) {
    return(list(status = "unavailable", values = NULL))
  }
  if (is.null(enclos)) {
    enclos <- frames$fit_env
  }
  values <- evaluate_fit_data_variables(
    variables,
    frames$data,
    enclos,
    rownames(frames$fit_frame)
  )
  if (is.null(values)) {
    return(list(status = "unavailable", values = NULL))
  }
  if (any(vapply(values, anyNA, logical(1)))) {
    return(list(status = "incomplete", values = NULL))
  }
  ord <- get_fit_frame_order(object, frames$fit_frame)
  if (is.null(ord)) {
    return(list(status = "unaligned", values = NULL))
  }
  list(status = "ok", values = lapply(values, function(x) x[ord]))
}


get_cached_value <- function(object, key, compute) {
  cache <- object$cache
  if (is.environment(cache)) {
    stored <- cache[[key]]
    if (!is.null(stored) &&
        identical(stored$coefficients, object$coefficients)) {
      return(stored$value)
    }
  }
  value <- compute()
  if (is.environment(cache)) {
    assign(
      key,
      list(coefficients = object$coefficients, value = value),
      envir = cache
    )
  }
  value
}


is_binary_response <- function(y) {
  !anyNA(y) && all(y == 0 | y == 1)
}


is_degenerate_probability <- function(mu) {
  !is.finite(mu) | mu <= 0 | mu >= 1
}


get_pass_iterations <- function(pass, control) {
  if (pass$one_step) {
    return(list(1L, 1L, 1))
  }
  list(control$maxiter, control$step_maxiter, control$step_multiplier)
}
