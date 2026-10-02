resolve_runs_order <- function(object, order_by, env = parent.frame()) {
  natural_order <- order(object$id, object$repeated)
  if (is.character(order_by)) {
    if (length(order_by) != 1L || is.na(order_by) || !nzchar(order_by)) {
      stop(
        paste0(
          "'order_by' must be a single non-empty character value or a ",
          "numeric vector"
        ),
        call. = FALSE
      )
    }

    if (identical(order_by, "natural")) {
      return(list(
        index = natural_order,
        label = "natural cluster/repeated order",
        natural = TRUE
      ))
    }
    x_names <- colnames(object$x)
    if (identical(order_by, "fitted")) {
      key <- object$fitted.values
      label <- "fitted values"
    } else if (!is.null(x_names) && order_by %in% x_names) {
      key <- object$x[, order_by]
      label <- sprintf("model-matrix column '%s'", order_by)
    } else {
      key <- extract_runs_order_variable(object, order_by, env = env)
      label <- sprintf("variable '%s'", order_by)
    }
  } else if (is.numeric(order_by)) {
    if (length(order_by) != length(object$fitted.values)) {
      stop(
        "numeric 'order_by' must have one value per fitted observation",
        call. = FALSE
      )
    }
    key <- as.numeric(order_by)
    label <- "supplied ordering vector"
  } else {
    stop(
      "'order_by' must be a character value or a numeric vector",
      call. = FALSE
    )
  }

  if (any(!is.finite(key))) {
    stop("ordering values must all be finite and non-missing", call. = FALSE)
  }
  index <- natural_order[
    order(key[natural_order], seq_along(natural_order))
  ]

  list(index = index, label = label, natural = identical(index, natural_order))
}


extract_runs_order_variable <- function(object, name, env = parent.frame()) {
  unknown_message <- sprintf(
    paste0(
      "unknown 'order_by' value '%s': use 'natural', 'fitted', a ",
      "model-matrix column name, a variable in the data used to fit the ",
      "model, or a numeric vector"
    ),
    name
  )
  if (is.null(object$call) || is.null(object$formula)) {
    stop(unknown_message, call. = FALSE)
  }
  extended_formula <- tryCatch(
    stats::update.formula(
      object$formula,
      stats::as.formula(
        call("~", quote(.), call("+", quote(.), as.name(name)))
      )
    ),
    error = function(e) NULL
  )
  if (is.null(extended_formula)) {
    stop(unknown_message, call. = FALSE)
  }
  model_call <- object$call
  model_call$formula <- extended_formula
  eval_env <- environment(object$formula)
  if (is.null(eval_env)) {
    eval_env <- env
  }
  model_frame <- tryCatch(
    build_geer_model_frame(model_call, env = eval_env),
    error = function(e) NULL
  )
  if (is.null(model_frame) || !(name %in% names(model_frame))) {
    stop(unknown_message, call. = FALSE)
  }
  if (nrow(model_frame) != object$obs_no) {
    stop(
      sprintf(
        paste0(
          "'order_by' variable '%s' cannot be aligned: including it changes ",
          "the number of complete cases; supply an aligned numeric vector ",
          "instead"
        ),
        name
      ),
      call. = FALSE
    )
  }
  id_repeated <- extract_geer_id_repeated(model_frame, nrow(model_frame))
  ord <- order(id_repeated$id, id_repeated$repeated)
  aligned <- !is.null(object$id) && !is.null(object$repeated) &&
    isTRUE(all.equal(
      as.numeric(id_repeated$id[ord]),
      as.numeric(object$id),
      check.attributes = FALSE
    )) &&
    isTRUE(all.equal(
      as.numeric(id_repeated$repeated[ord]),
      as.numeric(object$repeated),
      check.attributes = FALSE
    ))
  if (!aligned) {
    stop(
      sprintf(
        paste0(
          "'order_by' variable '%s' cannot be aligned with the observations ",
          "stored in the fitted model; supply an aligned numeric vector ",
          "instead"
        ),
        name
      ),
      call. = FALSE
    )
  }
  key <- model_frame[[name]]
  if (inherits(key, "POSIXlt")) {
    key <- as.POSIXct(key)
  }
  if (inherits(key, c("Date", "POSIXct", "difftime"))) {
    key <- as.numeric(key)
  } else if (is.factor(key) || is.character(key) || is.logical(key)) {
    key <- as.numeric(factor(key))
  }
  if (!is.numeric(key)) {
    stop(
      sprintf(
        "'order_by' variable '%s' must be numeric, a factor, or date-like",
        name
      ),
      call. = FALSE
    )
  }
  as.numeric(key)[ord]
}


compute_runs_statistics <- function(residual_values, alternative) {
  if (!is.numeric(residual_values) || length(residual_values) < 1L) {
    stop("residuals must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(residual_values))) {
    stop("residuals must all be finite and non-missing", call. = FALSE)
  }
  if (!is.character(alternative) || length(alternative) != 1L ||
      is.na(alternative)) {
    stop("'alternative' must be a single character value", call. = FALSE)
  }

  signs <- sign(residual_values)
  retained <- signs != 0
  zero_no <- sum(!retained)
  signs <- as.integer(signs[retained])

  positive_no <- sum(signs > 0)
  negative_no <- sum(signs < 0)
  if (positive_no == 0L || negative_no == 0L) {
    stop(
      "the runs test requires at least one positive and one negative residual",
      call. = FALSE
    )
  }

  runs_no <- 1L
  if (length(signs) > 1L) {
    runs_no <- runs_no + sum(signs[-1L] != signs[-length(signs)])
  }

  n_nonzero <- positive_no + negative_no
  expected_runs <- 2 * positive_no * negative_no / n_nonzero + 1
  variance_runs <-
    2 * positive_no * negative_no *
    (2 * positive_no * negative_no - positive_no - negative_no) /
    (n_nonzero^2 * (n_nonzero - 1))

  if (!is.finite(variance_runs) || variance_runs <= 0) {
    stop(
      "the runs-test variance is not positive for this residual sign sequence",
      call. = FALSE
    )
  }

  statistic <- (runs_no - expected_runs) / sqrt(variance_runs)
  p_value <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(statistic)),
    less = stats::pnorm(statistic),
    greater = stats::pnorm(statistic, lower.tail = FALSE),
    stop(
      "'alternative' must be one of 'two.sided', 'less' or 'greater'",
      call. = FALSE
    )
  )

  list(
    signs = signs,
    retained = retained,
    runs = as.integer(runs_no),
    expected_runs = expected_runs,
    variance_runs = variance_runs,
    statistic = statistic,
    p_value = p_value,
    positive = as.integer(positive_no),
    negative = as.integer(negative_no),
    zero = as.integer(zero_no),
    nonzero = as.integer(n_nonzero)
  )
}
