resolve_runs_order <- function(object, order_by) {
  n_obs <- length(object$fitted.values)
  natural_order <- order(
    object$id,
    object$repeated,
    seq_len(n_obs)
  )

  if (is.character(order_by)) {
    if (length(order_by) != 1L || is.na(order_by) || !nzchar(order_by)) {
      stop("'order_by' must be a single non-empty character value or a numeric vector", call. = FALSE)
    }

    if (identical(order_by, "natural")) {
      return(list(
        index = natural_order,
        label = "natural cluster/repeated order"
      ))
    }

    if (identical(order_by, "fitted")) {
      key <- object$fitted.values
      label <- "fitted values"
    } else {
      x_names <- colnames(object$x)
      if (is.null(x_names) || !(order_by %in% x_names)) {
        stop(
          sprintf(
            "unknown 'order_by' value '%s'; use 'natural', 'fitted', a model-matrix column name, or a numeric vector",
            order_by
          ),
          call. = FALSE
        )
      }
      key <- object$x[, order_by]
      label <- sprintf("model-matrix column '%s'", order_by)
    }
  } else if (is.numeric(order_by)) {
    if (length(order_by) != n_obs) {
      stop("numeric 'order_by' must have one value per fitted observation", call. = FALSE)
    }
    key <- as.numeric(order_by)
    label <- "supplied ordering vector"
  } else {
    stop("'order_by' must be a character value or a numeric vector", call. = FALSE)
  }

  if (anyNA(key) || any(!is.finite(key))) {
    stop("ordering values must all be finite and non-missing", call. = FALSE)
  }

  ordered_natural <- natural_order
  tie_break <- seq_along(ordered_natural)
  index <- ordered_natural[order(key[ordered_natural], tie_break)]

  list(index = index, label = label)
}


compute_runs_statistics <- function(residual_values) {
  if (!is.numeric(residual_values) || length(residual_values) < 1L) {
    stop("residuals must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(residual_values) || any(!is.finite(residual_values))) {
    stop("residuals must all be finite and non-missing", call. = FALSE)
  }

  signs <- sign(residual_values)
  zero_no <- sum(signs == 0)
  signs <- signs[signs != 0]

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
  p_value <- 2 * stats::pnorm(-abs(statistic))

  list(
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
