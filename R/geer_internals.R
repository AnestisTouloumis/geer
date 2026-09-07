format_percent <- function(probs, digits = 2, check = FALSE) {
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
    stop("'check' must be a single logical value", call. = FALSE)
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


is_positive_scalar <- function(x) {
  is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x > 0
}


is_positive_integer_scalar <- function(x, tol = geer_integer_tol) {
  is_positive_scalar(x) && abs(x - round(x)) <= tol
}


check_single_numeric <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric value", name), call. = FALSE)
  }
  invisible(x)
}


check_probability_open <- function(x, name) {
  check_single_numeric(x, name)
  if (x <= 0 || x >= 1) {
    stop(sprintf("'%s' must be strictly between 0 and 1", name), call. = FALSE)
  }
  invisible(x)
}


check_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be a single character value", name), call. = FALSE)
  }
  if (!(x %in% choices)) {
    stop(
      sprintf("'%s' must be one of: %s", name, paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(x)
}


format_test_label <- function(test) {
  out <- switch(
    test,
    wald = "Wald",
    score = "Score",
    `working-wald` = "Modified Working Wald",
    `working-score` = "Modified Working Score",
    `working-lrt` = "Modified Working LRT"
  )
  if (is.null(out)) {
    stop("invalid 'test' value", call. = FALSE)
  }
  out
}


refit_geer <- function(object, formula) {
  refit_call <- stats::update(object, formula = formula, evaluate = FALSE)
  env <- environment(object$formula)
  if (is.null(env)) {
    env <- parent.frame()
  }
  eval(refit_call, envir = env)
}


check_unused_dots <- function(dots, context) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep.int("", length(dots))
  }
  labels <- rep.int("<unnamed>", length(dots))
  named <- nzchar(dot_names)
  labels[named] <- sprintf("'%s'", dot_names[named])
  stop(
    sprintf(
      "%s does not use the argument(s): %s",
      context,
      paste(labels, collapse = ", ")
    ),
    call. = FALSE
  )
}


covariance_standard_errors <- function(vcov_matrix, parm = NULL, context) {
  variances <- diag(vcov_matrix)
  if (!is.null(parm)) {
    variances <- variances[parm]
  }
  labels <- names(variances)
  if (is.null(labels)) {
    labels <- as.character(seq_along(variances))
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
