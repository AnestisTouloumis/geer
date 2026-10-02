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


compute_covariance_standard_errors <- function(vcov_matrix, parm = NULL, context) {
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
