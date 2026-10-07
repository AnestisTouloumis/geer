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


## Offsets written inside the model formula, as character strings such as
## "offset(log(n))". stats::update.formula() drops them whenever the right-hand
## side is replaced, so refits that rebuild the right-hand side must add them
## back. Offsets supplied through the 'offset' argument stay in the call and
## are not returned here.
geer_formula_offset_terms <- function(object) {
  offset_index <- attr(object$terms, "offset")
  if (is.null(offset_index)) {
    return(character(0))
  }
  variables <- attr(object$terms, "variables")
  vapply(
    offset_index,
    function(i) paste(deparse(variables[[i + 1L]], width.cutoff = 500L), collapse = " "),
    character(1)
  )
}


refit_geer <- function(object, formula) {
  refit_call <- stats::update(object, formula = formula, evaluate = FALSE)
  ## A fixed-length 'beta_start' does not fit a refit with added or dropped
  ## terms, so the refit always computes its own starting values.
  refit_call$beta_start <- NULL
  env <- environment(object$formula)
  if (is.null(env)) {
    env <- parent.frame()
  }
  refit <- eval(refit_call, envir = env)
  ## geewa() only warns when the solver fails and returns the last accepted
  ## iterate, which must not enter a test table as if it were the fit.
  if (isFALSE(refit$converged)) {
    stop(
      "the refit with formula '",
      paste(deparse(formula, width.cutoff = 500L), collapse = " "),
      "' did not converge; the test cannot be computed",
      call. = FALSE
    )
  }
  refit
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


## Standard errors from variances that may be unusable. A negative or
## non-finite variance (the bias-corrected covariance estimator is not
## guaranteed to be positive semi-definite in small samples) gives NA and a
## warning that names the affected entries, instead of a silent NA or NaN.
## confint() is stricter and stops, see compute_covariance_standard_errors().
standard_errors_or_na <- function(variances, context, what = "coefficient") {
  bad <- !is.finite(variances) | variances < 0
  if (any(bad)) {
    labels <- names(variances)
    if (is.null(labels)) {
      labels <- as.character(seq_along(variances))
    }
    warning(
      sprintf(
        paste0(
          "%s: the variance is negative or non-finite for %s %s, so its ",
          "standard error is set to NA"
        ),
        context,
        what,
        paste(labels[bad], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  se <- rep.int(NA_real_, length(variances))
  se[!bad] <- sqrt(variances[!bad])
  names(se) <- names(variances)
  se
}
