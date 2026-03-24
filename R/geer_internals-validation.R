## internal validation helpers for geer user-facing functions

.geer_test_choices <- c(
  "wald", "score", "working-wald", "working-score", "working-lrt"
)

.geer_cov_type_choices <- c(
  "robust", "bias-corrected", "df-adjusted", "naive"
)

.geer_pmethod_choices <- c(
  "rao-scott", "satterthwaite"
)

.geer_direction_choices <- c(
  "backward", "forward", "both"
)


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


check_geer_object <- function(x, arg = "object") {
  if (!inherits(x, "geer")) {
    stop(sprintf("'%s' must be of 'geer' class", arg), call. = FALSE)
  }
  x
}


check_summary_geer_object <- function(x, arg = "x") {
  if (!inherits(x, "summary.geer")) {
    stop(sprintf("'%s' must be of 'summary.geer' class", arg), call. = FALSE)
  }
  x
}


check_nonnegative_integerish <- function(x, name) {
  check_single_numeric(x, name)
  if (x < 0 || abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    stop(sprintf("'%s' must be a single non-negative integer", name), call. = FALSE)
  }
  as.integer(round(x))
}


check_probability_open <- function(x, name) {
  check_probability(x, name, open = TRUE)
  invisible(x)
}


check_probability_closed <- function(x, name) {
  check_probability(x, name, open = FALSE)
  invisible(x)
}


match_test_type <- function(test) {
  if (missing(test) || is.null(test)) {
    test <- .geer_test_choices[1L]
  }
  check_choice(test, .geer_test_choices, "test")
  test
}


match_cov_type <- function(cov_type) {
  if (missing(cov_type) || is.null(cov_type)) {
    cov_type <- .geer_cov_type_choices[1L]
  }
  check_choice(cov_type, .geer_cov_type_choices, "cov_type")
  cov_type
}


match_pmethod <- function(pmethod) {
  if (missing(pmethod) || is.null(pmethod)) {
    pmethod <- .geer_pmethod_choices[1L]
  }
  check_choice(pmethod, .geer_pmethod_choices, "pmethod")
  pmethod
}


match_direction_type <- function(direction) {
  if (missing(direction) || is.null(direction)) {
    direction <- .geer_direction_choices[1L]
  }
  check_choice(direction, .geer_direction_choices, "direction")
  direction
}


is_working_test <- function(test) {
  test %in% c("working-wald", "working-score", "working-lrt")
}


check_working_lrt_allowed <- function(object, test = "working-lrt", arg = "object") {
  object <- check_geer_object(object, arg)
  if (identical(test, "working-lrt") &&
      !identical(object$association_structure, "independence")) {
    stop(
      "the modified working LRT can only be applied to an independence working model",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


check_test_cov_type_compatibility <- function(test, cov_type) {
  test <- match_test_type(test)
  cov_type <- match_cov_type(cov_type)
  ## all currently supported test/cov_type combinations are allowed
  invisible(TRUE)
}


normalize_test_options <- function(test,
                                   cov_type,
                                   pmethod = NULL,
                                   object = NULL) {
  test <- match_test_type(test)
  cov_type <- match_cov_type(cov_type)
  if (is_working_test(test)) {
    pmethod <- match_pmethod(pmethod %||% .geer_pmethod_choices[1L])
  } else {
    pmethod <- NULL
  }
  check_test_cov_type_compatibility(test, cov_type)
  if (identical(test, "working-lrt") && !is.null(object)) {
    check_working_lrt_allowed(object, test = test)
  }
  list(
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
}


validate_step_thresholds <- function(p_enter, p_remove) {
  check_probability_open(p_enter, "p_enter")
  check_probability_open(p_remove, "p_remove")
  list(
    p_enter = p_enter,
    p_remove = p_remove
  )
}


validate_step_direction <- function(direction) {
  match_direction_type(direction)
}


validate_step_count <- function(steps) {
  check_nonnegative_integerish(steps, "steps")
}
