geer_test_choices <- c(
  "wald", "score", "working-wald", "working-score", "working-lrt"
)

geer_cov_type_choices <- c(
  "robust", "bias-corrected", "df-adjusted", "naive"
)

geer_pmethod_choices <- c(
  "rao-scott", "satterthwaite"
)

geer_direction_choices <- c(
  "backward", "forward", "both"
)

geer_integer_tol <- sqrt(.Machine$double.eps)


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


check_class <- function(x, class_name, name) {
  if (!inherits(x, class_name)) {
    stop(sprintf("'%s' must be of '%s' class", name, class_name), call. = FALSE)
  }
  x
}

check_geer_object <- function(x, name = "object") {
  check_class(x, class_name = "geer", name = name)
}


check_summary_geer_object <- function(x, name = "x") {
  check_class(x, class_name = "summary.geer", name = name)
}


check_nonnegative_integerish <- function(x, name) {
  check_single_numeric(x, name)
  if (x < 0 || abs(x - round(x)) > geer_integer_tol) {
    stop(sprintf("'%s' must be a single non-negative integer", name), call. = FALSE)
  }
  as.integer(round(x))
}


check_probability_open <- function(x, name) {
  check_probability(x, name, open = TRUE)
  invisible(x)
}


match_choice_with_default <- function(x,
                                      choices,
                                      name,
                                      default = choices[1L]) {
  if (missing(x) || is.null(x)) {
    return(default)
  }

  check_choice(x, choices, name)
  x
}


match_test_type <- function(test) {
  match_choice_with_default(
    x = test,
    choices = geer_test_choices,
    name = "test"
  )
}


match_cov_type <- function(cov_type) {
  match_choice_with_default(
    x = cov_type,
    choices = geer_cov_type_choices,
    name = "cov_type"
  )
}


match_pmethod <- function(pmethod) {
  match_choice_with_default(
    x = pmethod,
    choices = geer_pmethod_choices,
    name = "pmethod"
  )
}


match_direction_type <- function(direction) {
  match_choice_with_default(
    x = direction,
    choices = geer_direction_choices,
    name = "direction"
  )
}


is_working_test <- function(test) {
  test %in% c("working-wald", "working-score", "working-lrt")
}


check_working_lrt_allowed <- function(object,
                                      test = "working-lrt",
                                      name = "object") {
  object <- check_geer_object(object, name = name)
  test <- match_test_type(test)
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
  invisible(list(test = test, cov_type = cov_type))
}


normalize_geer_test_options <- function(test,
                                   cov_type,
                                   pmethod = NULL,
                                   object = NULL) {
  test <- match_test_type(test)
  cov_type <- match_cov_type(cov_type)
  if (is_working_test(test)) {
    pmethod <- match_pmethod(pmethod %||% geer_pmethod_choices[1L])
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


validate_step_count <- function(steps) {
  check_nonnegative_integerish(steps, "steps")
}
