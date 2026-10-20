geer_method_choices <- c(
  "gee",
  "brgee-naive", "brgee-robust", "brgee-empirical",
  "bcgee-naive", "bcgee-robust", "bcgee-empirical",
  "pgee-jeffreys", "opgee-jeffreys", "hpgee-jeffreys"
)

geer_corstr_choices <- c(
  "independence", "exchangeable", "ar1", "toeplitz",
  "m-dependent", "unstructured", "fixed"
)

geer_orstr_choices <- c(
  "independence", "exchangeable", "unstructured", "fixed"
)

geer_family_choices <- c(
  "gaussian", "poisson", "binomial", "Gamma", "inverse.gaussian",
  "quasi", "quasibinomial", "quasipoisson"
)

geer_link_choices <- c(
  "logit", "probit", "cauchit", "cloglog", "identity", "log",
  "sqrt", "1/mu^2", "inverse"
)

geer_bcgee_methods <- c(
  "bcgee-naive", "bcgee-robust", "bcgee-empirical"
)

geer_onestep_methods <- c(
  geer_bcgee_methods, "opgee-jeffreys", "hpgee-jeffreys"
)

geer_test_choices <- c(
  "wald", "score", "working-wald", "working-score", "working-lrt"
)

geer_small_sample_cov_types <- c(
  "kauermann-carroll", "mancl-derouen", "ford-westgate", "fay-graubard",
  "fan-zhang-zhang", "pan", "gosho-sato-takeuchi", "wang-long",
  "westgate-burchett"
)

geer_pooled_cov_types <- c(
  "pan", "gosho-sato-takeuchi", "wang-long", "westgate-burchett"
)

geer_cov_type_choices <- c(
  "bias-corrected", "robust", "df-adjusted", "jackknife", "naive",
  geer_small_sample_cov_types
)

geer_pmethod_choices <- c(
  "rao-scott", "satterthwaite"
)

geer_direction_choices <- c(
  "backward", "forward", "both"
)

geer_step_criteria <- c(
  "QICu", "QIC", "QICC", "AGPC", "SGPC"
)

geer_integer_tol <- sqrt(.Machine$double.eps)


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


check_class <- function(x, class_name, name) {
  if (!inherits(x, class_name)) {
    stop(sprintf("'%s' must be a '%s' object", name, class_name), call. = FALSE)
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
    stop(
      sprintf("'%s' must be a non-negative integer", name),
      call. = FALSE
    )
  }
  as.integer(round(x))
}


check_integer_at_least <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > .Machine$integer.max ||
      abs(x - round(x)) > geer_integer_tol) {
    message_text <- if (lower == 1L) {
      "'%s' must be a positive integer"
    } else {
      paste0("'%s' must be an integer greater than or equal to ", lower)
    }
    stop(sprintf(message_text, name), call. = FALSE)
  }
  as.integer(round(x))
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


match_cov_type <- function(cov_type, default = geer_cov_type_choices[1L]) {
  match_choice_with_default(
    x = cov_type,
    choices = geer_cov_type_choices,
    name = "cov_type",
    default = default
  )
}


match_pmethod <- function(pmethod) {
  match_choice_with_default(
    x = pmethod,
    choices = geer_pmethod_choices,
    name = "pmethod"
  )
}


match_step_criterion <- function(criterion) {
  if (is.character(criterion) && length(criterion) == 1L &&
      !is.na(criterion) && !(criterion %in% geer_step_criteria) &&
      criterion %in% geer_criteria_columns) {
    stop(
      sprintf(
        paste0(
          "'criterion' must be one of: %s; %s is meant for selecting the ",
          "working association structure, not the mean model"
        ),
        paste(geer_step_criteria, collapse = ", "),
        criterion
      ),
      call. = FALSE
    )
  }
  match_choice_with_default(
    x = criterion,
    choices = geer_step_criteria,
    name = "criterion"
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
      "the modified working LRT requires an independence working structure",
      call. = FALSE
    )
  }
  invisible(TRUE)
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
  if (identical(test, "working-lrt") && !is.null(object)) {
    check_working_lrt_allowed(object, test = test)
  }
  list(
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
}


check_step_thresholds <- function(p_enter, p_remove) {
  check_probability_open(p_enter, "p_enter")
  check_probability_open(p_remove, "p_remove")
  list(
    p_enter = p_enter,
    p_remove = p_remove
  )
}


check_step_count <- function(steps) {
  check_nonnegative_integerish(steps, "steps")
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
    stop(
      sprintf("'%s' must be a finite number", name),
      call. = FALSE
    )
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
    stop(sprintf("'%s' must be a character string", name), call. = FALSE)
  }
  if (!(x %in% choices)) {
    stop(
      sprintf("'%s' must be one of: %s", name, paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(x)
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
