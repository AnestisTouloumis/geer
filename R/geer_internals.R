format_perc <- function(probs, digits = 2, check = FALSE) {
  if (!is.numeric(probs)) {
    stop("'probs' must be numeric", call. = FALSE)
  }
  if (!is_pos_int_scalar(digits)) {
    stop("'digits' must be a positive integer", call. = FALSE)
  }
  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("'check' must be a single logical value", call. = FALSE)
  }
  if (check && any(probs < 0 | probs > 1, na.rm = TRUE)) {
    warning("Some probabilities are outside [0, 1]", call. = FALSE)
  }
  formatted <- formatC(
    100 * probs,
    digits = digits,
    format = "fg",
    flag = "",
    drop0trailing = FALSE
  )
  paste0(trimws(formatted), "%")
}


is_pos_scalar <- function(x) {
  is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x > 0
}

is_pos_int_scalar <- function(x, tol = .Machine$double.eps^0.5) {
  is_pos_scalar(x) && abs(x - round(x)) < tol
}

check_single_numeric <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric value", name), call. = FALSE)
  }
}

check_positive_numeric <- function(x, name) {
  check_single_numeric(x, name)
  if (x <= 0) {
    stop(sprintf("'%s' must be positive", name), call. = FALSE)
  }
}

check_probability <- function(x, name, open = TRUE) {
  check_single_numeric(x, name)
  if (open) {
    if (x <= 0 || x >= 1) {
      stop(sprintf("'%s' must be strictly between 0 and 1", name), call. = FALSE)
    }
  } else {
    if (x < 0 || x > 1) {
      stop(sprintf("'%s' must be between 0 and 1", name), call. = FALSE)
    }
  }
}

check_choice <- function(x, choices, name) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !(x %in% choices)) {
    stop(
      sprintf("'%s' should be one of: %s", name, paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }
}

test_label <- function(test) {
  switch(
    test,
    wald = "Wald",
    score = "Score",
    `working-wald` = "Modified Working Wald",
    `working-score` = "Modified Working Score",
    `working-lrt` = "Modified Working LRT"
  )
}
