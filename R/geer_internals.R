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
