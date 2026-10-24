
normalize_pvalue_labels <- function(lines) {
  matches <- gregexpr(" *< ?1e-04", lines)
  replacement <- lapply(regmatches(lines, matches), function(x) {
    vapply(
      x,
      function(label) formatC("<0.0001", width = max(nchar(label), 7L)),
      character(1),
      USE.NAMES = FALSE
    )
  })
  regmatches(lines, matches) <- replacement
  lines
}


print_with_small_pvalues <- function(print_function) {
  lines <- utils::capture.output(print_function())
  cat(normalize_pvalue_labels(lines), sep = "\n")
  invisible(NULL)
}


check_print_dots <- function(dots, context, functions) {
  for (print_function in functions) {
    print_call <- as.call(c(list(as.name("print")), list(x = NULL), dots))
    matched_call <- match.call(
      definition = print_function, call = print_call, expand.dots = FALSE
    )
    dots <- as.list(matched_call$...)
    if (!length(dots)) {
      break
    }
  }
  check_unused_dots(dots, context)
}


#' @title
#' Print Methods for geer Test Objects
#'
#' @description
#' Prints test tables and hypothesis test results returned by geer, displaying
#' p-values below 0.0001 as \code{<0.0001}.
#'
#' @param x an object of class \code{"geer_anova"} or
#'   \code{"geer_htest"}.
#' @param ... supported formatting arguments passed to the corresponding
#'   standard print method. Unused arguments are an error.
#'
#' @details
#' Test tables accept arguments of \code{stats::printCoefmat()} and the
#' standard ANOVA and default print methods. Hypothesis test results accept
#' \code{digits}, \code{prefix}, and default print arguments used for their
#' estimates. Unnamed and partially named formatting arguments retain standard
#' R argument matching.
#'
#' @return
#' The input object \code{x} is returned invisibly.
#'
#' @seealso \code{\link{anova.geer}}, \code{\link{runs_test}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment, family = poisson(), data = epilepsy,
#'   id = id
#' )
#' print(drop1(fit), digits = 3)
#'
#' @name print.geer_htest
#' @export
print.geer_anova <- function(x, ...) {
  args <- list(...)
  check_print_dots(
    args, "print.geer_anova",
    list(utils::getS3method("print", "anova"), stats::printCoefmat,
         base::print.default)
  )
  if (is.null(args$eps.Pvalue)) {
    args$eps.Pvalue <- 1e-04
  }
  y <- x
  class(y) <- setdiff(class(y), "geer_anova")
  print_with_small_pvalues(function() do.call(print, c(list(y), args)))
  invisible(x)
}


#' @rdname print.geer_htest
#' @export
print.geer_htest <- function(x, ...) {
  check_print_dots(
    list(...), "print.geer_htest",
    list(utils::getS3method("print", "htest"), base::print.default)
  )
  y <- x
  class(y) <- setdiff(class(y), "geer_htest")
  lines <- utils::capture.output(print(y, ...))
  p_value <- x$p.value
  if (length(p_value) == 1L && is.finite(p_value) && p_value < 1e-04) {
    text <- paste(lines, collapse = "\n")
    text <- sub(
      "p-value\\s*(=|<)\\s*[-+0-9.eE]+",
      "p-value < 0.0001",
      text
    )
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  }
  cat(lines, sep = "\n")
  invisible(x)
}
