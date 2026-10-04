## Display helpers: p-values smaller than 0.0001 are printed as "<0.0001".
## printCoefmat()/format.pval() can only produce "<1e-04" (or "< 1e-04"), so
## the printed text is post-processed. Replacement keeps the column width, so
## the table stays aligned.

fix_small_pvalue_labels <- function(lines) {
  matches <- gregexpr(" *< ?1e-04", lines)
  replacement <- lapply(regmatches(lines, matches), function(x) {
    vapply(
      x,
      function(label) formatC("<0.0001", width = max(nchar(label), 7L)),
      character(1L),
      USE.NAMES = FALSE
    )
  })
  regmatches(lines, matches) <- replacement
  lines
}


print_with_small_pvalues <- function(print_function) {
  lines <- utils::capture.output(print_function())
  cat(fix_small_pvalue_labels(lines), sep = "\n")
  invisible(NULL)
}


#' @export
print.geer_anova <- function(x, ...) {
  args <- list(...)
  if (is.null(args$eps.Pvalue)) {
    args$eps.Pvalue <- 1e-04
  }
  y <- x
  class(y) <- setdiff(class(y), "geer_anova")
  print_with_small_pvalues(function() do.call(print, c(list(y), args)))
  invisible(x)
}


## Test results are "htest" objects. print.htest() writes the p-value as
## "p-value = 3e-08" or "p-value < 2.2e-16" and offers no way to change that, so
## the printed text is post-processed: any p-value below 0.0001 is shown as
## "p-value < 0.0001". Stored p-values are untouched.
#' @export
print.geer_htest <- function(x, ...) {
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
