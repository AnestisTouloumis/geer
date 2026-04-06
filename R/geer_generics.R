geer <- function(x, ...) {
  UseMethod("geer")
}

#' @exportS3Method
geer.default <- function(x, ...) {
  x <- new_geer(x)
  validate_geer(x)
}
