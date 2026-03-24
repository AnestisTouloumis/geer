#' @export
geer <- function(x, ...) {
  UseMethod("geer")
}

#' @export
geer.default <- function(x, ...) {
  x <- new_geer(x)
  validate_geer(x)
}
