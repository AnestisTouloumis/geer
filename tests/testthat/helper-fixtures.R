get_cached_test_fixtures <- function(topic) {
  builders <- c(
    test_native_fixtures, test_estimation_fixtures,
    test_methods_fixtures, test_diagnostics_fixtures,
    test_models_fixtures
  )
  if (!is.character(topic) || length(topic) != 1L || is.na(topic) ||
        !(topic %in% names(builders))) {
    stop("'topic' must name an existing test fixture", call. = FALSE)
  }
  if (is.null(test_data[[topic]])) {
    test_data[[topic]] <<- builders[[topic]]()
  }
  test_data[[topic]]
}
