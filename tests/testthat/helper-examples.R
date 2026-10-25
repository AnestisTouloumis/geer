compute_example_output <- function(topic) {
  if (!is.character(topic) || length(topic) != 1L || is.na(topic) ||
        !grepl("^[A-Za-z0-9._-]+$", topic)) {
    stop("'topic' must be a help topic name", call. = FALSE)
  }
  filename <- paste0(topic, ".Rd")
  path <- testthat::test_path("..", "..", "man", filename)
  rd <- if (file.exists(path)) {
    tools::parse_Rd(path)
  } else {
    tools::Rd_db("geer")[[filename]]
  }
  if (is.null(rd)) {
    stop("compute_example_output failed: help topic is missing", call. = FALSE)
  }
  path <- tempfile(fileext = ".R")
  on.exit(unlink(path), add = TRUE)
  tools::Rd2ex(rd, out = path, commentDonttest = TRUE)
  environment <- new.env(parent = asNamespace("geer"))
  capture.output({
    for (expression in parse(path)) {
      value <- withVisible(eval(expression, envir = environment))
      if (value$visible) {
        print(value$value)
      }
    }
  })
}
