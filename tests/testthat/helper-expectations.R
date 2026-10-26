expect_test_result <- function(x) {
  values <- unlist(x[c("test_stat", "test_df", "test_p")])
  expect_length(values, 3L)
  expect_true(all(is.finite(values)))
  expect_true(x$test_stat >= 0 && x$test_df > 0)
  expect_true(x$test_p >= 0 && x$test_p <= 1)
}

jackknife_delete_estimates_by_refit <- function(fit, data, refit) {
  ids <- unique(data$id)
  out <- t(vapply(
    ids,
    function(id_out) coef(refit(data[data$id != id_out, ])),
    numeric(length(coef(fit)))
  ))
  dimnames(out) <- list(as.character(ids), names(coef(fit)))
  out
}

expect_jackknife_vcov <- function(fit, delete_estimates) {
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  expected <- ((nrow(delete_estimates) - 1L) / nrow(delete_estimates)) *
    crossprod(centered)
  dimnames(expected) <- list(names(coef(fit)), names(coef(fit)))
  observed <- vcov(fit, cov_type = "jackknife")
  expect_equal(compute_jackknife_deletes(fit), delete_estimates,
    tolerance = 1e-7
  )
  expect_equal(observed, expected, tolerance = 1e-7)
  expect_equal(observed, t(observed), tolerance = 1e-12)
}

expect_published <- function(object, published, digits) {
  if (!is.numeric(object) || !length(object) || any(!is.finite(object))) {
    stop("'object' must contain finite numeric values", call. = FALSE)
  }
  if (!is.numeric(published) || any(!is.finite(published))) {
    stop("'published' must contain finite numeric values", call. = FALSE)
  }
  if (length(object) != length(published) ||
        !identical(dim(object), dim(published)) ||
        !identical(names(object), names(published)) ||
        !identical(dimnames(object), dimnames(published))) {
    stop("'object' and 'published' must have matching shapes", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L ||
        !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("'digits' must be a nonnegative integer", call. = FALSE)
  }
  tolerance <- 0.5 * 10^(-digits)
  if (tolerance <= 0) {
    stop("'digits' must give a positive rounding unit", call. = FALSE)
  }
  for (index in seq_along(object)) {
    expect_equal(as.numeric(object[index] - published[index]), 0,
      tolerance = tolerance, info = paste("published coordinate", index)
    )
  }
}


expect_published_table <- function(fit, estimate, std_error, p_value) {
  coefs <- summary(fit)$coefficients
  shown <- !is.na(p_value)
  expect_published(unname(coefs[, "Estimate"]), estimate, 4)
  expect_published(unname(coefs[, "Std. Error"]), std_error, 4)
  expect_published(unname(coefs[shown, "Pr(>|z|)"]), p_value[shown], 4)
}

collect_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = w)
}

in_subprocess <- function(code) {
  if (!is.character(code) || length(code) != 1L || is.na(code)) {
    stop("'code' must be a character string", call. = FALSE)
  }
  parse(text = code)
  dev <- requireNamespace("pkgload", quietly = TRUE) &&
    isTRUE(pkgload::is_dev_package("geer"))
  startup <- if (dev) {
    path <- getNamespaceInfo(asNamespace("geer"), "path")
    paste0(
      "suppressMessages(pkgload::load_all(", deparse(path),
      ", quiet = TRUE, compile = FALSE, helpers = FALSE))"
    )
  } else {
    "suppressMessages(library(geer))"
  }
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(startup, code), script)
  withr::local_envvar(
    R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)
  )
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(script), stdout = TRUE, stderr = TRUE
  ))
  list(
    status = attr(out, "status") %||% 0L,
    output = paste(out, collapse = "\n")
  )
}
