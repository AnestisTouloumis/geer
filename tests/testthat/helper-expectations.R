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
  expect_equal(
    object,
    published,
    tolerance = 0.5 * 10^(-digits) / mean(abs(published))
  )
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
  dev <- requireNamespace("pkgload", quietly = TRUE) &&
    isTRUE(pkgload::is_dev_package("geer"))
  if (dev) {
    out <- capture.output(
      eval(parse(text = code), envir = new.env(parent = asNamespace("geer")))
    )
    return(list(status = 0L, output = paste(out, collapse = "\n")))
  }
  old_libs <- Sys.getenv("R_LIBS", unset = NA)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  restore_libs <- function() {
    if (is.na(old_libs)) {
      Sys.unsetenv("R_LIBS")
    } else {
      Sys.setenv(R_LIBS = old_libs)
    }
  }
  on.exit(restore_libs(), add = TRUE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(paste0("suppressMessages(library(geer)); ", code))),
    stdout = TRUE,
    stderr = TRUE
  ))
  list(
    status = attr(out, "status") %||% 0L,
    output = paste(out, collapse = "\n")
  )
}
