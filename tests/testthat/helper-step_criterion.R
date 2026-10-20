step_criterion_main_terms <- c(
  "baseline", "treatment", "gender", "visit", "age"
)

criterion_value <- function(fit, criterion, cov_type = "robust") {
  out <- geecriteria(
    fit,
    criteria = criterion,
    cov_type = cov_type,
    digits = 15
  )
  out[[criterion]]
}

expect_greedy_path <- function(out,
                               start,
                               criterion,
                               direction,
                               lower = character(),
                               upper = step_criterion_main_terms,
                               cov_type = "robust") {
  current <- start
  current_value <- criterion_value(current, criterion, cov_type)
  aod <- out$anova
  expect_equal(aod[[criterion]][1L], current_value, tolerance = 1e-10)
  for (i in seq_len(nrow(aod))) {
    terms <- attr(terms(formula(current)), "term.labels")
    drops <- if (direction != "forward") setdiff(terms, lower) else character()
    adds <- if (direction != "backward") setdiff(upper, terms) else character()
    changes <- c(sprintf("- %s", drops), sprintf("+ %s", adds))
    fits <- lapply(changes, function(change) {
      call <- update(
        current,
        as.formula(paste(". ~ .", change)),
        evaluate = FALSE
      )
      eval(call, environment(formula(start)))
    })
    values <- vapply(
      fits,
      criterion_value,
      numeric(1),
      criterion = criterion,
      cov_type = cov_type
    )
    best <- which.min(values)
    if (i == nrow(aod)) {
      expect_true(length(values) == 0L || values[best] >= current_value)
    } else {
      expect_identical(
        rownames(aod)[i + 1L],
        paste0("Step ", i, ": ", changes[best])
      )
      expect_equal(aod[[criterion]][i + 1L], values[best], tolerance = 1e-10)
      current <- fits[[best]]
      current_value <- values[best]
    }
  }
  expect_setequal(
    attr(terms(formula(out)), "term.labels"),
    attr(terms(formula(current)), "term.labels")
  )
}
