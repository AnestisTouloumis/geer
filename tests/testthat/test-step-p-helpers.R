testthat::local_edition(3)

test_that("backward selector prioritizes zero-df drops", {
  aod <- data.frame(
    Df = c(NA, 1, 0, 1),
    Chi = c(NA, 2, 0, 3),
    `Pr(>Chi)` = c(NA, 0.40, NA, 0.90),
    check.names = FALSE
  )
  rownames(aod) <- c("<none>", "- x1", "- x2", "- x3")

  expect_identical(.step_p_selected_row_backward(aod), 3L)
})

test_that("backward selector otherwise uses the first maximum p-value", {
  aod <- data.frame(
    Df = c(NA, 1, 1, 1),
    Chi = c(NA, 2, 3, 4),
    `Pr(>Chi)` = c(NA, 0.70, 0.70, 0.20),
    check.names = FALSE
  )
  rownames(aod) <- c("<none>", "- a", "- b", "- c")

  expect_identical(.step_p_selected_row_backward(aod), 2L)
})

test_that("backward stop helper stops only when the selected p-value is below threshold", {
  aod_stop <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.01),
    check.names = FALSE
  )
  rownames(aod_stop) <- c("<none>", "- x1")

  aod_continue <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.40),
    check.names = FALSE
  )
  rownames(aod_continue) <- c("<none>", "- x1")

  row_stop <- .step_p_selected_row_backward(aod_stop)
  row_continue <- .step_p_selected_row_backward(aod_continue)

  expect_true(.step_p_backward_should_stop(aod_stop, row_stop, pvalue = 0.15))
  expect_false(.step_p_backward_should_stop(aod_continue, row_continue, pvalue = 0.15))
})

test_that("forward selector uses the first minimum p-value", {
  aod <- data.frame(
    Df = c(NA, 1, 1, 1),
    Chi = c(NA, 4, 5, 6),
    `Pr(>Chi)` = c(NA, 0.03, 0.03, 0.20),
    check.names = FALSE
  )
  rownames(aod) <- c("<none>", "+ x1", "+ x2", "+ x3")

  expect_identical(.step_p_selected_row_forward(aod), 2L)
})

test_that("forward stop helper stops only when the selected p-value is above threshold", {
  aod_stop <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.30),
    check.names = FALSE
  )
  rownames(aod_stop) <- c("<none>", "+ x1")

  aod_continue <- data.frame(
    Df = c(NA, 1),
    Chi = c(NA, 2),
    `Pr(>Chi)` = c(NA, 0.01),
    check.names = FALSE
  )
  rownames(aod_continue) <- c("<none>", "+ x1")

  row_stop <- .step_p_selected_row_forward(aod_stop)
  row_continue <- .step_p_selected_row_forward(aod_continue)

  expect_true(.step_p_forward_should_stop(aod_stop, row_stop, pvalue = 0.15))
  expect_false(.step_p_forward_should_stop(aod_continue, row_continue, pvalue = 0.15))
})

test_that("step label and term helpers return the expected values", {
  aod <- data.frame(Df = 1, row.names = "+ x1")

  expect_identical(.step_p_selected_step_label(aod, 1L), "+ x1")
  expect_identical(.step_p_selected_term("+ x1"), "x1")
  expect_identical(.step_p_selected_term("- x2"), "x2")
})

test_that("cycle detection depends on the underlying term", {
  expect_true(.step_p_is_cycle("+ x1", "- x1"))
  expect_true(.step_p_is_cycle("- x1", "+ x1"))
  expect_false(.step_p_is_cycle("+ x1", "- x2"))
  expect_false(.step_p_is_cycle("+ x1", NA_character_))
})
